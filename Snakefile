import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####
configfile: "config.yaml"

samples_table = pd.read_csv(config["samples"], sep='\t')


##### rules #####

wildcard_constraints:
        sample="|".join(samples_table['samples'])

rule all:
        input:
                "AFD/BSAF7_pwc",
                "stats/BSAF7_w500.fst",
                expand("stats/BSAF7_{min}.gwas", min = config["CMH_min"])

rule q20_filter_bam:
        input:
                "dedup_bam/{sample}_dedup.bam"
        output:
                "q20_filter/{sample}_q20.bam"
        shell:
                "samtools view -q 20 -b -S {input} > {output}"

rule mpileup:
        input:
                expand("q20_filter/{sample}_q20.bam", sample=samples_table["samples"])
        output:
                "pileup/BSAF7.mpileup"
        shell:
                """
                samtools mpileup -B {input} > {output}
                """
rule sync_file:
        input:
                "pileup/BSAF7.mpileup"
        output:
                "sync/BSAF7.sync"
        shell:
                """
                java -ea -Xmx7g -jar /workdir/jcf236/popoolation2_1201/mpileup2sync.jar \
                  --input {input} --output {output}.sync --fastq-type sanger --min-qual 20 --threads 24
                """

rule AF_diff:
        input:
                "sync/BSAF7.sync"
        output:
                "AFD/BSAF7_pwc"
        shell:
                """
                out_pref=`echo {input} | sed 's/.sync//' | sed 's ^.*/  '`;
                perl /workdir/jcf236/popoolation2_1201/snp-frequency-diff.pl --input {input} --output-prefix ./AFD/${{out_pref}} \
                        --min-count 6 --min-coverage 20 --max-coverage 200
                """

rule Fst:
        input:
                "sync/BSAF7.sync"
        output:
                "stats/BSAF7_w500.fst"
        shell:
                """
                perl /workdir/jcf236/popoolation2_1201/fst-sliding.pl --input {input} \
                --output {output} --min-count 6 --min-coverage 10 --max-coverage 200 \
                --min-covered-fraction 1 --window-size 5000 --step-size 5000 --pool-size 100 \
                """

rule CMH_test:
        input:
                "sync/BSAF7.sync"
        output:
                expand("stats/BSAF7_{min}.cmh", min = config["CMH_min"])
        params:
                min_cov    = config["CMH_min"],
                max_cov    = 100,
                min_count  = 6
        shell:
                """
                perl /workdir/jcf236/popoolation2_1201/cmh-test.pl --input {input} --output {output} --min-count {params.min_count}  \
                --min-coverage {params.min_cov} --max-coverage {params.max_cov} --population 1-2,4-5,7-8
                """

rule CMH_gwas:
        input:
                "stats/BSAF7_{CMH_min}.cmh"
        output:
                "stats/BSAF7_{CMH_min}.gwas"
        params:
                min_pvalue  = 1.0e-20
        shell:
                "perl /workdir/jcf236/popoolation2_1201/export/cmh2gwas.pl --input {input} --output {output} --min-pvalue {params.min_pvalue}"

rule merge:
        input:
                expand("q20_filter/{sample}_q20.bam", sample=samples_table["samples"])
        output:
                "merge.txt"
        shell:
                "echo hello > {output}"

rule clean:
        shell:
                "rm merge.txt"
                                        
