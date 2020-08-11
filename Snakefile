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
                "merge.txt"

rule q20_filter_bam:
        input:  
                "dedup_bam/{sample}_dedup.bam"
        output: 
                "q20_filter/{sample}_q20.bam"
        shell:  
                "samtools view -q 20 -b -S {input} > {output}"

rule sync_file:
        input:
                "q20_filter/{sample}_q20.bam"
        output:
                mpileup: "pileup/BSAF7.mpileup"
                sync_file: "sync/BSAF7.sync" 
        shell:
                """
                samtools mpileup -B {input} > {output.mpileup};
                java -ea -Xmx7g -jar /workdir/jcf236/popoolation2_1201/mpileup2sync.jar \
                --input {output.mpileup} --output {output.sync_file} --fastq-type anger --min-qual 20 --threads 24
                """        
                
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
