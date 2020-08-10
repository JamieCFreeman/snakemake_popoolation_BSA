

rule q20_filter_bam:
    input:
        "dedup_bam/BSAF7_1_RES_abcd_dedup.bam"
    output:
        "q20_filter/BSAF7_1_RES_q20.bam"
    shell:
    """
    samtools view -q 20 -b -S {input} > {output}
    """

rule pileup:
    input:
        "q20_filter/BSAF7_1_RES_q20.bam"
    output:
        "pileup/BSAF7_1_RES.pileup"
