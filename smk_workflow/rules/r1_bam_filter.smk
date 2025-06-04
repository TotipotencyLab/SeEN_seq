# /////////////////////////////////////////////////////////////////////////////////////////////////
#region MARK: BAM processing
# /////////////////////////////////////////////////////////////////////////////////////////////////

bt2_sort_primary_map = False
if bt2_sort_primary_map:
    # Sort the primary mapped reads before converting to BAM
    rule samtools_sam_to_bam:
        input:
            sam = f"{proj_dir}/1_bam/{{sample}}.sam",
        output:
            bam = f"{proj_dir}/1_bam/{{sample}}.bam",
        log:
            f"{proj_dir}/1_bam/log/{{sample}}_sam_to_bam.log"
        threads: 4
        shell:
            """
            samtools sort -@ {threads} {input.sam} | samtools view -o {output.bam} -bS - 2> {log}
            """
else:
    # Convert SAM to BAM without sorting
    rule samtools_sam_to_bam:
        input:
            sam = f"{proj_dir}/1_bam/{{sample}}.sam"
        output:
            bam = f"{proj_dir}/1_bam/{{sample}}.bam"
        log:
            f"{proj_dir}/1_bam/log/{{sample}}_sam_to_bam.log"
        shell:
            """
            samtools view -o {output.bam} -bS {input.sam} 2> {log}
            """


# Filter BAM based on the flag and quality, then sort
# TODO: Test this rule on linux
rule samtools_filter_sort:
    input:
        bam = f"{proj_dir}/1_bam/{{sample}}.bam"
    output:
        bam = f"{proj_dir}/1_bam/{{sample}}_sort.bam"
    log: f"{proj_dir}/1_bam/log/{{sample}}_sort.log"
    params:
        flags = "-F 0x04",  # include only mapped reads
        quality = "-q 30",
        extra = "",
    threads: 4
    shell:
        """
        samtools view -b {params.flags} {params.quality} {params.extra} {input.bam} | samtools sort -@ {threads} -o {output} - > {log} 2>&1
        """
# samtools view -Sb -F {params.exclude_flag} -q {params.min_mapq} {input} | samtools sort -@ {threads} -o {output} 2> {log.err} 1> {log.out}


# macOS specific rule
# NOTE: For mac users:
# It seems that samtools on mac os doesn't support thread option (-@) 
#    and the piping of samtools view to samtools sort command doesn't work the way I used to do it on linux platform
# This pipeline is developed on
# system: Mac OS
# samtools version: 0.1.18 (r982:295)

rule samtools_filter_sort_mac:
    input:
        bam = f"{proj_dir}/1_bam/{{sample}}.bam"
    output:
        bam = f"{proj_dir}/1_bam/{{sample}}_sort.bam"
    log: f"{proj_dir}/1_bam/log/{{sample}}_sort.log"
    params:
        flags = "-F 0x04",  # include only mapped reads
        quality = "-q 30",
        extra = "",
    threads: 1
    shell:
        """
        out_prefix=$(echo {output.bam} | sed 's/\.bam$//')
        samtools view -b {params.flags} {params.quality} {params.extra} {input.bam} | samtools sort - $out_prefix > {log} 2>&1
        """

if(platform.system() == "Darwin"):
    ruleorder: samtools_filter_sort_mac > samtools_filter_sort
