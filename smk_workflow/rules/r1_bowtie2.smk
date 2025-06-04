rule bowtie2_index_ref:
    input:
        fastq = config["ref_seq_path"]
    output:
        index_dir = directory(f"{proj_dir}/bt2_index")
    log:
        f"{proj_dir}/bt2_index/bowtie2_index.log"
    params:
        index_name = config["project_name"],
    threads: 1
    shell:
        """
        bowtie2-build --verbose --threads {threads} {input.fastq} {output.index_dir}/{params.index_name} > {log} 2>&1
        """
# /////////////////////////////////////////////////////////////////////////////////////////////////
#region MARK: mapping
# /////////////////////////////////////////////////////////////////////////////////////////////////

bt2_params = config["bowtie2"]["extra_params"]
# bt2_params = "-N 1 -L 20 -k 1 --no-unal --no-mixed -I 1 --no-discordant --dovetail --no-contain --no-overlap"
bt2_params = re.sub(r"\s+", " ", bt2_params).strip()
bt2_params_SE = bt2_params

# Processing the single-end specific parameters, if needed
# Here are the bowtie2 flags that are specific to paired-end reads:
# -I/--minins <int>  minimum fragment length (0)
# -X/--maxins <int>  maximum fragment length (500)
# --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
# --no-mixed         suppress unpaired alignments for paired reads
# --no-discordant    suppress discordant alignments for paired reads
# --dovetail         concordant when mates extend past each other
# --no-contain       not concordant when one mate alignment contains other
# --no-overlap       not concordant when mates overlap at all
if any(FASTQ_TABLE["layout"]=="SE"):
    bt2_params_PE_flag_regex = [r"\-I\s+\d+", r"\-\-minins\s+\d+", r"\-X\s+\d+", r"\-\-maxins\s+\d+" r"\-\-((fr)|(rf)|(ff))", r"\-\-no\-mixed", r"\-\-no\-discordant", r"\-\-dovetail", r"\-\-no\-contain", r"\-\-no\-overlap"]
    # Remove the paired-end specific flags if we found one
    for regex in bt2_params_PE_flag_regex:
        bt2_params_SE = re.sub(regex, "", bt2_params_SE)
    bt2_params_SE = re.sub(r"\s+", " ", bt2_params_SE).strip()
    if bt2_params_SE != bt2_params:
        print(f"INFO: The single-end samples will use these Bowtie2 parameters\n\t: {bt2_params_SE}\n")


rule bowtie2_SE:
    input: 
        index = rules.bowtie2_index_ref.output.index_dir,
        fastq = f"{proj_dir}/0_fastq_trimmed/{{sample}}_trimmed.fq.gz"
    output:
        sam = temp(f"{proj_dir}/1_bam/{{sample}}.sam")
    log:
        f"{proj_dir}/1_bam/log/{{sample}}.log"
    params:
        index = f"{proj_dir}/bt2_index/{config['project_name']}",
        extra = bt2_params_SE,
    threads: 4
    shell:
        """
        bowtie2 -x {params.index} -U {input.fq} -p {params.threads} > {output.sam} 2> {log}
        """

rule bowtie2_PE:
    input:
        index = rules.bowtie2_index_ref.output.index_dir,
        fastq1 = f"{proj_dir}/0_fastq_trimmed/{{sample}}_val_1.fq.gz",
        fastq2 = f"{proj_dir}/0_fastq_trimmed/{{sample}}_val_2.fq.gz",
    output:
        sam = temp(f"{proj_dir}/1_bam/{{sample}}.sam")
    log:
        f"{proj_dir}/1_bam/log/{{sample}}.log"
    params:
        index = f"{proj_dir}/bt2_index/{config['project_name']}",
        extra = bt2_params,
    threads: 4
    shell:
        """
        bowtie2 -p {threads} {params.extra} -x {params.index} -1 {input.fastq1} -2 {input.fastq2} > {output.sam} 2> {log}
        """

ruleorder: bowtie2_PE > bowtie2_SE
