# The wildcards functions for retrieving input fastq file are defined in the Snakefile

trim_params = config["trim_galore_params"]
trim_params_extra = f"--quality {trim_params['quality']} --length {trim_params['length']} {trim_params['extra']}"

rule r0_trim_galore_SE:
    input: r0_wc_get_fastqc_raw_SE
    output:
        trimmed = f"{proj_dir}/0_fastq_trimmed/{{sample}}_trimmed.fq.gz",
    log: f"{proj_dir}/0_fastq_trimmed/log/{{sample}}.log"
    params:
        outdir = f"{proj_dir}/0_fastq_trimmed",
        extra = trim_params_extra,
    threads: 4
    shell:
        """
        trim_galore --cores {threads} {params.extra} --basename {wildcards.sample} -o {params.outdir} {input} > {log} 2>&1
        """

rule r0_trim_galore_PE:
    input: 
        # Only run this rule if both R1 and R2 are present
        R1 = r0_wc_get_fastqc_raw_R1,
        R2 = r0_wc_get_fastqc_raw_R2,
    output:
        R1_trimmed = f"{proj_dir}/0_fastq_trimmed/{{sample}}_val_1.fq.gz",
        R2_trimmed = f"{proj_dir}/0_fastq_trimmed/{{sample}}_val_2.fq.gz",
    log: f"{proj_dir}/0_fastq_trimmed/log/{{sample}}.log"
    params:
        outdir = f"{proj_dir}/0_fastq_trimmed",
        extra = trim_params_extra,
    threads: 4
    shell:
        """
        trim_galore --cores {threads} --paired {params.extra} --basename {wildcards.sample} -o {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """
