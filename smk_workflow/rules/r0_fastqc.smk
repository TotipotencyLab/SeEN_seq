
# Note on zcat: 
# Somehow zcat always add Z at the end of the file, so using gunzip --stdout is a better option
# The wildcards functions for retrieving input fastq file are defined in the Snakefile

# MARK: Single-end
rule r0_fastqc_raw_SE:
    input: r0_wc_get_fastqc_raw_SE
    output:
        html = f"{proj_dir}/0_fastqc_raw/{{sample}}_fastqc.html",
        zip = temp(f"{proj_dir}/0_fastqc_raw/{{sample}}_fastqc.zip"),
    log: f"{proj_dir}/0_fastqc_raw/log/{{sample}}.log"
    params:
        outdir = f"{proj_dir}/0_fastqc_raw",
        extra = "--quiet",
    shell:
        """
        [[ "{input}" == *.gz ]] && cat_cmd="gunzip --stdout" || cat_cmd="cat"  # If end with gz use zcat, otherwise use cat
        cmd="$cat_cmd {input} | fastqc stdin:{wildcards.sample} --threads 1 {params.extra} -o {params.outdir}"
        echo $cmd > {log}; eval $cmd
        """


# MARK: Paired-end
rule r0_fastqc_raw_PE:
    input: 
        # Only run this rule if both R1 and R2 are present
        R1 = r0_wc_get_fastqc_raw_R1,
        R2 = r0_wc_get_fastqc_raw_R2,
    output:
        R1_html = f"{proj_dir}/0_fastqc_raw/{{sample}}_R1_fastqc.html",
        R2_html = f"{proj_dir}/0_fastqc_raw/{{sample}}_R2_fastqc.html",
        R1_zip = temp(f"{proj_dir}/0_fastqc_raw/{{sample}}_R1_fastqc.zip"),
        R2_zip = temp(f"{proj_dir}/0_fastqc_raw/{{sample}}_R2_fastqc.zip"),
    log: f"{proj_dir}/0_fastqc_raw/log/{{sample}}.log"
    params:
        outdir = f"{proj_dir}/0_fastqc_raw",
        extra = "--quiet",
    shell:
        """
        [[ {input.R1} == *.gz ]] && cat_cmd="gunzip --stdout" || cat_cmd="cat"
        cmd="$cat_cmd {input.R1} | fastqc stdin:{wildcards.sample}_R1 --threads 1 {params.extra} -o {params.outdir}"
        echo $cmd > {log}; eval $cmd >> {log} 2>&1

        [[ {input.R2} == *.gz ]] && cat_cmd="gunzip --stdout" || cat_cmd="cat"
        cmd="$cat_cmd {input.R2} | fastqc stdin:{wildcards.sample}_R2 --threads 1 {params.extra} -o {params.outdir}"
        printf '\n\n' >> {log}; echo $cmd >> {log}; eval $cmd >> {log} 2>&1
        """

ruleorder: r0_fastqc_raw_PE > r0_fastqc_raw_SE
