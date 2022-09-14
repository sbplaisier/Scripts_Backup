import os

configfile: "/scratch/splaisie/gbm_imaging/test_mutect2/TrimmedFastqFilePaths_local.rdgroups.withnormals.json"

strelka_path = "/home/splaisie/software/Strelka/strelka-2.9.10.centos6_x86_64/bin"
gatk_path = "gatk"

alignment_directory = "/scratch/splaisie/gbm_imaging/test_mutect2/"  # this is where your alignments/ should be

refpath_XX ="T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
#refpath_XY = "T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",
refpath_XY = "/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",

#ruleorder: configure_tumor_with_matched_normal_XY > run_tumor_with_matched_normal_XY 

rule all:
   input: # variant calling and filtering with strelka
        expand("variant_calling_strelka/{sample}/runWorkflow.py", sample=config["test_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.snvs.vcf.gz", sample=config["test_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.snvs.pass.vcf.gz", sample=config["test_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.indels.vcf.gz", sample=config["test_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.indels.pass.vcf.gz", sample=config["test_tumors"])


rule configure_tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join(alignment_directory,"alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join(alignment_directory,"alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam")
    output:
        "variant_calling_strelka/{sample}/runWorkflow.py"
    params:
        strelka = strelka_path,
        run_dir = "variant_calling_strelka/{sample}",
    shell:
        """
        {params.strelka}/configureStrelkaSomaticWorkflow.py --referenceFasta {input.ref} --tumorBam {input.tumor_bam} --normalBam {input.normal_bam} --runDir {params.run_dir} --exome
        """

rule run_tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join(alignment_directory,"alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join(alignment_directory, "alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam"),
        workflow = "variant_calling_strelka/{sample}/runWorkflow.py"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.indels.vcf.gz"
    params:
        run = "variant_calling_strelka/{sample}/runWorkflow.py",
    shell:
        """
        {params.run} -m local -j 20
        """

rule select_pass_variants_XY:
    input:
        ref = refpath_XY,
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.indels.vcf.gz"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.snvs.pass.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.indels.pass.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.snvs} --exclude-filtered -O {output.snvs}
        {params.gatk} SelectVariants -R {input.ref} -V {input.indels} --exclude-filtered -O {output.indels}
        """
