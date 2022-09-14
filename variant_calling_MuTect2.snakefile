import os

configfile: "TrimmedFastqFilePaths_local.rdgroups.withnormals.json"

gatk_path = "gatk"

refpath_XX ="T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
#refpath_XY = "T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",
refpath_XY = "/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",


rule all:
   input: # variant calling and filtering with MuTect2
        expand("variant_calling_mutect2/{sample}.somatic.vcf.gz", sample=config["test_tumors"]),
        expand("variant_calling_mutect2/{sample}.somatic.filtered.vcf.gz", sample=config["test_tumors"]),
        expand("variant_calling_mutect2/{sample}.somatic.filtered.pass.vcf.gz", sample=config["test_tumors"])

#ref ="T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",

rule tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam")
    output:
        os.path.join("variant_calling_mutect2/", "{sample}.somatic.vcf.gz")
    params:
        gatk = gatk_path,
        sm = lambda wildcards: config[wildcards.sample]["normal"]
    shell:
        """
        {params.gatk} Mutect2 -R {input.ref} -I {input.tumor_bam} -normal {params.sm} -O {output}
        """

rule filter_XY:
    input:
        ref = refpath_XY,
        unfiltered = os.path.join("variant_calling_mutect2/", "{sample}.somatic.vcf.gz")
    output:
        filtered = os.path.join("variant_calling_mutect2/", "{sample}.somatic.filtered.vcf.gz")
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} FilterMutectCalls -R {input.ref} -V {input.unfiltered} -O {output.filtered}
        """

rule select_pass_variants_XY:
    input:
        ref = refpath_XY,
        vcf = os.path.join("variant_calling_mutect2/", "{sample}.somatic.filtered.vcf.gz")
    output:
        os.path.join("variant_calling_mutect2/", "{sample}.somatic.filtered.pass.vcf.gz")
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.vcf} --exclude-filtered -O {output}
        """
