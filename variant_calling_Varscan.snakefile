import os

configfile: "/scratch/splaisie/gbm_imaging/test_mutect2/TrimmedFastqFilePaths_local.rdgroups.withnormals.json"

varscan_path = "/home/splaisie/software/Varscan/VarScan.v2.3.9.jar"
gatk_path = "gatk"

alignment_directory = "/scratch/splaisie/gbm_imaging/test_mutect2/"  # this is where your alignments/ should be

refpath_XX ="T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
#refpath_XY = "T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",
refpath_XY = "/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa"


rule all:
   input: # variant calling and filtering with varscan
        expand("variant_calling_varscan/pileups/{sample}.T2Tv2.XY.pileup", sample=config["test_tumors"]),
        #expand("variant_calling_varscan/pileups/{sample}.pileup", sample=config["all_germline_samples"]),
        expand("variant_calling_varscan/pileups/{sample}.T2Tv2.XY.pileup", sample="G92BL_CKDN200007839-1A_trimmed"),
        expand("variant_calling_varscan/{sample}.varscan.snp", sample=config["test_tumors"]),
        expand("variant_calling_varscan/{sample}.varscan.snp.Somatic.hc", sample=config["test_tumors"]),
        expand("variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter", sample=config["test_tumors"]),
        expand("variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter.vcf", sample=config["test_tumors"]),
        expand("variant_calling_varscan/{sample}.varscan.indel", sample=config["test_tumors"])

rule bam_pileup: #for both normal and tumor
    input:
        ref = refpath_XY,
        bam = os.path.join(alignment_directory,"alignments/", "{sample}.T2Tv2.XY.dedup.bam")
    output:
        pileup = "variant_calling_varscan/pileups/{sample}.T2Tv2.XY.pileup"
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output.pileup}
        """

rule run_varscan:
    input:
        ref = refpath_XY,
        normal_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.pileup"),
        tumor_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.pileup"),
    output:
        snp = "variant_calling_varscan/{sample}.varscan.snp",
        indel = "variant_calling_varscan/{sample}.varscan.indel"
    params:
        varscan = varscan_path,
        basename = "variant_calling_varscan/{sample}.varscan"
    shell:
        """
        java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05
        """

rule isolate_calls_by_type_and_confidence:
    input: 
        snp = "variant_calling_varscan/{sample}.varscan.snp",
    output: 
        snp_somatic_hc = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc",
    params: 
        varscan = varscan_path
    shell:
        """
        java -jar {params.varscan} processSomatic {input.snp}
        """

rule somatic_filter:
    input: 
        snp_somatic_hc = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc",
        indel = "variant_calling_varscan/{sample}.varscan.indel",
    output: 
        snp_somatic_hc_filter = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter",
    params: 
        varscan = varscan_path
    shell:
        """
        java -jar {params.varscan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_vcf:
    input:
        snp_somatic_hc_filter = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter",
    output:
        snp_somatic_hc_filter_vcf = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter.vcf",
    shell:
        """
        cat {input.snp_somatic_hc_filter} | awk '{{print $1"\t" $2"\t" "." "\t" $3 "\t" $4}}' > {output.snp_somatic_hc_filter_vcf}
        """

