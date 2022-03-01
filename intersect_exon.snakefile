# This pipeline produces converts vcfs for heterozygous biallelic snps
#    to bed format and then makes an exon filtered version

chromosomes = ["chrX","chr8"]
samples = ["MW-11","MW-15","MW-21","MW-31","MW-33","MW-43","MW-53","OBG0055-D5","OBG0055-P1"]
exonref = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.exon.bed"

rule all:
    input:
        expand("{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf", sample=samples, chr=chromosomes),
        expand("{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.exon.bed", sample=samples, chr=chromosomes)
 
rule het_variants_bed:
    input:
        variant_vcf = "{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf"
    output:
        variant_bed = "{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.bed"
    shell: 
        "vcf2bed < {input.variant_vcf} > {output.variant_bed}"

rule intersect_exon:
    input:
        variant_bed = "{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.bed"
    output:
        exon_filtered_bed = "{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.exon.bed"
    params:
        ref = exonref
    shell:
        "bedtools intersect -wa -u -a {input.variant_bed} -b {params.ref} > {output.exon_filtered_bed}"


