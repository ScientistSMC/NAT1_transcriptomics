###### Snakefile for fastqc & alignment to hg38 for human RNA-seq data ######

from os.path import join

configfile:"rna_samples.json"


#Load DIR
GENOME="GRCh38.primary_assembly.genome.fa"
REF_GENOME="ref/gencode_hg38/"+GENOME
FASTQ="fastq_files/"
FASTQC_DIR="fastqc/"
FASTQ_TRIM="fastq_files_trimmed/"
FASTQC_TRIM_DIR="fastqc_trimmed/"
GTF="gencode.v31.annotation.gtf"
REF_GTF="ref/gencode_hg38/"+GTF
STAR_REF= "ref/STARref_hg38/"
STAR_OUT="STAR_hg38_alignment/"

###RULES###

#rule clean:
#    shell:"rm -rf .snakemake"


rule all:
	input:
		expand(FASTQC_TRIM_DIR+"{sample}_ALL_R1_trimmed_fastqc.zip", sample=config["samples"]),
        expand(FASTQC_DIR+"{sample}_ALL_R1_fastqc.html", sample=config["samples"]),
        expand(FASTQ_TRIM+"{sample}_ALL_R1_trimmed.fq.gz", sample=config["samples"]),
        expand(STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam", sample=config["samples"]),
        expand(STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
        "multiqc_fastqc/multiqc_report.html",
        "multiqc_align/multiqc_report.html",
        "multiqc_fastqc_trimmed/multiqc_report.html"



###FASTQC###

rule fastqc:
    input:
        r1=FASTQ+"{sample}_ALL_R1.fastq.gz"
    output:
        html=FASTQC_DIR+"{sample}_ALL_R1_fastqc.html",
        zip=FASTQC_DIR+"{sample}_ALL_R1_fastqc.zip"
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'fastqc -f fastq {input.r1} --outdir fastqc '


#generating multiqc report on fastqc

rule multiqc_fastqc:
    input:
        fastqc=FASTQC_DIR,
        r1=expand(FASTQC_DIR+"{sample}_ALL_R1_fastqc.html", sample=config["samples"])
    output:
        "multiqc_fastqc/multiqc_report.html"
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'multiqc -f {input.fastqc} -o multiqc_fastqc'

###ADAPTER TRIMMING###
# Trimming adaptor sequences; val in output file stands for validated

rule fastq_trim:
    input:
        r1=FASTQ+"{sample}_ALL_R1.fastq.gz",
    output:
        r1=FASTQ_TRIM+"{sample}_ALL_R1_trimmed.fq.gz",
        r1_trim_report=FASTQ_TRIM+"{sample}_ALL_R1.fastq.gz_trimming_report.txt"
    conda:
        "envs/trim-galore_conda_env.yaml"
    shell:
        "trim_galore -o fastq_files_trimmed/ {input.r1} "

rule fastqc_trimmed_paired:
    input:
        FASTQ_TRIM+"{sample}_ALL_R1_trimmed.fq.gz"
    output:
        html=FASTQC_TRIM_DIR+"{sample}_ALL_R1_trimmed_fastqc.html",
        zip=FASTQC_TRIM_DIR+"{sample}_ALL_R1_trimmed_fastqc.zip",
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'fastqc -f fastq {input} --outdir fastqc_trimmed'

#generating multiqc report on trimmed fastqc

rule multiqc_fastqc_trimmed:
    input:
        fastqc=FASTQC_TRIM_DIR,
        html=expand(FASTQC_TRIM_DIR+"{sample}_ALL_R1_trimmed_fastqc.html", sample=config["samples"]),
        zip=expand(FASTQC_TRIM_DIR+"{sample}_ALL_R1_trimmed_fastqc.zip", sample=config["samples"])
    output:
        "multiqc_fastqc_trimmed/multiqc_report.html"
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'PYTHONPATH= && '
        'multiqc -f {input.fastqc} -o multiqc_fastqc_trimmed'


###DOWNLOAD REF GENOME & GTF (GENCODE V31)###

rule download_genome:
    output:
        fasta=REF_GENOME
    shell:
        'wget -O "{output.fasta}.gz" "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.primary_assembly.genome.fa.gz" && '
        'gunzip {output.fasta}".gz" '


rule download_gtf:
    output:
        gtf=REF_GTF
    shell:
        'wget -O "{output.gtf}.gz" "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz" && '
        'gunzip {output.gtf}".gz" '


rule download_genome_md5:
    output:
        md5="MD5SUMS.txt"
    shell:
        'wget -O "{output.md5}.gz" "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/MD5SUMS" && '


###GENERATE GENOME INDEX###

rule generate_genome_index:
    input:
        genome=REF_GENOME,
        gtf=REF_GTF
    output:
        STAR_REF+"SAindex"
    threads: 12
    params:
        star_genome_dir=STAR_REF
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.star_genome_dir} '
        '--genomeFastaFiles {input.genome} '
        '--sjdbGTFfile {input.gtf} '



###ALIGNMENT###

rule star_2pass_alignment:
    input:
        r1=FASTQ_TRIM+"{sample}_ALL_R1_trimmed.fq.gz",
        gtf=REF_GTF,
        indexed_ref= STAR_REF+"SAindex"
    output:
        STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam",
        STAR_OUT+"trimmed/{sample}.SJ.out.tab",
        STAR_OUT+"trimmed/{sample}.Log.out",
        STAR_OUT+"trimmed/{sample}.Log.final.out",
        STAR_OUT+"trimmed/{sample}.Log.progress.out",
    threads: 12
    params:
        genomeDir= STAR_REF
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genomeDir} '
        '--sjdbGTFfile {input.gtf} '
        '--readFilesCommand zcat '
        '--readFilesIn {input.r1} '
        # BAM file in transcript coords, in addition to genomic BAM file.
        '--quantMode GeneCounts '
        # Basic 2-pass mapping, with all 1st pass junctions inserted
        # into the genome indices on the fly.
        '--twopassMode Basic '
        # By default, this prefix is "./".
        '--outFileNamePrefix {STAR_OUT}trimmed/{wildcards.sample}. '
        # If exceeded, the read is considered unmapped.
        '--outFilterMultimapNmax 20 '
        # Minimum overhang for unannotated junctions.
        '--alignSJoverhangMin 8 '
        # Minimum overhang for annotated junctions.
        '--alignSJDBoverhangMin 1 '
        # Maximum number of mismatches per pair.
        '--outFilterMismatchNmax 999 '
        # Minimum intron length.
        '--alignIntronMin 1 '
        # Maximum intron length.
        '--alignIntronMax 1000000 '
        # Maximum genomic distance between mates.
        '--alignMatesGapMax 1000000 '
	    # output files
        #' 2> {output.progress}'
        '--outSAMtype BAM SortedByCoordinate '



#index bam

rule index_bam:
    input:
        bam=STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bai=STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'samtools index {input.bam} > {output.bai} '

#multiqc

rule multiqc_align:
    input:
        fastqc=FASTQC_TRIM_DIR,
        star_align=expand(STAR_OUT+"trimmed/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=config["samples"])
    output:
        report="multiqc_align/multiqc_report.html"
    params:
        input_dir=STAR_OUT
    conda:
        "envs/rna-seq_conda_env.yaml"
    shell:
        'multiqc -f {params.input_dir} {input.fastqc} -o multiqc_align '
