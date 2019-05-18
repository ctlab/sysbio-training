# getting fastqs
mkdir fastqs
ln -s /mnt/rnaseq/fastqs/*.fastq.gz fastqs/

# strting with only one sample
TAG=MPH_untr_rep1

# fastqc

OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR"
fastqc -o "$OUTDIR" "fastqs/$TAG.fastq.gz" |& tee "$OUTDIR/$TAG.fastqc.log"

#======= Hisat2 pipeline ==========

HISAT_IDX=/mnt/reference/Gencode_mouse/release_M20/GRCm38.primary_assembly

# aligning to the genome reference

OUTDIR="hisat2/$TAG"; mkdir -p "$OUTDIR"
date
hisat2 -p 4 --new-summary -x  ${HISAT_IDX} \
  -U "fastqs/$TAG.fastq.gz" \
  2> "$OUTDIR/$TAG.hisat2.log" \
  | samtools view -b - > "$OUTDIR/$TAG.raw.bam"
date

ls $OUTDIR

# post-processing the alignments

date
samtools sort -@ 4 -O bam "$OUTDIR/$TAG.raw.bam" > "$OUTDIR/$TAG.bam" && \
  samtools index "$OUTDIR/$TAG.bam" && \
  rm -v "$OUTDIR/$TAG.raw.bam"
date

# calculating coverage for vizualization
bamCoverage -b "$OUTDIR/$TAG.bam" -o "$OUTDIR/$TAG.cov.bw" |& tee "$OUTDIR/$TAG.bamcov.log"

# QC

REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
infer_experiment.py -i "$OUTDIR/$TAG.bam" \
  -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.infer_experiment.txt"

REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
read_distribution.py -i "$OUTDIR/$TAG.bam" \
  -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.read_distribution.txt"

# REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10.HouseKeepingGenes.bed
# geneBody_coverage.py \
#   -i $OUTDIR/$TAG.bam \
#   -o $OUTDIR/$TAG \
#   -r $REFGENE_MODEL 
  

# Counting reads

GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf

OUTDIR="featureCounts/$TAG"; mkdir -p "$OUTDIR"
date
featureCounts -a "$GTF" -s 0 -o "$OUTDIR/$TAG.fc.txt" \
  "hisat2/$TAG/$TAG.bam" |& tee "$OUTDIR/$TAG.fc.log"
date

head "$OUTDIR/$TAG.fc.txt"
wc -l "$OUTDIR/$TAG.fc.txt"


#========== Kallisto ======================

mkdir kallisto

KALLISTO_IDX=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.transcripts.kalliso.idx

OUTDIR="kallisto/$TAG"; mkdir -p "$OUTDIR"
date
# --single -l and -s option should be set for each dataset separately, 200+-50 is most common for single end
kallisto quant -i $KALLISTO_IDX -t 4 \
  --single -l 200 -s 50 \
  --plaintext \
  -o $OUTDIR \
  fastqs/$TAG.fastq.gz |& tee $OUTDIR/$TAG.kallisto.log
date

#========== multiqc for everything =============

multiqc -x .Rproj.user -f .


#========== mmquant ============
OUTDIR="mmquant/$TAG"; mkdir -p "$OUTDIR"
GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf
date
mmquant -a "$GTF" -s U -o "$OUTDIR/$TAG.mmq.txt" \
  -r "hisat2/$TAG/$TAG.bam" |& tee "$OUTDIR/$TAG.mmq.log"
date