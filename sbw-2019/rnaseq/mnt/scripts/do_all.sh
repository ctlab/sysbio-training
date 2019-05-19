# getting fastqs
mkdir fastqs
ln -sf /mnt/RNAseq/GSE120762/downsampled/*.fastq.gz fastqs/

TAGS=$(ls fastqs/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//')
echo $TAGS

# strting with only one sample
TAG=SRR7956038

# fastqc
for TAG in $TAGS; do
  OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR"
  fastqc -o "$OUTDIR" "fastqs/${TAG}_1.fastq.gz" |& tee "$OUTDIR/${TAG}_1.fastqc.log"
  
  OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR"
  fastqc -o "$OUTDIR" "fastqs/${TAG}_2.fastq.gz" |& tee "$OUTDIR/${TAG}_1.fastqc.log"
done
  
#======= Hisat2 pipeline ==========  
for TAG in $TAGS; do
  echo Processing $TAG
  
  HISAT_IDX=/mnt/RNAseq/reference/Gencode_mouse/release_M20/GRCm38.primary_assembly
  
  # aligning to the genome reference
  OUTDIR="hisat2/$TAG"; mkdir -p "$OUTDIR"
  date
  hisat2 -p 8 --new-summary -x  ${HISAT_IDX} \
    -1 "fastqs/${TAG}_1.fastq.gz" -2 "fastqs/${TAG}_2.fastq.gz" \
    2> "$OUTDIR/$TAG.log" \
    | samtools view -b - > "$OUTDIR/$TAG.raw.bam"
  date

  # post-processing the alignments
  
  date
  samtools sort -@ 8 -O bam "$OUTDIR/$TAG.raw.bam" > "$OUTDIR/$TAG.bam" && \
    samtools index "$OUTDIR/$TAG.bam" && \
    rm -v "$OUTDIR/$TAG.raw.bam"
  date
done  

for TAG in $TAGS; do
  echo Processing $TAG
  OUTDIR="hisat2/$TAG"
  
  # calculating coverage for vizualization
  bamCoverage -b "$OUTDIR/$TAG.bam" -o "$OUTDIR/$TAG.cov.bw" |& tee "$OUTDIR/$TAG.bamcov.log"
done

# QC  
for TAG in $TAGS; do
  echo Processing $TAG
  OUTDIR="hisat2/$TAG"

  REFGENE_MODEL=/mnt/RNAseq/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  infer_experiment.py -i "$OUTDIR/$TAG.bam" \
    -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.infer_experiment.txt"
  
  REFGENE_MODEL=/mnt/RNAseq/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  read_distribution.py -i "$OUTDIR/$TAG.bam" \
    -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.read_distribution.txt"
done

# Counting reads    
for TAG in $TAGS; do
  echo Processing $TAG
  
  GTF=/mnt/RNAseq/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf
  
  OUTDIR="featureCounts/$TAG"; mkdir -p "$OUTDIR"
  date
  featureCounts -a "$GTF" -s 0 -o "$OUTDIR/$TAG.fc.txt" \
    "hisat2/$TAG/$TAG.bam" |& tee "$OUTDIR/$TAG.log"
  date
done
  
  
#========== Kallisto ======================
for TAG in $TAGS; do
  echo Processing $TAG
  
  mkdir kallisto
  
  KALLISTO_IDX=/mnt/RNAseq/reference/Gencode_mouse/release_M20/gencode.vM20.transcripts.kalliso.idx
  
  OUTDIR="kallisto/$TAG"; mkdir -p "$OUTDIR"
  date
  # --single -l and -s option should be set for each dataset separately, 200+-50 is most common for single end
  kallisto quant -i $KALLISTO_IDX -t 8 \
    --rf-stranded \
    --plaintext \
    -o $OUTDIR \
    fastqs/${TAG}_1.fastq.gz fastqs/${TAG}_2.fastq.gz |& tee $OUTDIR/$TAG.log
  date

done

#========== multiqc for everything =============

multiqc -x .Rproj.user -f .

