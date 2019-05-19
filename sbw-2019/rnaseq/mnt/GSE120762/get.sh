mkdir raw
prefetch SRR79560{38,39,40,41,42,43
cd raw
fastq-dump -v --gzip --split-files SRR79560{38,39,40,41,42,43}
cd -
mkdir downsampled

for f in raw/*.gz; do
    seqtk sample  $f 0.1 | gzip > downsampled/$(basename $f)
done
