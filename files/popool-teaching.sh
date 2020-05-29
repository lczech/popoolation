# trim reads
perl $1/basic-pipeline/trim-fastq.pl --input1 read_1.fastq --input2 read_2.fastq --output trim --quality-threshold 20 --min-length 50

# prepare the reference
mkdir wg
mv dmel-2R-chromosome-r5.22.fasta wg
awk '{print $1}' wg/dmel-2R-chromosome-r5.22.fasta > wg/dmel-2R-short.fa
bwa index wg/dmel-2R-short.fa

# map the reads
bwa aln wg/dmel-2R-short.fa trim_1 > trim_1.sai
# (bwa aln -l 100 -o 2 -d 12 -e 12 -n 0.01 -t 3 wg/dmel-2R-short.fa trim_1 > trim_1.sai)
bwa aln wg/dmel-2R-short.fa trim_2 > trim_2.sai
bwa sampe wg/dmel-2R-short.fa trim_1.sai trim_2.sai trim_1 trim_2 > maped.sam

# extract reads with a mapping quality of >= 20 (unambiguously mapped reads) and create a sorted bam file
samtools view -q 20 -bS maped.sam| samtools sort - maped.sort

# create a pileup file
samtools pileup maped.sort.bam > cyp6g1.pileup

# calculate Tajima's pi
perl $1/Variance-sliding.pl --measure pi --input cyp6g1.pileup --min-count 2 --min-qual 20 --min-coverage 4 --max-coverage 70 --pool-size 500 --window-size 1000 --step-size 1000 --output cyp6g1.varslid.pi --region 2R:7800000-8300000

# a short overview
perl $1/Visualise-output.pl --input cyp6g1.varslid.pi --output cyp6g1.pdf --chromosomes "2R" --ylab pi

# display results in IGV
perl $1/VarSliding2Wiggle.pl --input cyp6g1.varslid.pi --output cyp6g1.pi.wig --trackname "nat-pop-pi"
samtools index maped.sort.bam
# open IGV and load maped.sort.bam, cyp6g1.gtf and cyp6gi.pi.wig

# calculate genewise pi
perl $1/Variance-at-position.pl --measure pi --pileup cyp6g1.pileup --gtf cyp6g1.gtf --output genewise.pi --pool-size 500 --min-count 2 --min-coverage 4 --max-coverage 70 --min-qual 20
# sort (cyp6g1 is CG8453-RA)
sort -k 4,4n genewise.pi> sorted-genewise.pi

# Display the results
# less sorted-genewise.pi


