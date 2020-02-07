# Minute pipeline test data

Created from these files:

P1_H3K4m3_R1.fastq.gz
P1_H3K4m3_R2.fastq.gz
P1_Inp-PE_R1.fastq.gz
P1_Inp-PE_R2.fastq.gz


## Reference

See the `Snakefile`


# Test data set

Map all reads to mm9 with BWA-MEM. It will soft-clip the UMI and barcode, so they can be left in:

    bwa mem -t 16 /sw/data/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/BWAIndex/genome.fa fastq/P1_H3K4m3_R* | samtools view -b - > h3k4m3.bam
    bwa mem -t 16 /sw/data/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/BWAIndex/genome.fa fastq/P1_Inp-PE_R* | samtools view -b - > inp.bam

Pick names of all reads that map to chr19, position <= 3.7 Mbp:

    samtools view h3k4m3.bam | awk '$3==19 && $4 <= 3700000 {print $1}' | sort -u > h3k4m3.target.txt
    samtools view inp.bam | awk '$3==19 && $4 <= 3700000 {print $1}' | sort -u > inp.target.txt

A couple of unmapped reads:

    for ds in h3k4m3 inp; do samtools view -f 4 $ds.bam | cut -f1|uniq|head -n 100000|sort -u|tail -n 10000 > $ds.unmapped.txt; done

    cat h3k4m3.{target,unmapped}.txt | sort -u > h3k4m3.readnames.txt
    cat inp.{target,unmapped}.txt | sort -u > inp.readnames.txt

Pick reads from original FASTQ:

    for r in 1 2; do zgrep -A 3 --no-group-separator -F -f h3k4m3.readnames.txt fastq/P1_H3K4m3_R$r.fastq.gz | pigz > h3k4m3_R$r.fastq.gz; done
    for r in 1 2; do zgrep -A 3 --no-group-separator -F -f inp.readnames.txt fastq/P1_Inp-PE_R$r.fastq.gz | pigz > inp_R$r.fastq.gz; done

