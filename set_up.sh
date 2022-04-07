#link the ocean dir
ln -s /ocean/projects/mcb180074p/sdai clip_proj

#get data eclip TARDBP K562
cd clip_proj
xargs -L 1 curl -O -J -L < files.txt

#rename the sequencing files
mv ENCFF734UEC.fastq.gz rep1.IP.umi.r1.fq.gz
mv ENCFF147JYD.fastq.gz rep1.IP.umi.r2.fq.gz
mv ENCFF661TYX.fastq.gz rep2.IP.umi.r1.fq.gz
mv ENCFF218BOC.fastq.gz rep2.IP.umi.r2.fq.gz