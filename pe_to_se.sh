ls $1/*bam | sed 's/\.bam//g' > list

while read name
do


#Mate 1

samtools view -F4 -F1024 -F16 -f64 $name.bam | awk '{first=$1; second=0; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m1.f.std.sam
samtools view -F4 -f1024 -F16 -f64 $name.bam | awk '{first=$1; second=1024; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m1.f.dup.sam
samtools view -F4 -F1024 -f16 -f64 $name.bam | awk '{first=$1; second=16; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m1.r.std.sam
samtools view -F4 -f1024 -f16 -f64 $name.bam | awk '{first=$1; second=1040; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m1.r.dup.sam


#Mate 2: all reads need to be reversed here

samtools view -F4 -F1024 -F16 -f128 $name.bam | awk '{first=$1; second=16; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m2.f.std.sam
samtools view -F4 -f1024 -F16 -f128 $name.bam | awk '{first=$1; second=1040; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m2.f.dup.sam
samtools view -F4 -F1024 -f16 -f128 $name.bam | awk '{first=$1; second=0; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m2.r.std.sam
samtools view -F4 -f1024 -f16 -f128 $name.bam | awk '{first=$1; second=1024; $1=""; $2=""; print first, second, $0}' OFS="\t" | cut -f1,2,5- > $name.m2.r.dup.sam


cat <(samtools view -H $name.bam) $name.m1.f.std.sam $name.m1.f.dup.sam $name.m1.r.std.sam $name.m1.r.dup.sam $name.m2.f.std.sam $name.m2.f.dup.sam $name.m2.r.std.sam $name.m2.r.dup.sam | samtools sort | samtools view -bh > $name.se.bam

rm $name*sam

done <list

rm list

