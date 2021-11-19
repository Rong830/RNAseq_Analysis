mkdir B200735
cp /localdisk/data/BPSM/AY21/fastq/* B200735
cp /localdisk/data/BPSM/AY21/Tcongo_genome/* B200735
cp /localdisk/data/BPSM/AY21/*.bed B200735
cp ask B200735
cd B200735
gunzip -r *.gz
cp /localdisk/data/BPSM/AY21/fastq/* .
#Unzipe all the file in directoty B200735.


cp /localdisk/data/BPSM/AY21/fastq/100k.fqfiles names.txt

mkdir quality_control
sed -i '1d' names.txt
#Make a txt file that contains all the samples' describtion, and remove the headline.
cut -f 6,7 names.txt  | while read file1 file2; do fastqc -t 12 -o quality_control --extract $file1 $file2 ; done
#Wrote a pipe line with a loop in it, using the names of samples to analyze every sample.
#Using fastqc to generate the quality repoerts of all samples, and output all of the result in a directory called quality_control.

bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta index
#Built a bowtie index with the gene given.
cut -f 1,6,7 names.txt  | while read name file1 file2;do bowtie2 -p 10 -x index -1 $file1 -2 $file2 | samtools view -u | samtools sort -@ 10 > $name.bam; done
#Using bowtie2 to algib sequences, and generate sam files.
#Using samtools to convert the formats form sam to bam, and then sort them.


cut -f 1 names.txt | while read name; do samtools index $name.bam; bedtools multicov -bams $name.bam -bed TriTrypDB-46_TcongolenseIL3000_2019.bed > $name.reads_count; echo $name been processing; done
#Using bedtools to counts the number of reads(compare in gene level)


touch Report-quality.txt
cd quality_control
mkdir result
mv *.html result
mv *.zip result
for name in `ls`; do
    cd $name
    for r in `cat summary.txt | cut -f 1`; do
        n=0;
        m=0;
        if [ $r = "FAIL" ]; then
            let n++
        fi
        if [ $r = "WARN" ]; then
            let m++
        fi
        if [ $n != 0 ] || [ $m != 0 ]; then
            echo "There is $n FAIL and $m WARN in $name! Please take this serious" >> ../Report-quality.txt
        fi
    done
    cd ..
done
cd ..
cp quality_control/Report-quality.txt .
#Generate a report for the samples' quality in summary.

cut -f2 names.txt | sort | uniq > type
cut -f4 names.txt | sort | uniq > time
cut -f5 names.txt | sort | uniq > treatment
cat type | while read type; do grep $type names.txt > $type.txt; done
cut -f 1,2,4,5 names.txt | while read name type time treatment; do cut -f6 $name.reads_count > $name-$type-$time-$treatment.txt; done
cat type | while read type; do paste *$type*.txt > $type.txt; done
cat type | while read type; do awk '{FS="\t"; if(FNR=1)columns=NF}{for(i=1;i<=columns;i++)sum+=$i;}{print sum/columns}{sum=0}' $type.txt > $type.totalmean.txt; done
for name in `cut -f1 names.txt`;do cut -f1,2,3,4,5 $name.reads_count > normal.txt;     break; done


clear
echo "***********************************************************************************************"
echo "**************************************Final report time!***************************************"
echo "***********************************************************************************************"

for type in `cat type`; do
    for treatment in `cat treatment`; do
        for group1 in `ls | find *$type*$treatment*`; do
            awk '{FS="\t"; if(FNR=1)columns=NF}{for(i=1;i<=columns;i++)sum+=$i;}{print sum/columns}{sum=0}' $group1 > $type-$treatment-AverageReads.txt
            echo $type-$treatment >> type-treatment-List.txt
        done
    done
done

cat type-treatment-List.txt | sort | uniq > type-treatment-List-sorted.txt

is_true=false
for type1 in `cat type-treatment-List-sorted.txt`; do
    for type2 in `cat type-treatment-List-sorted.txt`; do
        if [ $is_true = 'true' ]; then
            paste <(cut -f1 $type1-AverageReads.txt) <(cut -f1 $type2-AverageReads.txt) | awk '{FS="\t";{print FNR"\t"log($1+1)/log(2)-log($2+1)/log(2)}}' > $type1-$type2-FoldChangeMean.txt
            awk '{FS="\t";if($2 >= 1 || $2 <= -1){print $1,$2}}' $type1-$type2-FoldChangeMean.txt > ExtremeDiff-$type1-$type2-FoldChangeMean.txt
            echo -e "the number of extreme different expression level gene betweent $type1 and $type2 is: `cat ExtremeDiff-$type1-$type2-FoldChangeMean.txt | wc -l`" >> Report-count.txt
        fi
        if [ $type1 = $type2 ]; then
            is_true=true
        fi
    done
    is_true=false
done
#Caculate the mean expression level and the fold-change(plus 1 to each counts, and caculate the fc by log2) between groups.
#Count the number of extreme different expressed gene between different groups.

echo "***********************************************************************************************"


echo Compare them between different treatment!
for type in `cat type`; do
    for treatment in `cat treatment`; do
        for time in `cat time`; do
            for group1 in `ls | find *$type*$time*$treatment*`; do
                paste *$type*$time*$treatment* | awk '{FS="\t"; if(FNR=1)columns=NF}{for(i=1;i<=columns;i++)sum+=$i;}{print sum/columns}{sum=0}' > $type-$time-$treatment-AverageReads.txt
                echo $type-$time-$treatment >> type-time-treatment-List.txt
            done
        done
    done
done

cat type-time-treatment-List.txt | sort | uniq > type-time-treatment-List-sorted.txt

is_true=false
for time in `cat time`; do
    for type in `cat type`; do
        for treatment1 in `cat treatment`; do
            for treatment2 in `cat treatment`; do
                if [ $is_true = 'true' ]; then
                    if [ $time != 0 ]; then
                        paste <(cut -f1 $type-$time-$treatment1-AverageReads.txt) <(cut -f1 $type-$time-$treatment2-AverageReads.txt) | awk '{FS="\t";{print FNR"\t"log($1+1)/log(2)-log($2+1)/log(2)}}' > $type-$time-$treatment1-$treatment2-FoldChangeMean.txt
                        awk '{FS="\t";if($2 >= 1 || $2 <= -1){print $1,$2}}' $type-$time-$treatment1-$treatment2-FoldChangeMean.txt > ExtremeDiff-$type-$time-$treatment1-$treatment2-FoldChangeMean.txt
                        echo -e "the number of extreme different expression level gene betweent $type-$time-$treatment1 and $type-$time-$treatment2 is: `cat ExtremeDiff-$type-$time-$treatment1-$treatment2-FoldChangeMean.txt | wc -l `" >> Report-count.txt
                        paste normal.txt <(cut -f2 $type-$time-$treatment1-$treatment2-FoldChangeMean.txt) | awk '{FS="\t"; {print $6*$6,$0;}}' | sort -k1 -r -g | cut -d " " -f 2- > Report-$type-$time-$treatment1-$treatment2-FoldChangeMean.txt
                        sed -i '1i Chrom \t Start \t End \t Name \t Description \t log2(the mean of reads +1)-log2(the mean of reads +1)' Report-$type-$time-$treatment1-$treatment2-FoldChangeMean.txt
                    fi
                fi
                if [ $treatment1 = $treatment2 ]; then
                    is_true=true
                fi
            done
            is_true=false
        done
    done
done
echo "***********************************************************************************************"


echo Compare them between treating time!
is_true=false
for type in `cat type`; do
    for treatment in `cat treatment`; do
        for time1 in `cat time`; do
            for time2 in `cat time`;  do
                if [ $is_true = 'true' ]; then
                    if [ $time1 != 0 ] && [ $time2 != 0 ]; then
                        paste <(cut -f1 $type-$time1-$treatment-AverageReads.txt) <(cut -f1 $type-$time2-$treatment-AverageReads.txt) | awk '{FS="\t";{print FNR"\t"log($1+1)/log(2)-log($2+1)/log(2)}}' > $type-$time1-$time2-$treatment-FoldChangeMean.txt
                        awk '{FS="\t";if($2 >= 1 || $2 <= -1){print $1,$2}}' $type-$time1-$time2-$treatment-FoldChangeMean.txt > ExtremeDiff-$type-$time1-$time2-$treatment-FoldChangeMean.txt
                        echo -e "the number of extreme different expression level gene betweent $type-$time1-$treatment and $type-$time2-$treatment is: `cat ExtremeDiff-$type-$time1-$time2-$treatment-FoldChangeMean.txt | wc -l`" >> Report-count.txt
                        paste normal.txt <(cut -f2 $type-$time1-$time2-$treatment-FoldChangeMean.txt) | awk '{FS="\t"; {print $6*$6,$0;}}' | sort -k1 -r -g | cut -d " " -f 2-  > Report-$type-$time1-$time2-$treatment-FoldChangeMean.txt
                        sed -i '1i Chrom \t Start \t End \t Name \t Description \t log2(the mean of reads +1)-log2(the mean of reads +1)' Report-$type-$time1-$time2-$treatment-FoldChangeMean.txt
                    fi
                fi
                if [ $time1 = $time2 ]; then
                    is_true=true
                fi
            done
            is_true=false
        done
    done
done
echo "***********************************************************************************************"

echo Compare them between groups!
is_true=false
for time in `cat time`; do
    for treatment in `cat treatment`; do
        for type1 in `cat type`; do
            for type2 in `cat type`;  do
                if [ $is_true = 'true' ]; then
                    if [ $time != 0 ]; then
                        paste <(cut -f1 $type1-$time-$treatment-AverageReads.txt) <(cut -f1 $type2-$time-$treatment-AverageReads.txt) | awk '{FS="\t";{print FNR"\t"log($1+1)/log(2)-log($2+1)/log(2)}}' > $type1-$type2-$time-$treatment-FoldChangeMean.txt
                        awk '{FS="\t";if($2 >= 1 || $2 <= -1){print $1,$2}}' $type1-$type2-$time-$treatment-FoldChangeMean.txt > ExtremeDiff-$type1-$type2-$time-$treatment-FoldChangeMean.txt
                        echo -e "the number of extreme different expression level gene betweent $type1-$time-$treatment and $type2-$time-$treatment is: `cat ExtremeDiff-$type1-$type2-$time-$treatment-FoldChangeMean.txt | wc -l`" >> Report-count.txt
                        paste normal.txt <(cut -f2 $type1-$type2-$time-$treatment-FoldChangeMean.txt) | awk '{FS="\t"; {print $6*$6,$0;}}' | sort -k1 -r -g | cut -d " " -f 2-  > Report-$type1-$type2-$time-$treatment-FoldChangeMean.txt
                        sed -i '1i Chrom \t Start \t End \t Name \t Description \t log2(the mean of reads +1)-log2(the mean of reads +1)' Report-$type1-$type2-$time-$treatment-FoldChangeMean.txt
                    fi
                fi
                if [ $type1 = $type2 ]; then
                    is_true=true
                fi
            done
            is_true=false
        done
    done
done

mkdir data
mkdir reads_count
mv *.reads_count reads_count
cd ..
chmod 700 ask
cp ask B200735
./ask

