#Extract DNA WGS from summary file from CGHUB for batch downloading. 

cd /Users/salpienowinski/Documents/prostate_cancer/TCGA/RNA_seq/gleason_6

mkdir ../../DNA_seq

mkdir ../../DNA_seq/gleason_6

scp /Users/salpienowinski/Documents/prostate_cancer/TCGA/RNA_seq/gleason_6/summary.tsv ../../DNA_seq/gleason_6

scp /Users/salpienowinski/Documents/prostate_cancer/TCGA/RNA_seq/gleason_6/manifest.xml ../../DNA_seq/gleason_6

cd /Users/salpienowinski/Documents/prostate_cancer/TCGA/DNA_seq/gleason_6/

cat summary.tsv | awk '$6 == "TP" && $8 == "DNA"  && $9 == "WGS" {print}' > DNA_gleason_6_unaligned_samples.txt

cat DNA_gleason_6_unaligned_samples.txt | cut -d- -f1-3 | uniq > nm.txt

cat nm.txt | awk '{print $2}' | while read line;  do grep $line < summary.tsv; done > all_out_match.txt

cat all_out_match.txt | awk '$8 == "DNA" && $9 == "WGS" {print}' | uniq > unique.txt

cat unique.txt | grep -v duplicat | grep Live > live_only.txt

cat live_only.txt | grep '01A' > tum.txt

cat live_only.txt | grep '10A\|11A' > norm.txt

cat tum.txt | awk '{print $2 "\t" $18}' > tum_sample_info_forCGHUB_download_dna.txt

cat norm.txt | awk '{print $2 "\t" $18}' > norm_sample_info_forCGHUB_download_dna.txt

paste tum_sample_info_forCGHUB_download_dna.txt norm_sample_info_forCGHUB_download_dna.txt > look_up_table.txt

cat tum_sample_info_forCGHUB_download_dna.txt | awk '{print "<analysis_id>"$2}' > tum_summary_ready_dna.txt

cat norm_sample_info_forCGHUB_download_dna.txt | awk '{print "<analysis_id>"$2}' > norm_summary_ready_dna.txt

awk '{print $1}' tum_summary_ready_dna.txt | while read line;  do grep -B 1 -A 24 $line < manifest.xml > manifest_tum_$line ; done

awk '{print $1}' norm_summary_ready_dna.txt | while read line;  do grep -B 1 -A 24 $line < manifest.xml > manifest_norm_$line ; done

for i in `ls manifest_tum_*`; do awk 'BEGIN{print "<ResultSet date=\"2015-15-06 08:33:03\">"}{print}END{print "</ResultSet>"}' $i > new_$i; done

for i in `ls manifest_norm_*`; do awk 'BEGIN{print "<ResultSet date=\"2015-15-06 08:33:03\">"}{print}END{print "</ResultSet>"}' $i > new_$i; done

# move files onto the cluster. the new_manifest* files to this path

scp -C new_manifest* nowinskicb@apollo.kcl.ac.uk:/home/nowinskicb/TCGA/DNA_seq/WGS/gleason_6/
/home/nowinskicb/TCGA/DNA_seq/WGS/gleason_6


### make qsub ready .sh files for downloading. 

cd $HOME/TCGA/DNA_seq/WGS/gleason_6/

for i in new_manifest*
do
    mv "$i" "`echo $i | sed 's/>//'`"
done

for i in new_manifest*
do
    mv "$i" "`echo $i | sed 's/<//'`"
done

awk '{print $2 "\t" $4}' look_up_table.txt |  while read -r a b; do printf '#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe threaded 8

module load cgpBattenberg/1.3.3

$HOME/gene_torrent/cghub/bin/gtdownload -k 454545 -vv -c $HOME/gene_torrent/cghub/cghub.key -d $HOME/TCGA/DNA_seq/WGS/gleason_6/new_manifest_tum_analysis_id'"$a"' --max-children 8\n
$HOME/gene_torrent/cghub/bin/gtdownload -k 454545 -vv -c $HOME/gene_torrent/cghub/cghub.key -d $HOME/TCGA/DNA_seq/WGS/gleason_6/new_manifest_norm_analysis_id'"$b"' --max-children 8\n

#is data paired end?
/apps/samtools/1.2/samtools flagstat '"$a"'/*.sorted.bam | grep 'duplicates' | awk '$1 + $3 > 0 {print}' > out_tum.txt
/apps/samtools/1.2/samtools flagstat '"$b"'/*.sorted.bam | grep 'duplicates' | awk '$1 + $3 > 0 {print}' > out_norm.txt


#remove duplicates first with samtools after checking if duplicates exist.

#if read -r && read -r

if wc -l out_tum.txt > 0
then
	/apps/samtools/1.2/samtools rmdup -S $a/*.sorted.bam $a/*.sorted.bam.dups_removed.bam
else
	echo " No Duplicates "
fi

if wc -l out_norm.txt > 0
then
	/apps/samtools/1.2/samtools rmdup -S $a/*.sorted.bam $a/*.sorted.bam.dups_removed.bam
else
	echo " No Duplicates "
fi < out.txt

/apps/cgpBattenberg/1.3.3/bin/battenberg.pl -r /home/ATRES/hg19/hg19.fa.fai -tb '"$a"'/*.sorted.bam -nb '"$a"'/*.sorted.bam -e /home/nowinskicb/ref/impute/impute_info.txt -o results_'"$a"'/ -u /home/nowinskicb/ref/1000genomesloci/


' > $HOME/gene_torrent/Scripts/WGS/gleason_6/$a.$b.sh; done



cd $HOME/gene_torrent/Scripts/WGS/gleason_6/
for i in `ls new*sh`; do echo qsub $i; done

