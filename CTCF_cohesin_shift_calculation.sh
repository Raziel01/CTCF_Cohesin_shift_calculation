#!/bin/sh

read -p 'CTCF BED name: ' file1 ;
read -p 'Cohesin 1. BED name: ' file2 ;
read -p 'Cohesin 2. BED name: ' file3 ;
read -p 'CTCF BEDGRAPH.GZ name: ' file_graph1 ;
read -p 'Cohesin 1. BEDGRAPH.GZ name: ' file_graph2 ;
read -p 'Cohesin 2. BEDGRAPH.GZ name: ' file_graph3 ;
read -p 'Motif file name: ' motif_file ;

bed1_base=`basename $file1 .bed`;
bed2_base=`basename $file2 .bed`;
bed3_base=`basename $file3 .bed`;
bedgraph1_base=`basename $file_graph1 .bedgraph.gz`;
bedgraph2_base=`basename $file_graph2 .bedgraph.gz`;
bedgraph3_base=`basename $file_graph3 .bedgraph.gz`;

if [ ! -d tmp ] ; then
        /bin/mkdir -p tmp0120231
fi

for ((i=1;i<=3;i+=1)); do
        if [ ! -d bed${i}_base ] ; then
                 /bin/mkdir -p bed${i}_base
fi
        PeakSplitter -p file${i}  -w file_graph$i -f -o bed${i}_base;
        awk '{OFS="\t"; if (NR>1) print $1,$5,$5+1,$1":"$5"-"$5+1,$5}' bed${i}_base/*.subpeaks.bed > bed${i}_base_summits.bed;
        awk '{OFS="\t"; print $1,$2-50,$3+50,$4,$5}' bed${i}_base_summits.bed > tmp0120231/bed${i}_base_summits_motker.tmp
done;


intersectBed -a $motif_file -b  tmp/bed1_base_summits_motker.tmp -wa |  intersectBed -a stdin  -b  tmp/bed2_base_summits_motker.tmp -wa |  intersectBed -a stdin  -b  tmp/bed3_base_summits_motker.tmp -wa |   sort -k1,1 -k2,2n >   tmp0120231/CTCF_motif.bed;


for ((i=1;i<=3;i+=1)); do
        sort -k1,1 -k2,2n  bed${i}_base_summits.bed |  closestBed -a tmp/CTCF_motif.bed  -b stdin  > tmp0120231/bed${i}_base_closest.bed;
        awk -F"\t" '{OFS="\t"; if ($6 == "+") print $4,$8-$2-8; else print $4,$2-$8+8}' tmp0120231/bed${i}_base_closest.bed  | sort -k1,1n  | awk '!n[$1]++' > tmp0120231/bed${i}_base_closest.tbl
done

if [  -d shift_table.tbl ] ; then
        /bin/cp -f shift_table.tbl shift_table_old.tbl
fi
echo -e "motifID\tCTCF_position\tCoh1_position\tCoh2_position" > shift_results.tbl
paste tmp0120231/bed1_base_closest.tbl tmp0120231/bed2_base_closest.tbl tmp0120231/bed3_base_closest.tbl | awk -F"\t" '{OFS="\t"; if ($1 == $3 && $1 == $5 &&  && $2 < 200 && $4 < 200 && $6 < 200 && $2 > -200 && $4 > -200 && $6 > -200) print $1,$2,$4,$6}' >> shift_table.tbl;
rm -rf tmp0120231
