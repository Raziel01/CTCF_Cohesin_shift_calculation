
#!/bin/sh

read -p 'CTCF BED name: ' file1 ;
read -p 'Cohesin 1. BED name: ' file2 ;
read -p 'Cohesin 2. BED name: ' file3 ;
read -p 'CTCF BEDGRAPH.GZ name: ' file_graph1 ;
read -p 'Cohesin 1. BEDGRAPH.GZ name: ' file_graph2 ;
read -p 'Cohesin 2. BEDGRAPH.GZ name: ' file_graph3 ;
read -p 'Motif file name: ' motif_file ;

bed_base1=`basename $file1 .bed`;
bed_base2=`basename $file2 .bed`;
bed_base2=`basename $file3 .bed`;
bedgraph_base1=`basename $file_graph1 .bedgraph.gz`;
bedgraph_base2=`basename $file_graph2 .bedgraph.gz`;
bedgraph_base2=`basename $file_graph3 .bedgraph.gz`;

if [ ! -d tmp0120231 ] ; then
        /bin/mkdir -p tmp0120231
fi

for ((i=1;i<=3;i+=1)); do
        if [ ! -d bed_base${i} ] ; then
                 /bin/mkdir -p bed_base${i}
fi
        v=`echo 'file'$i`;
        peak=${!v};
        bname=`basename $peak .bed`
        v2=`echo 'file_graph'$i`;
        graph=${!v2};
        echo  $graph;
        cat ${peak} | sort -k1,1 -k2,2n  > tmp0120231/bed_base${i}_sorted.bed;
        PeakSplitter -p tmp0120231/bed_base${i}_sorted.bed  -w $graph -f -o bed_base${i};
        awk '{OFS="\t"; if (NR>1) print $1,$5,$5+1,$1":"$5"-"$5+1,$5}' bed_base${i}/*.subpeaks.bed > bed_base${i}_summits.bed;
        awk '{OFS="\t"; if (NR>1) print $1,$5,$5+1,$1":"$5"-"$5+1,$5}' bed_base${i}/*.subpeaks.bed > ${bname}_summits.bed
        awk '{OFS="\t"; print $1,$2-50,$3+50,$4,$5}' bed_base${i}_summits.bed > tmp0120231/bed_base${i}_summits_motker.tmp
        rm -rf bed_base${i};
done;


intersectBed -a $motif_file -b  tmp0120231/bed_base1_summits_motker.tmp -wa |  intersectBed -a stdin  -b  tmp0120231/bed_base2_summits_motker.tmp -wa |  intersectBed -a stdin  -b  tmp0120231/bed_base3_summits_motker.tmp -wa |   sort -k1,1 -k2,2n >   tmp0120231/CTCF_motif.bed;


for ((i=1;i<=3;i+=1)); do
        sort -k1,1 -k2,2n  bed_base${i}_summits.bed |  closestBed -a tmp0120231/CTCF_motif.bed  -b stdin  > tmp0120231/bed_base${i}_closest.bed;
        awk -F"\t" '{OFS="\t"; if ($6 == "+") print $4,$8-$2-8; else print $4,$2-$8+8}' tmp0120231/bed_base${i}_closest.bed  | sort -k1,1n  | awk '!n[$1]++' > tmp0120231/bed_base${i}_closest.tbl
done

if [  -d shift_table.tbl ] ; then
        /bin/cp -f shift_table.tbl shift_table_old.tbl
fi
echo -e "motifID\tCTCF_position\tCoh1_position\tCoh2_position" > shift_table.tbl;
paste tmp0120231/bed_base1_closest.tbl tmp0120231/bed_base2_closest.tbl tmp0120231/bed_base3_closest.tbl | awk -F"\t" '{OFS="\t"; if ($1 == $3 && $1 == $5  && $2 < 200 && $4 < 200 && $6 < 200 && $2 > -200 && $4 > -200 && $6 > -200) print $1,$2,$4,$6}' >> shift_table.tbl;
rm -rf tmp0120231 bed_base*_summits.bed
