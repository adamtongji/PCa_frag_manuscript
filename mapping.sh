## env frag_3090

index=PATH_TO_BOWTIE2_GENOME
adapter=PATH_TO_ADAPTER

prefix=$1

outpath=$2

fq1=$3
fq2=$4
threads=10

java -jartrimmomatic-0.36.jar PE -phred33 $fq1 $fq2 ${outpath}/${prefix}_R1.fastq.gz ${outpath}/${prefix}_R1.unpair.fq.gz $outpath/${prefix}_R2.fastq.gz ${outpath}/${prefix}_R2.unpair.fq.gz  -threads $threads ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100


bowtie2 -q  -p ${threads} -x ${index} -1 ${outpath}/${prefix}_R1.fastq.gz -2 $outpath/${prefix}_R2.fastq.gz -S ${outpath}/${prefix}.sam
sambamba view -t $threads -S -f bam ${outpath}/${prefix}.sam > ${outpath}/${prefix}.bam 
sambamba index -t $threads ${outpath}/${prefix}.bam 
# bam_stat.py -i ${outpath}/${prefix}.bam > ${outpath}/${prefix}.stat 
sambamba sort -t $threads ${outpath}/${prefix}.bam -o ${outpath}/${prefix}.sort.bam

# samtools rmdup ${rmdup_mod} ${outpath}/${prefix}.sort.bam ${outpath}/${prefix}.sort.dedup.bam 2>  ${outpath}/${prefix}.flagstat.txt 
sambamba markdup -r  -t $threads ${outpath}/${prefix}.sort.bam ${outpath}/${prefix}.sort.dedup.bam 2>  ${outpath}/${prefix}.flagstat.txt

sambamba index -t $threads ${outpath}/${prefix}.sort.dedup.bam 
# bam_stat.py -i ${outpath}/${prefix}.sort.dedup.bam  1> ${outpath}/${prefix}.unique.count

rm  ${outpath}/${prefix}.sam ${outpath}/${prefix}.bam*  ${outpath}/${prefix}.sort.bam*

# sambamba view ${outpath}/${prefix}.sort.dedup.bam | awk '{if ($9>0 && $9<300 && $6=="150M"){printf "%s\t%d\t%d\t%d\n",$3,$4,$4+150,$9}}' | sort -k1,1 -k2,2g | gzip -c > ${outpath}/${prefix}.bed.gz

sambamba view -t 8 ${outpath}/${prefix}.sort.dedup.bam | awk '{if ($9>0 && $9<400){printf "%s\t%d\t%d\t%d\n",$3,$4,$4+150,$9}}' | sort -k1,1 -k2,2g | gzip -c > ${outpath}/${prefix}.bed.gz

# zcat ${outpath}/${prefix}.bed.gz | awk '{if ($4>260 || $4<150 && $4<400){print $0}}' | gzip -c > ${outpath}/${prefix}_fsc.bed.gz


window_file="/Data2/fusl/cnv_frag_analysis/fsc/ref/hg19_5m_window.sort.filter.bed"
intersectBed -a ${outpath}/${prefix}.bed.gz -b $window_file -wa -wb -sorted -f 0.51| gzip -c > ${outpath}/${prefix}_fsc_window.bed.gz

python /home/fusl/software/danbei_fragment/projects/early_ldt_v1/fs_table2.py ${prefix} $outpath

PATH_TO/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
${outpath}/${prefix}.sort.dedup.bam  > ${outpath}/${prefix}.wig

ichorCNA_path="PATH_TO/ichorCNA/"

Rscript ${ichorCNA_path}/scripts/runIchorCNA.R --id tumor_sample \
  --WIG ${outpath}/${prefix}.wig --ploidy "c(2,3)" --normal "c(0.95, 0.99, 0.995, 0.999)" --maxCN 5 \
  --gcWig ${ichorCNA_path}/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig ${ichorCNA_path}/inst/extdata/map_hg19_1000kb.wig \
  --centromere ${ichorCNA_path}/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --normalPanel ${ichorCNA_path}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $outpath

