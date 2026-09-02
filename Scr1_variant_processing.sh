#!/bin/bash
# ==============================================================================
# Pipeline: CHIP Variant Processing, Filtering, ANNOVAR Annotation, and VCF Generation
# Description: Processes MosaicHunter (SNV) and Pisces (Indel) outputs, applies
#              depth filtering, Fisher's exact test filtering, ANNOVAR functional
#              annotation, and generates final VCF/summary files without sbatch.
#              "AD_panel.sample.clinical.tsv" has age information before updates.
#              Use age information in "SupplTable2_Somatic_Variant_Calls.xlsx".
# ==============================================================================

set -euo pipefail

# Configuration & Paths
REFERENCE_DIR="${REFERENCE_DIR}"
ANNOVAR_DIR="${ANNOVAR_DIR}"
CUTOFF_DEPTH=100

echo "=== Step 1: Generating Sample Lists and Coverage Filtering ==="
mkdir -p logs ANNOVAR MosaicHunter Pisces

# Calculate coverage summaries from BED coverage files
grep "done" AD_panel.sample.list | cut -f1 | while read -r sample; do
    echo "$sample"
    awk '$1=="all"{sum+=$2*$5} END{print sum}' finalBam/$sample.final.coverageBed
    awk '$1=="all"{if($2>200){sum+=$5}} END{print sum}' finalBam/$sample.final.coverageBed
done | paste - - - > AD_panel.coverage.summary

# Filter samples by target depth cutoff and merge clinical data
awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2>cut{print $1}' summary_on_target_depth_Batch1+2.tsv | awk -F'[/.]' '{print $2}' > AD_panel.sample.included.list.txt
awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2>cut{print $1}' summary_on_target_depth_Batch3.tsv | awk -F'[/.]' '{print $2}' >> AD_panel.sample.included.list.txt
awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2>cut{print $1}' summary_on_target_depth_Batch4.tsv | awk -F'[/.]' '{print $2}' >> AD_panel.sample.included.list.txt

awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2<cut{print $1}' summary_on_target_depth_Batch1+2.tsv | awk -F'[/.]' '{print $2}' > AD_panel.sample.excluded.list.txt
awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2<cut{print $1}' summary_on_target_depth_Batch3.tsv | awk -F'[/.]' '{print $2}' >> AD_panel.sample.excluded.list.txt
awk -v cut="$CUTOFF_DEPTH" -F'\t' '$2<cut{print $1}' summary_on_target_depth_Batch4.tsv | awk -F'[/.]' '{print $2}' >> AD_panel.sample.excluded.list.txt

cat <(awk -F'\t' '{print $0"\tBatch1+2"}' ../Batch1+2_SMC_panel/AD_panel.sample.clinical.tsv) \
    <(awk -F'\t' '{print $0"\tBatch3"}' ../Batch3_SMC_panel/AD_panel.sample.clinical.tsv) \
    <(awk -F'\t' '{print $0"\tBatch4"}' ../Batch4_SMC_panel/AD_panel.sample.clinical.tsv) | \
    awk -F'\t' 'NR==FNR {valid[$1]=1; next} valid[$1]' AD_panel.sample.included.list.txt - > AD_panel.sample.clinical.tsv

cat AD_panel.sample.list | awk -F'\t' 'NR==FNR {valid[$1]=1; next} valid[$1]' AD_panel.sample.included.list.txt - > AD_panel.sample.filtered.list

echo "=== Step 2: Processing SNVs (MosaicHunter Raw & Filtered Checks) ==="
grep "done" AD_panel.sample.filtered.list | cut -f1 | while read -r sample; do
    # Raw check script logic (scr_summarize_div_p1.sh)
    awk '$7>3 && $7/($6+$7+$8)<0.3 && 10^$12>0.5' MosaicHunter/$sample.MH.tsv | cut -f1,2 | while read -r chr pos; do
        echo "$chr"; echo "$pos"
        for group in "CTRL" "AD" "C" ""; do
            cat AD_panel.sample.clinical.tsv | awk -F "\t" -v g="$group" '$6==g' | cut -f1 | while read -r s3; do
                awk '$7>0' MosaicHunter/$s3.MH.tsv
            done | awk -v c="$chr" -v p="$pos" '$1==c && $2==p' | wc -l
        done
    done | paste - - - - - - > MosaicHunter/$sample.raw.check

    # Filtered check script logic (scr_summarize_div_p2.sh)
    awk '$7>3 && $7/($6+$7+$8)<0.3 && 10^$12>0.5' MosaicHunter/$sample.MH.tsv | cut -f1,2 | while read -r chr pos; do
        echo "$chr"; echo "$pos"
        for group in "CTRL" "AD" "C" ""; do
            cat AD_panel.sample.clinical.tsv | awk -F "\t" -v g="$group" '$6==g' | cut -f1 | while read -r s3; do
                awk '$7/($6+$7+$8)<0.3 && 10^$12>0.5' MosaicHunter/$s3.MH.tsv
            done | awk -v c="$chr" -v p="$pos" '$1==c && $2==p' | wc -l
        done
    done | paste - - - - - - > MosaicHunter/$sample.filtered.check
done

echo "=== Step 3: Statistical Filtering & Annotation for SNVs ==="
N_CTRL=$(awk -F "\t" '$6=="CTRL"||$6=="C"' AD_panel.sample.clinical.tsv | wc -l)
N_AD=$(awk -F "\t" '$6=="AD"' AD_panel.sample.clinical.tsv | wc -l)
echo "Processing SNVs with N_CTRL=$N_CTRL and N_AD=$N_AD"

cat AD_panel.sample.clinical.tsv | cut -f1 | while read -r sample; do
    python3 fisher_filter_v2_2sided_0.1proportion.py MosaicHunter/$sample.filtered.check "$N_CTRL" "$N_AD" 0.05 0.1 > fisher_snv_tmp.out
    my.grep -k 1,2 -c 1,2 -q fisher_snv_tmp.out -f MosaicHunter/$sample.MH.tsv > MosaicHunter/$sample.MH.loose.tsv
done
rm -f fisher_snv_tmp.out
echo "SNV Loose counting completed."

# Generate SNV summary files
myjoin -m -F1 <(cat AD_panel.sample.clinical.tsv | awk -F "\t" '$6!=""' | cut -f1 | while read -r f; do echo "$f"; wc -l MosaicHunter/$f.MH.loose.tsv | cut -f1 -d " "; done | paste - -) -f1 AD_panel.sample.clinical.tsv | cut -f1,2,9 > AD_panel.MH.loose.summary

# SNV ANNOVAR Annotations
mkdir -p ANNOVAR
cat AD_panel.sample.clinical.tsv | cut -f1 | while read -r f; do awk -v ID=$f '{OFS="\t"; print $0,ID}' MosaicHunter/$f.MH.loose.tsv; done > ANNOVAR/AD_panel.MH.loose.input
annotate_variation.pl --geneanno --dbtype refgene --buildver hg19 --outfile ANNOVAR/AD_panel.MH.loose <(awk '{OFS="\t"; print $1,$2,$2,$3,$4}' ANNOVAR/AD_panel.MH.loose.input) "${ANNOVAR_DIR}/"
annotate_variation.pl --filter --dbtype gnomad_genome --buildver hg19 --outfile ANNOVAR/AD_panel.MH.loose <(awk '{OFS="\t"; print $1,$2,$2,$3,$4}' ANNOVAR/AD_panel.MH.loose.input) "${ANNOVAR_DIR}/"
annotate_variation.pl --filter --dbtype ljb26_all -otherinfo --buildver hg19 --outfile ANNOVAR/AD_panel.MH.loose <(awk '{OFS="\t"; print $1,$2,$2,$3,$4}' ANNOVAR/AD_panel.MH.loose.input) "${ANNOVAR_DIR}/"
annotate_variation.pl --filter --dbtype cosmic70 --buildver hg19 --outfile ANNOVAR/AD_panel.MH.loose <(awk '{OFS="\t"; print $1,$2,$2,$3,$4}' ANNOVAR/AD_panel.MH.loose.input) "${ANNOVAR_DIR}/"

myjoin -F1,2,3,4 ANNOVAR/AD_panel.MH.loose.input -f3,5,6,7 ANNOVAR/AD_panel.MH.loose.hg19_gnomad_genome_dropped | awk '{OFS="\t"; if($1=="="){print $2,$3,$4,$5,$6,$7,$8,$9,$17,$19} if($1=="+"){print $2,$3,$4,$5,$6,$7,$8,$9,$17,0}}' | uniq > ANNOVAR/AD_panel.MH.loose.tmp1
myjoin -F1,2,3,4 ANNOVAR/AD_panel.MH.loose.tmp1 -f3,5,6,7 ANNOVAR/AD_panel.MH.loose.variant_function | awk '$1=="="||$1=="+"' | cut -f2-13 | uniq > ANNOVAR/AD_panel.MH.loose.tmp2
myjoin -F1,2,3,4 ANNOVAR/AD_panel.MH.loose.tmp2 -f4,6,7,8 ANNOVAR/AD_panel.MH.loose.exonic_variant_function | awk '$1=="="||$1=="+"' | cut -f2-13,15,16 | sed -e 's/ SNV//g' | uniq > AD_panel.MH.loose.tsv

cat <(cat ANNOVAR/AD_panel.MH.loose.hg19_ljb26_all_dropped | awk '{split($2,array,","); if(array[2]=="D"||array[4]=="D"){print $3"\t"$4}}') \
    <(awk '$1=="splicing"' ANNOVAR/AD_panel.MH.loose.variant_function | cut -f3,4) \
    <(awk '$2~"stop"' ANNOVAR/AD_panel.MH.loose.exonic_variant_function | cut -f4,5) | sort | uniq > ANNOVAR/AD_panel.MH.loose.deleterious

myjoin -F1,2 <(myjoin -F3,5,7 ANNOVAR/AD_panel.MH.loose.variant_function -f4,6,8 ANNOVAR/AD_panel.MH.loose.exonic_variant_function | awk -F "\t" '{OFS="\t";split($3,a,"(");symbol=a[1];if($1=="="){print $4,$6,$7,$8,symbol,$2,$10}if($1=="+"){print $4,$6,$7,$8,symbol,$2,"NA"}}') -f1,2 ANNOVAR/AD_panel.MH.loose.deleterious | awk -F "\t" '{OFS="\t";if($1=="="){print $2,$3,$4,$5,$6,$7,$8,"deleterious"}if($1=="+"){print $2,$3,$4,$5,$6,$7,$8,"neutral"}}' | sed -e 's/ SNV//g' | sort | uniq > AD_panel.MH.loose.gene_anno.tsv
rm -f ANNOVAR/AD_panel.MH.loose.tmp* ANNOVAR/AD_panel.MH.loose.log

cat <(echo "#CHROM_POS_ID_REF_ALT_QUAL_FILTER" | sed -e 's/_/\t/g') \
    <(cut -f1-4 AD_panel.MH.loose.tsv | sort -k1,1n -k2,2n | uniq | awk '{OFS="\t";print "chr"$1,$2,$1":"$2"_"$3"/"$4,$3,$4,".","PASS"}') > AD_panel.MH.loose.vcf

echo "=== Step 4: Processing Indels ==="
grep "done" AD_panel.sample.filtered.list | cut -f1 | while read -r sample; do
    grep PASS Pisces/$sample.final.vcf | awk '$10~"0/1"' | awk 'length($4)>1||length($5)>1' | awk '{split($10,a,":");split(a[3],b,",");if(b[2]>3&&a[5]<0.3){print $0}}' | cut -f1,2 | while read -r f1 f2; do
        echo "$f1"; echo "$f2"
        for group in "CTRL" "C" "AD" ""; do
            cat AD_panel.sample.clinical.tsv | awk -F "\t" -v g="$group" '$6==g' | cut -f1 | while read -r f3; do
                awk 'length($4)>1||length($5)>1' Pisces/$f3.final.vcf
            done | awk -v chr=$f1 -v pos=$f2 '$1==chr&&$2==pos' | wc -l
        done
    done | paste - - - - - - > Pisces/$sample.raw.check
done

cat AD_panel.sample.clinical.tsv | cut -f1 | while read -r f; do
    python3 fisher_filter_indel_v2_2sided_0.1proportion.py Pisces/$f.raw.check "$N_CTRL" "$N_AD" 0.05 0.1 > fisher_indel_tmp.out
    my.grep -k 1,2 -c 1,2 -q fisher_indel_tmp.out -f Pisces/$f.final.vcf > Pisces/$f.loose.vcf
done
rm -f fisher_indel_tmp.out
echo "Indel loose counting completed."

# Indel ANNOVAR Annotations & Summaries
myjoin -m -F1 <(cat AD_panel.sample.clinical.tsv | awk -F "\t" '$6!=""' | cut -f1 | while read -r f; do echo "$f"; wc -l Pisces/$f.loose.vcf | cut -f1 -d " "; done | paste - -) -f1 AD_panel.sample.clinical.tsv | cut -f1,2,8 > AD_panel.indel.loose.summary

cat AD_panel.sample.clinical.tsv | awk -F "\t" '$6!=""' | cut -f1 | while read -r f; do
    convert2annovar.pl -format vcf4 -includeinfo Pisces/$f.loose.vcf | awk -v ID=$f '{OFS="\t"; print $1,$2,$3,$4,$5,$15,ID}'
done | awk '$4=="-"||$5=="-"' | awk '{split($6,a,":");split(a[3],b,",");if(b[2]>3&&a[5]<0.3){print $0}}' > ANNOVAR/AD_panel.indel.loose.input

annotate_variation.pl --geneanno --dbtype refgene --buildver hg19 --outfile ANNOVAR/AD_panel.indel.loose <(cut -f1-5 ANNOVAR/AD_panel.indel.loose.input) "${ANNOVAR_DIR}/"
annotate_variation.pl --filter --dbtype gnomad_genome --buildver hg19 --outfile ANNOVAR/AD_panel.indel.loose <(cut -f1-5 ANNOVAR/AD_panel.indel.loose.input) "${ANNOVAR_DIR}/"

myjoin -F1,2,3,4,5 ANNOVAR/AD_panel.indel.loose.input -f3,4,5,6,7 ANNOVAR/AD_panel.indel.loose.hg19_gnomad_genome_dropped | awk '{OFS="\t"; if($1=="="){print $2,$3,$4,$5,$6,$7,$8,$10} if($1=="+"){print $2,$3,$4,$5,$6,$7,$8,0}}' | uniq > ANNOVAR/AD_panel.indel.loose.tmp1
myjoin -F1,2,3,4,5 ANNOVAR/AD_panel.indel.loose.tmp1 -f3,4,5,6,7 ANNOVAR/AD_panel.indel.loose.variant_function | awk '$1=="="||$1=="+"' | cut -f2-11 | uniq > ANNOVAR/AD_panel.indel.loose.tmp2
myjoin -F1,2,3,4,5 ANNOVAR/AD_panel.indel.loose.tmp2 -f4,5,6,7,8 ANNOVAR/AD_panel.indel.loose.exonic_variant_function | awk '$1=="="||$1=="+"' | cut -f2-11,13,14 | uniq > AD_panel.indel.loose.tsv

cat <(awk '$1=="splicing"' ANNOVAR/AD_panel.indel.loose.variant_function | cut -f3,4,5) \
    <(awk '$2~"^frameshift"||$2~"stop"' ANNOVAR/AD_panel.indel.loose.exonic_variant_function | cut -f4,5,6) | sort | uniq > ANNOVAR/AD_panel.indel.loose.deleterious

myjoin -F1,2,3 <(myjoin -F3,4,5,6,7 ANNOVAR/AD_panel.indel.loose.variant_function -f4,5,6,7,8 ANNOVAR/AD_panel.indel.loose.exonic_variant_function | awk -F "\t" '{OFS="\t";split($3,a,"(");symbol=a[1];if($1=="="){print $4,$5,$6,$7,$8,symbol,$2,$10}if($1=="+"){print $4,$5,$6,$7,$8,symbol,$2,"NA"}}') -f1,2,3 ANNOVAR/AD_panel.indel.loose.deleterious | awk -F "\t" '{OFS="\t";if($1=="="){print $2,$3,$4,$5,$6,$7,$8,$9,"deleterious"}if($1=="+"){print $2,$3,$4,$5,$6,$7,$8,$9,"neutral"}}' | sort | uniq > AD_panel.indel.loose.gene_anno.tsv
rm -f ANNOVAR/AD_panel.indel.loose.tmp* ANNOVAR/AD_panel.indel.loose.log

echo "Complete variant processing pipeline finished successfully."
