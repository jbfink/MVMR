#!/bin/bash

# /genetics/PAREG/perrotn/scripts/IVs_MVMR_PURE.sh



# DO NOT USE, USE THE ONE THAT INCLUDE BM AS ADDITIONNAL EXPOSURE






























########################################
########## Multivariate MR IV ##########
########################################

# MVMR IV FOR PURE BM ONLY. DOES NOT RETRIEVE IV FOR THE OTHER EXPOSURES
# FROM SIGNIFICANT pQTLs EXTRACT SAME SNPs IN OTHER EXPOSURES
# EXTRACT SAME SNPs IN OUTCOME
# CLUMPING
# ...

echo "Started IVs MVMR at: "$(date)
SECONDS=0

SUFFIX=$1
PANEL=$2
ASSAY=$3
GENE=$4
LIST_EXPOSURES=$5
PVALUE_THRESHOLD=$6
CIS_WINDOW=$7
LD_THRESHOLD_PRUNING=$8
OUTCOME=$9
ROOT_DIR=${10}
GRS_FORMATING_FILE=${11}
MODEL=${12}
POPULATION=${13}

#### TEST ####
# CIS_WINDOW=200000
# SUFFIX="NO_MHC_MISSENSE_SPLICING"
# PVALUE_THRESHOLD=0.01
# LD_THRESHOLD_PRUNING=0.1
# PANEL="CMET"
# ASSAY="F11"
# GENE="F11"
# POPULATION="EUROPEAN"
# ROOT_DIR="/storage/genetics_work2/perrotn/"
# GRS_FORMATING_FILE="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/GRS_FORMATTING_FILE.K2Uwt9"
# LIST_EXPOSURES="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/FXI_CES_VS_ISCHEMIC_STROKE.txt"
# OUTCOME="ischemic_stroke_european_GIGASTROKE_2022"
# MODEL="FXI_CES_VS_ISCHEMIC_STROKE"
##############


echo "SUFFIX:" $SUFFIX
echo "PANEL:" $PANEL
echo "ASSAY:" $ASSAY
echo "GENE:" $GENE
echo "LIST_EXPOSURES:" $LIST_EXPOSURES
echo "PVALUE_THRESHOLD:" $PVALUE_THRESHOLD
echo "CIS_WINDOW:" $CIS_WINDOW
echo "LD_THRESHOLD_PRUNING:" $LD_THRESHOLD_PRUNING
echo "OUTCOME:" $OUTCOME
echo "ROOT_DIR:" $ROOT_DIR
echo "GRS_FORMATING_FILE:" $GRS_FORMATING_FILE
echo "MODEL:" $MODEL
echo "POPULATION:" $POPULATION


export TMPDIR=${ROOT_DIR}tmp

ROOT_OUTPUT_DIR=${ROOT_DIR}PURE/MVMR_consortia/${SUFFIX}/
IVs_PURE=/storage/genetics_work2/perrotn/PURE/MR_INSTRUMENTS_V5/PURE_MR_IVs_PANEL_${PANEL}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${SUFFIX}
GRS_DIR="/genetics/PAREG/common/Public_Data_Organized/GWAS_FORMATTING_2022_02_02/"
plink2="/genetics/PAREG/common/PURE_DATA/MR_INSTRUMENTS/2020_05_06/plink"

mkdir -p $ROOT_OUTPUT_DIR

echo "##### OUTCOME:" $OUTCOME "#####"

OUTPUT_DIR=${ROOT_OUTPUT_DIR}${OUTCOME}/
BIOMARKER_OUTCOME_DIR=${OUTPUT_DIR}OUTCOME_FORMATTING_FOR_MRBASE
BIOMARKER_EXPOSURES_DIR=${OUTPUT_DIR}EXPOSURES_FORMATTING_FOR_MRBASE

mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR} || exit 1

mkdir -p ${BIOMARKER_OUTCOME_DIR}
mkdir -p ${BIOMARKER_EXPOSURES_DIR}

if [ ! -f $IVs_PURE ]; then
  echo "NO PURE IV FILE"
  exit 1
fi

#### Format Exposure & Outcome Summary Statistics ######
unset COL_NUM_PURE
declare -A COL_NUM_PURE
count=0
for i in $(head -n 1 $IVs_PURE); do
count=$((count+1))
COL_NUM_PURE[$i]=$count
done

awk 'BEGIN{FS=OFS="\t"}NR>1 &&  $('${COL_NUM_PURE["Phenotype"]}')=="'$POPULATION'""_""'$PANEL'""_""'$ASSAY'""_""'$GENE'" && (length($('${COL_NUM_PURE["other_allele"]}'))+length($('${COL_NUM_PURE["effect_allele"]}'))==2) {print $('${COL_NUM_PURE["chr"]}'), $('${COL_NUM_PURE["pos"]}'), $('${COL_NUM_PURE["other_allele"]}'), $('${COL_NUM_PURE["effect_allele"]}')}' $IVs_PURE | sort -u > list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}


echo $(echo $(wc -l list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} | awk '{print $1}' - ) "SIGNIFICANT SNPs" $POPULATION $PANEL $ASSAY $GENE)


echo "" | awk 'BEGIN{FS=OFS="\t";print "Phenotype","SNP","chr","pos","other_allele","effect_allele", "beta","se","eaf","pval","units","gene","samplesize","ncase","ncontrol"}' - > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}

unset COL_NUM
count=0
declare -A COL_NUM
for name in $(head -n 1 ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}) ; do
count=$((count+1))
COL_NUM[$name]=$count
done

unset COL_NUM_FOR
declare -A COL_NUM_FOR
count=0
for i in $(head -n 1 $GRS_FORMATING_FILE); do
count=$((count+1))
COL_NUM_FOR[$i]=$count
done

EXPOSURE_MODEL="ESUS_CCSc_EUR_SiGN_2019"
for EXPOSURE_MODEL in $(awk 'BEGIN{FS=OFS="\t"}NR>1{print}' $LIST_EXPOSURES)
do
  echo $EXPOSURE_MODEL

  if [ -f ${GRS_DIR}${EXPOSURE_MODEL}/4_final.txt.gz ]; then
      SUMSTAT_FILE=${GRS_DIR}${EXPOSURE_MODEL}/4_final.txt.gz
  else
    echo $EXPOSURE_MODEL | cat - >> ${ROOT_OUTPUT_DIR}MISSING_EXPOSURE_MODEL
    cd $ROOT_OUTPUT_DIR || exit 1
    rm -r $OUTPUT_DIR
    exit 1
  fi

  HEAD_SS=$(gzip -cd $SUMSTAT_FILE | head -n 1 -)

  unset COL_NUM_EXP_MOD
  declare -A COL_NUM_EXP_MOD
  count=0
  for i in $HEAD_SS; do
  count=$((count+1))
  COL_NUM_EXP_MOD[$i]=$count
  done

  SAMPLE_SIZE_EXPOSURE_MODEL=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$EXPOSURE_MODEL'"{print $('${COL_NUM_FOR["N_total"]}')}' $GRS_FORMATING_FILE)
  N_CASE_EXPOSURE_MODEL=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$EXPOSURE_MODEL'"{print $('${COL_NUM_FOR["N_case"]}')}' $GRS_FORMATING_FILE)
  N_CONTROL_EXPOSURE_MODEL=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$EXPOSURE_MODEL'"{print $('${COL_NUM_FOR["N_ctrl"]}')}' $GRS_FORMATING_FILE)
  
  ETH_EXPOSURE_MODEL=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$EXPOSURE_MODEL'"{print $('${COL_NUM_FOR["Ancestry"]}')}' $GRS_FORMATING_FILE)


  # USE GNOMAD EAF IF MISSING
  if [ $(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$EXPOSURE_MODEL'" {print $('${COL_NUM_FOR["EAF"]}')}' $GRS_FORMATING_FILE) == "NA" ] && [ $ETH_EXPOSURE_MODEL != "TRANS_ETHNIC" ]; then
    gzip -cd $SUMSTAT_FILE | awk 'BEGIN{FS=OFS="\t"}NR>1 && (length($('${COL_NUM_EXP_MOD["ref"]}'))+length($('${COL_NUM_EXP_MOD["alt"]}'))==2) && $('${COL_NUM_EXP_MOD["effect_size"]}')!=0{print "'$EXPOSURE_MODEL'",$('${COL_NUM_EXP_MOD["chr_hg19"]}')"_"$('${COL_NUM_EXP_MOD["pos_hg19"]}')"_"$('${COL_NUM_EXP_MOD["ref"]}')"_"$('${COL_NUM_EXP_MOD["alt"]}'),$('${COL_NUM_EXP_MOD["chr_hg19"]}'), $('${COL_NUM_EXP_MOD["pos_hg19"]}'),$('${COL_NUM_EXP_MOD["ref"]}'),$('${COL_NUM_EXP_MOD["alt"]}'),$('${COL_NUM_EXP_MOD["effect_size"]}'),$('${COL_NUM_EXP_MOD["standard_error"]}'),$('${COL_NUM_EXP_MOD["gnomAD_AF_"${ETH_EXPOSURE_MODEL}]}'),$('${COL_NUM_EXP_MOD["pvalue"]}'),"NA","NA","'$SAMPLE_SIZE_EXPOSURE_MODEL'", "'$N_CASE_EXPOSURE_MODEL'", "'$N_CONTROL_EXPOSURE_MODEL'"}' - | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$3,$4,$5,$6]==1{print}' list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} - > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp
  else 
    gzip -cd $SUMSTAT_FILE | awk 'BEGIN{FS=OFS="\t"}NR>1 && (length($('${COL_NUM_EXP_MOD["ref"]}'))+length($('${COL_NUM_EXP_MOD["alt"]}'))==2) && $('${COL_NUM_EXP_MOD["effect_size"]}')!=0{print "'$EXPOSURE_MODEL'",$('${COL_NUM_EXP_MOD["chr_hg19"]}')"_"$('${COL_NUM_EXP_MOD["pos_hg19"]}')"_"$('${COL_NUM_EXP_MOD["ref"]}')"_"$('${COL_NUM_EXP_MOD["alt"]}'),$('${COL_NUM_EXP_MOD["chr_hg19"]}'), $('${COL_NUM_EXP_MOD["pos_hg19"]}'),$('${COL_NUM_EXP_MOD["ref"]}'),$('${COL_NUM_EXP_MOD["alt"]}'),$('${COL_NUM_EXP_MOD["effect_size"]}'),$('${COL_NUM_EXP_MOD["standard_error"]}'),$('${COL_NUM_EXP_MOD["EAF"]}'),$('${COL_NUM_EXP_MOD["pvalue"]}'),"NA","NA","'$SAMPLE_SIZE_EXPOSURE_MODEL'", "'$N_CASE_EXPOSURE_MODEL'", "'$N_CONTROL_EXPOSURE_MODEL'"}' - | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$3,$4,$5,$6]==1{print}' list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} - > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp
  fi 

  # IF SOME MISSING EAF
  if [ $(awk 'BEGIN{FS=OFS="\t"}NR>1 && ($('${COL_NUM["eaf"]}') == NA || $('${COL_NUM["eaf"]}') == "NA" || $('${COL_NUM["eaf"]}') == ""){print $0}' ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp | wc -l) -gt 0 ] && [ $ETH_EXPOSURE_MODEL != "TRANS_ETHNIC" ];then
    gzip -cd $SUMSTAT_FILE | awk 'BEGIN{FS=OFS="\t"}NR>1{print $('${COL_NUM_EXP_MOD["chr_hg19"]}')"_"$('${COL_NUM_EXP_MOD["pos_hg19"]}')"_"$('${COL_NUM_EXP_MOD["ref"]}')"_"$('${COL_NUM_EXP_MOD["alt"]}'),$('${COL_NUM_EXP_MOD["gnomAD_AF_"${ETH_EXPOSURE_MODEL}]}')}' - | awk 'BEGIN{FS=OFS="\t"}NR==FNR{MAF[$1]=$5;next} MAF[$('${COL_NUM["SNP"]}')]!="" && MAF[$('${COL_NUM["SNP"]}')]!=NA && MAF[$('${COL_NUM["SNP"]}')]!="NA" && ($('${COL_NUM["eaf"]}')== "" ||  $('${COL_NUM["eaf"]}')=="NA" ||  $('${COL_NUM["eaf"]}') == NA) {$('${COL_NUM["eaf"]}')=MAF[$('${COL_NUM["SNP"]}')]; print $0}{print $0}' - ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp2
    mv ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp2 ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp
  fi

  # UPDATED LIST SNPS
  awk 'BEGIN{FS=OFS="\t"}NR==FNR{SNP[$1,$2,$3,$4]=1;SNP[$1,$2,$4,$3]=1;next}SNP[$('${COL_NUM["chr"]}'), $('${COL_NUM["pos"]}'), $('${COL_NUM["other_allele"]}'), $('${COL_NUM["effect_allele"]}')]==1{print $('${COL_NUM["chr"]}'), $('${COL_NUM["pos"]}'), $('${COL_NUM["other_allele"]}'), $('${COL_NUM["effect_allele"]}')}' list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp | sort -u > list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp
  mv list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}


  echo $(echo $(wc -l list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} | awk '{print $1}' - ) "SNPs IN COMMON")

  cat ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp2
  mv ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp2 ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
  rm -f ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp

  # UPDATE ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_
  awk 'BEGIN{FS=OFS="\t"}NR==FNR{SNP[$1"_"$2"_"$3"_"$4]=1;SNP[$1"_"$2"_"$4"_"$3]=1;next}FNR==1 || SNP[$('${COL_NUM["SNP"]}')]==1{print $0}' list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} > ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp
  mv ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_temp ${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}

done # END EXPOSURE_MODEL LOOP


####### OUTCOME #######
echo "" | awk 'BEGIN{FS=OFS="\t";print "Phenotype","SNP","chr","pos","beta","se","effect_allele","other_allele","eaf","pval","units","gene","samplesize","ncase","ncontrol"}' > $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}


if [ -f ${GRS_DIR}${OUTCOME}/4_final.txt.gz ]; then
    SUMSTAT_FILE=${GRS_DIR}${OUTCOME}/4_final.txt.gz
else
  echo $OUTCOME | cat - >> ${ROOT_OUTPUT_DIR}MISSING_OUTCOMES
  cd $ROOT_OUTPUT_DIR || exit 1
  rm -r $OUTPUT_DIR
  exit 1
fi

HEAD_SS=$(gzip -cd $SUMSTAT_FILE | head -n 1 -)

unset COL_NUM_OUT
declare -A COL_NUM_OUT
count=0
for i in $HEAD_SS; do
count=$((count+1))
COL_NUM_OUT[$i]=$count
done


SAMPLE_SIZE_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["N_total"]}')}' $GRS_FORMATING_FILE)
N_CASE_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["N_case"]}')}' $GRS_FORMATING_FILE)
N_CONTROL_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["N_ctrl"]}')}' $GRS_FORMATING_FILE)
ETH_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["Ancestry"]}')}' $GRS_FORMATING_FILE)

gzip -cd $SUMSTAT_FILE | awk 'BEGIN{FS=OFS="\t"}NR>1 && (length($('${COL_NUM_OUT["ref"]}'))+length($('${COL_NUM_OUT["alt"]}'))==2) && $('${COL_NUM_OUT["effect_size"]}')!=0{print "'$OUTCOME'",$('${COL_NUM_OUT["chr_hg19"]}')"_"$('${COL_NUM_OUT["pos_hg19"]}')"_"$('${COL_NUM_OUT["ref"]}')"_"$('${COL_NUM_OUT["alt"]}'),$('${COL_NUM_OUT["chr_hg19"]}'), $('${COL_NUM_OUT["pos_hg19"]}'), $('${COL_NUM_OUT["effect_size"]}'),$('${COL_NUM_OUT["standard_error"]}'),$('${COL_NUM_OUT["alt"]}'),$('${COL_NUM_OUT["ref"]}'),$('${COL_NUM_OUT["EAF"]}'),$('${COL_NUM_OUT["pvalue"]}'),"NA","NA","'$SAMPLE_SIZE_OUTCOME'", "'$N_CASE_OUTCOME'", "'$N_CONTROL_OUTCOME'"}' - | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$3,$4,$7,$8]==1{print}' list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} - >> $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}

echo "### N SNPs OUTCOME_MRBASE.FORMAT" $(awk 'BEGIN{FS=OFS="\t"}NR>1{print}' $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} | wc -l ) "###"

if [ $(awk 'BEGIN{FS=OFS="\t"}NR>1{print}' $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} | wc -l) -eq 0 ];then
  echo "NO OUTCOME_MRBASE.FORMAT OUTCOME" $OUTCOME "PANEL" $PANEL "CIS" ${CIS_WINDOW} "P-VALUE" ${PVALUE_THRESHOLD} "LD THRESHOLD" ${LD_THRESHOLD_PRUNING} > ${ROOT_OUTPUT_DIR}exit_${PANEL}_${OUTCOME}_${CIS_WINDOW}_${PVALUE_THRESHOLD}_${LD_THRESHOLD_PRUNING}
  rm -f list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL} \
  $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
  exit 1
fi

rm -f list_SNPs_exposures_PANEL_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}

mkdir -p 1_MATCH_TO_CONSORTIA

mkdir -p 2_LD_PRUNE

mkdir -p 3_EXTRACT_OUTCOME

mkdir -p 4_EXTRACT_EXPOSURES

mkdir -p 5_MR_RESULTS

unset COL_NUM_OUT
count=0
declare -A COL_NUM_OUT
for name in $(head -n 1 $BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}) ; do
count=$((count+1))
COL_NUM_OUT[$name]=$count
done

unset col_num
declare -A col_num
count=0
for i in $(head -n 1 $IVs_PURE); do
count=$((count+1))
col_num[$i]=$count
done

##### 1) MATCH EXPOSURES SNPS WITH OUTCOME SNPS #####
echo "#### 1) MATCH EXPOSURES SNPS WITH OUTCOME SNPS ####"
1_MATCH_TO_CONSORTIA_BIOMARKERS(){
  EXPOSURES=$1
  INPUT_FILE=$2
  OUTPUT_FILE=$3
  
  head -n 1 $EXPOSURES > $OUTPUT_FILE
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $('${COL_NUM_OUT["chr"]}'), $('${COL_NUM_OUT["pos"]}'), $('${COL_NUM_OUT["other_allele"]}'), $('${COL_NUM_OUT["effect_allele"]}')}' $INPUT_FILE | sort -u | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$('${col_num["chr"]}'),$('${col_num["pos"]}'),$('${col_num["other_allele"]}'),$('${col_num["effect_allele"]}')]==1 && $('${COL_NUM_PURE["Phenotype"]}')=="'$POPULATION'""_""'$PANEL'""_""'$ASSAY'""_""'$GENE'"{$('${col_num["SNP"]}')=$('${col_num["chr"]}')"_"$('${col_num["pos"]}')"_"$('${col_num["other_allele"]}')"_"$('${col_num["effect_allele"]}'); print $0}' - $EXPOSURES >> $OUTPUT_FILE
}

EXPOSURES=$IVs_PURE
INPUT_FILE=$BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
OUTPUT_FILE=${OUTPUT_DIR}1_MATCH_TO_CONSORTIA/EXPOSURE.MATCHED_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
1_MATCH_TO_CONSORTIA_BIOMARKERS $EXPOSURES $INPUT_FILE $OUTPUT_FILE


##### 2) PERFORM CLUMPING #####
echo "#### 2) PERFORM CLUMPING ####"
2_CLUMPING(){
  INPUT_FILE=$1
  LD_THRESHOLD_PRUNING=$2
  PVALUE_THRESHOLD=$3
  OUTPUT_FILE=$4
  ETH_OUTCOME=$5
  
  #### TEST ####
  #  INPUT_FILE=${OUTPUT_DIR}1_MATCH_TO_CONSORTIA/EXPOSURES.MATCHED_PANEL_${PANEL}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}
  #  OUTPUT_FILE=${OUTPUT_DIR}2_LD_PRUNE/EXPOSURES.LD_PRUNED_${LD_THRESHOLD_PRUNING}_PANEL_${PANEL}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}
  #  ETH_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["Ancestry"]}')}' $GRS_FORMATING_FILE)
  ##############

  rm -f ${OUTPUT_FILE}*
  
  head -n 1 $INPUT_FILE > $OUTPUT_FILE
  count=0

  # phenotype=$(awk 'BEGIN {FS=OFS="\t"}NR>1{print $('${col_num["Phenotype"]}')}' $INPUT_FILE | sort -u | head -n 1)
  # phenotype="PERSIAN_CVDII_PRSS2_PRSS2"
  # phenotype="AFRICAN_CVDII_ANG_ANG"
  for phenotype in $(awk 'BEGIN {FS=OFS="\t"}NR>1{print $('${col_num["Phenotype"]}')}' $INPUT_FILE | sort -u); do
    count=$((count+1))

    echo "########## $phenotype $count/$(awk 'BEGIN {FS=OFS="\t"}NR>1{print $('${col_num["Phenotype"]}')}' $INPUT_FILE | sort -u | wc -l) ##########"
    
    awk 'BEGIN {FS=OFS="\t"} NR>1 && $('${col_num["Phenotype"]}')=="'$phenotype'"{print}' $INPUT_FILE | sort -u > ${OUTPUT_FILE}_${phenotype}_temp_input.extract
    
    chr=$(awk '{print $('${col_num["chr"]}')}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | sort -u)
    
    POPULATION=$(echo $phenotype | sed 's/'_${PANEL}\.*'//g')
    

    if [ "$POPULATION" == "METAL_LATIN_EUROPEAN_PERSIAN" ]; then
      declare -a LIST_POP=("LATIN" "EUROPEAN" "PERSIAN")
    elif [ "$POPULATION" == "METAL_ALL_ETHNICITIES" ]; then
      declare -a LIST_POP=("LATIN" "EUROPEAN" "PERSIAN" "ARAB" "AFRICAN" "SOUTH_ASIAN" "EAST_ASIAN")
    else
      declare -a LIST_POP=($POPULATION)
    fi

    if [ "$ETH_OUTCOME" == "TRANS_ETHNIC" ]; then
      LIST_POP+=("LATIN" "EUROPEAN" "PERSIAN" "ARAB" "AFRICAN" "SOUTH_ASIAN" "EAST_ASIAN")
    elif [ "$ETH_OUTCOME" == "AFR" ]; then
      LIST_POP+=("AFRICAN")
    elif [ "$ETH_OUTCOME" == "EAS" ]; then
      LIST_POP+=("EAST_ASIAN")
    elif [ "$ETH_OUTCOME" == "SAS" ]; then
      LIST_POP+=("SOUTH_ASIAN")
    elif [ "$ETH_OUTCOME" == "EUR" ]; then
      LIST_POP+=("EUROPEAN")
    elif [ "$ETH_OUTCOME" == "AMR" ]; then
      LIST_POP+=("LATIN")
    fi

    # remove duplicate value in LIST_POP
    IFS=" " read -r -a LIST_POP <<< "$(echo "${LIST_POP[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"

    # KEEP SNPS IN ALL THE GENOTYPING FILES
    # TEST
    # POP="AFRICAN"
    for POP in ${LIST_POP[@]}; do
      echo "###########################################################"
      echo "########### KEEP SNPS IN GENOTYPING FILE" $POP "###########"
      echo "###########################################################"
      PURE_GENO_DIR=/genetics/PAREG/perrotn/PURE/IMPUTED_GENOTYPES_TOPMED_HG19/${POP}/
      PURE_input=${PURE_GENO_DIR}${chr}.IMP_QC_HG19

      # awk 'BEGIN{FS=OFS="\t"}{print "chr"$3":"$4":"$5":"$6;print "chr"$3":"$4":"$6":"$5}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | grep -wf - "${PURE_input}".bim | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$4,$5,$6]=1;snp[$1,$4,$6,$5]=1;next}snp[$3,$4,$5,$6]==1{print}' - ${OUTPUT_FILE}_${phenotype}_temp_input.extract > ${OUTPUT_FILE}_${phenotype}_temp_input.extract_temp
      awk 'BEGIN{FS=OFS="\t"}{print "chr"$3":"$4":"$5":"$6;print "chr"$3":"$4":"$6":"$5}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | awk 'BEGIN{FS=OFS="\t"}NR==FNR{SNP[$1]=1;next}SNP[$2]==1{print $0}' - "${PURE_input}".bim | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$4,$5,$6]=1;snp[$1,$4,$6,$5]=1;next}snp[$3,$4,$5,$6]==1{print}' - ${OUTPUT_FILE}_${phenotype}_temp_input.extract > ${OUTPUT_FILE}_${phenotype}_temp_input.extract_temp
      mv ${OUTPUT_FILE}_${phenotype}_temp_input.extract_temp ${OUTPUT_FILE}_${phenotype}_temp_input.extract

      # awk 'BEGIN{FS=OFS="\t"}{print $1, $4, $5, $6}' "${PURE_input}".bim | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$3,$4,$5,$6]==1{print}' - ${OUTPUT_FILE}_${phenotype}_temp_input.extract > ${OUTPUT_FILE}_${phenotype}_temp_input.extract_temp
      # mv ${OUTPUT_FILE}_${phenotype}_temp_input.extract_temp ${OUTPUT_FILE}_${phenotype}_temp_input.extract

    done # END POP LOOP

    if [ $(awk 'BEGIN{FS=OFS="\t"}{print}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | wc -l) -eq 0 ]; then
      echo "###################################"
      echo "######### NO SNPS REMAINS #########"
      echo "###################################"
      rm -f ${OUTPUT_FILE}_${phenotype}_temp_input.extract
      continue
    fi

    awk 'BEGIN {FS=OFS="\t"} ((length($('${col_num["other_allele"]}'))+length($('${col_num["effect_allele"]}'))==2) && $('${col_num["chr"]}')=="'$chr'") {print "chr"$('${col_num["chr"]}')":" $('${col_num["pos"]}') ":" $('${col_num["effect_allele"]}') ":" $('${col_num["other_allele"]}');print "chr" $('${col_num["chr"]}')":" $('${col_num["pos"]}') ":" $('${col_num["other_allele"]}') ":" $('${col_num["effect_allele"]}')}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | sort -gk 2 > ${OUTPUT_FILE}_${phenotype}_temp.ids
    awk 'BEGIN {FS=OFS="\t"} ((length($('${col_num["other_allele"]}'))+length($('${col_num["effect_allele"]}'))==2) && $('${col_num["chr"]}')=="'$chr'") {print "chr" $('${col_num["chr"]}') ":" $('${col_num["pos"]}') ":" $('${col_num["effect_allele"]}') ":" $('${col_num["other_allele"]}'), $0;print "chr" $('${col_num["chr"]}') ":" $('${col_num["pos"]}') ":" $('${col_num["other_allele"]}') ":" $('${col_num["effect_allele"]}'), $0}' ${OUTPUT_FILE}_${phenotype}_temp_input.extract | sort -gk 4,4 > ${OUTPUT_FILE}_${phenotype}_temp.input
    
    
    if [ $(wc -l ${OUTPUT_FILE}_${phenotype}_temp.ids | awk '{print $1}') -eq 2 ]; then # ONLY 1 SNP
    
      grep -wf ${OUTPUT_FILE}_${phenotype}_temp.ids ${OUTPUT_FILE}_${phenotype}_temp.input | head -n1 | cut -d$'\t' -f 2- >> ${OUTPUT_FILE}
    
    else

    # TEST
    # POP="EUROPEAN"
    for POP in ${LIST_POP[@]}; do
      echo "######################################"
      echo "########### --CLUMP" $POP "###########"
      echo "######################################"

      PURE_GENO_DIR=/genetics/PAREG/perrotn/PURE/IMPUTED_GENOTYPES_TOPMED_HG19/${POP}/
      PURE_input=${PURE_GENO_DIR}${chr}.IMP_QC_HG19

      echo "" | awk 'BEGIN{FS=OFS="\t"}{print "SNP", "P"}' - > ${OUTPUT_FILE}_${phenotype}_temp.assoc
      # awk 'BEGIN{FS=OFS="\t"}{print $1, $11}' ${OUTPUT_FILE}_${phenotype}_temp.input >> ${OUTPUT_FILE}_${phenotype}_temp.assoc
      cut -d$'\t' -f 2- ${OUTPUT_FILE}_${phenotype}_temp.input | awk 'BEGIN{FS=OFS="\t"}NR==FNR{SNP[$('${col_num["SNP"]}')]=1;PVAL[$('${col_num["SNP"]}')]=$('${col_num["pval"]}');next}SNP[$3]==1{print $1, PVAL[$3]}' - ${OUTPUT_FILE}_${phenotype}_temp.input >> ${OUTPUT_FILE}_${phenotype}_temp.assoc

      $plink2 \
      --bfile $PURE_input \
      --keep /genetics/PAREG/common/PURE_DATA/GENOTYPING/PURE_PMRA_COMBINED_2020_01_06.FINAL.fam \
      --extract ${OUTPUT_FILE}_${phenotype}_temp.ids \
      --clump ${OUTPUT_FILE}_${phenotype}_temp.assoc \
      --clump-p1 "$PVALUE_THRESHOLD" \
      --clump-p2 "$PVALUE_THRESHOLD" \
      --clump-r2 "$LD_THRESHOLD_PRUNING" \
      --clump-kb 1000 \
      --out ${OUTPUT_FILE}_${phenotype}_temp.clumped

      if [ -f ${OUTPUT_FILE}_${phenotype}_temp.clumped.clumped ]
      then
        # update ${OUTPUT_FILE}_${phenotype}_temp.ids, ${OUTPUT_FILE}_${phenotype}_temp.input
        awk 'BEGIN{FS=" ";OFS="\t"}NR>1 && $3!="" {print $3}' ${OUTPUT_FILE}_${phenotype}_temp.clumped.clumped | grep -wf - ${OUTPUT_FILE}_${phenotype}_temp.ids > ${OUTPUT_FILE}_${phenotype}_temp.ids_temp
        mv ${OUTPUT_FILE}_${phenotype}_temp.ids_temp ${OUTPUT_FILE}_${phenotype}_temp.ids

        awk 'BEGIN{FS=" ";OFS="\t"}NR>1 && $3!="" {print $3}' ${OUTPUT_FILE}_${phenotype}_temp.clumped.clumped | grep -wf - ${OUTPUT_FILE}_${phenotype}_temp.input > ${OUTPUT_FILE}_${phenotype}_temp.input_temp
        mv ${OUTPUT_FILE}_${phenotype}_temp.input_temp ${OUTPUT_FILE}_${phenotype}_temp.input

        if [ $(wc -l ${OUTPUT_FILE}_${phenotype}_temp.ids | awk '{print $1}') -eq 1 ]; then # ONLY 1 SNP
          break
        fi
      else
        > ${OUTPUT_FILE}_${phenotype}_temp.ids
        break
      fi
    done # END POP LOOP

      if [ $(wc -l ${OUTPUT_FILE}_${phenotype}_temp.ids | awk '{print $1}' - ) -eq 1 ]; then # ONLY 1 SNP
        grep -wf ${OUTPUT_FILE}_${phenotype}_temp.ids ${OUTPUT_FILE}_${phenotype}_temp.input | head -n1 | cut -d$'\t' -f 2- >> ${OUTPUT_FILE}
        rm ${OUTPUT_FILE}_${phenotype}_temp*
        continue
      fi
      
      if [ $(wc -l ${OUTPUT_FILE}_${phenotype}_temp.ids | awk '{print $1}' - ) -eq 0 ]; then
        echo "NO SNPS REMAINS"
        rm ${OUTPUT_FILE}_${phenotype}_temp*
        continue
      fi

      # grep -wf ${OUTPUT_FILE}_${phenotype}_temp.ids ${OUTPUT_FILE}_${phenotype}_temp.input | cut -d$'\t' -f 2- >> ${OUTPUT_FILE}
      # rm ${OUTPUT_FILE}_${phenotype}_temp*

      rm -f ${OUTPUT_FILE}_${phenotype}_temp.r2
      # # TEST
      # # POP="AFRICAN"
      for POP in ${LIST_POP[@]}; do
        echo "###################################"
        echo "########### --R2" $POP "###########"
        echo "###################################" 
        PURE_GENO_DIR=/genetics/PAREG/perrotn/PURE/IMPUTED_GENOTYPES_TOPMED_HG19/${POP}/
        PURE_input=${PURE_GENO_DIR}${chr}.IMP_QC_HG19
        
        # remove existing file
        rm -f ${OUTPUT_FILE}_${phenotype}_temp.ids.*

        # generate genotyping file with --keep
        $plink2 --bfile $PURE_input \
        --keep /genetics/PAREG/common/PURE_DATA/GENOTYPING/PURE_PMRA_COMBINED_2020_01_06.FINAL.fam \
        --extract ${OUTPUT_FILE}_${phenotype}_temp.ids \
        --make-bed \
        --out ${OUTPUT_FILE}_${phenotype}_temp.ids

        # generate r2 matrix
        $plink2 \
        --bfile ${OUTPUT_FILE}_${phenotype}_temp.ids \
        --r2 yes-really inter-chr \
        --ld-window-r2 0 \
        --out ${OUTPUT_FILE}_${phenotype}_temp.r2

        if [ ! -s ${OUTPUT_FILE}_${phenotype}_temp.r2.ld ]; then
          echo "ERROR" ${OUTPUT_FILE}_${phenotype}_temp.r2.ld "DOES NO EXISTS"
          echo "phenotype:" $phenotype
          echo "POP:" $POP
          exit 1
        fi

        awk 'BEGIN{OFS="\t"}NR>1{print}' ${OUTPUT_FILE}_${phenotype}_temp.r2.ld >> ${OUTPUT_FILE}_${phenotype}_temp.r2
        
      done # END POP LOOP
      
      awk 'BEGIN{OFS="\t"}NR>1{print $3; print $6}' ${OUTPUT_FILE}_${phenotype}_temp.r2 | sort -u | grep -wf - ${OUTPUT_FILE}_${phenotype}_temp.input > ${OUTPUT_FILE}_${phenotype}_temp
      mv ${OUTPUT_FILE}_${phenotype}_temp ${OUTPUT_FILE}_${phenotype}_temp.input
      
      
      until [[ ! -s "${OUTPUT_FILE}_${phenotype}_temp.input" ]] ; do
      var=$(awk 'NR==1{print $1}'  ${OUTPUT_FILE}_${phenotype}_temp.input)
      # Find proxy SNPs based on r2 LD_THRESHOLD_PRUNING
      grep -w $var ${OUTPUT_FILE}_${phenotype}_temp.r2 | awk 'BEGIN {OFS="\t"}$3=="'$var'"{if ($NF>='${LD_THRESHOLD_PRUNING}')print $6}$6=="'$var'"{if ($NF>='${LD_THRESHOLD_PRUNING}')print $3}' | sort -u | grep -w -v $var > ${OUTPUT_FILE}_${phenotype}_temp_input.proxy_SNPS
      grep -w -v -f ${OUTPUT_FILE}_${phenotype}_temp_input.proxy_SNPS ${OUTPUT_FILE}_${phenotype}_temp.r2 > ${OUTPUT_FILE}_${phenotype}_temp
      mv ${OUTPUT_FILE}_${phenotype}_temp ${OUTPUT_FILE}_${phenotype}_temp.r2
      # Output query variant to final file
      awk 'BEGIN{FS=OFS="\t"}{print}' ${OUTPUT_FILE}_${phenotype}_temp.input | grep -w $var - | cut -d$'\t' -f 2- >> ${OUTPUT_FILE}
      # Remove query variant from SNP info file
      # awk 'NR>1' ${OUTPUT_FILE}_${phenotype}_temp.input > ${OUTPUT_FILE}_${phenotype}_temp # Remove first variant (query)=
      grep -w -v $var ${OUTPUT_FILE}_${phenotype}_temp.input > ${OUTPUT_FILE}_${phenotype}_temp
      # Remove proxy SNPs from SNP info file
      grep -w -v -f ${OUTPUT_FILE}_${phenotype}_temp_input.proxy_SNPS ${OUTPUT_FILE}_${phenotype}_temp > ${OUTPUT_FILE}_${phenotype}_temp.input # Remove variants that are in LD.
      done

    fi # END IF 1 SNP
    
    rm -f ${OUTPUT_FILE}_${phenotype}_temp*
      
  done # phenotype loop
}


INPUT_FILE=${OUTPUT_DIR}1_MATCH_TO_CONSORTIA/EXPOSURE.MATCHED_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
OUTPUT_FILE=${OUTPUT_DIR}2_LD_PRUNE/EXPOSURE.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
ETH_OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1 && $('${COL_NUM_FOR["trait_name"]}')=="'$OUTCOME'"{print $('${COL_NUM_FOR["Ancestry"]}')}' $GRS_FORMATING_FILE)
2_CLUMPING $INPUT_FILE $LD_THRESHOLD_PRUNING $PVALUE_THRESHOLD $OUTPUT_FILE $ETH_OUTCOME


##### 3) EXTRACT OUTCOME SUMMARY STATISTICS #####
echo "#### 3) EXTRACT OUTCOME SUMMARY STATISTICS ####"
3_EXTRACT_OUTCOME(){
  EXPOSURE=$1
  OUTCOME=$2
  OUTPUT_FILE=$3
  head -n 1 $OUTCOME > $OUTPUT_FILE
  
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $('${col_num["chr"]}'), $('${col_num["pos"]}'), $('${col_num["other_allele"]}'), $('${col_num["effect_allele"]}')}' $EXPOSURE | sort -u | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$('${COL_NUM_OUT["chr"]}'),$('${COL_NUM_OUT["pos"]}'),$('${COL_NUM_OUT["other_allele"]}'),$('${COL_NUM_OUT["effect_allele"]}')]==1' - $OUTCOME >> $OUTPUT_FILE
}

EXPOSURE=${OUTPUT_DIR}2_LD_PRUNE/EXPOSURE.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
OUTCOME=$BIOMARKER_OUTCOME_DIR/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
OUTPUT_FILE=${OUTPUT_DIR}3_EXTRACT_OUTCOME/OUTCOME.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
3_EXTRACT_OUTCOME $EXPOSURE $OUTCOME $OUTPUT_FILE


##### 4) EXTRACT EXPOSURES #####
echo "#### 4) EXTRACT EXPOSURES ####"
4_EXTRACT_EXPOSURES(){
  EXPOSURE=$1
  EXPOSURES=$2
  OUTPUT_FILE=$3
  head -n1 $EXPOSURE | awk 'BEGIN{FS=OFS="\t"}{print $0, "ncase","ncontrol"}' - > $OUTPUT_FILE
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $0, NA, NA}' $EXPOSURE >> $OUTPUT_FILE
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $('${col_num["chr"]}'), $('${col_num["pos"]}'), $('${col_num["other_allele"]}'), $('${col_num["effect_allele"]}')}' $EXPOSURE | sort -u | awk 'BEGIN{FS=OFS="\t"} NR==FNR {snp[$1,$2,$3,$4]=1;snp[$1,$2,$4,$3]=1;next}snp[$('${col_num["chr"]}'),$('${col_num["pos"]}'),$('${col_num["other_allele"]}'),$('${col_num["effect_allele"]}')]==1' - $EXPOSURES >> $OUTPUT_FILE
}

EXPOSURE=${OUTPUT_DIR}2_LD_PRUNE/EXPOSURE.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
EXPOSURES=${BIOMARKER_EXPOSURES_DIR}/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
OUTPUT_FILE=${OUTPUT_DIR}4_EXTRACT_EXPOSURES/EXPOSURES.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}
4_EXTRACT_EXPOSURES $EXPOSURE $EXPOSURES $OUTPUT_FILE


######################################################################################################################################################################################
SCRIPT=$1
CIS_WINDOW=$2
LD_THRESHOLD_PRUNING=$3
PANEL=$4
ASSAY=$5
GENE=$6
ROOT_DIR=$7
GRS_FORMATING_FILE=$8
OUTCOME=$9
LIST_EXPOSURE_MODEL=${10}
MODEL=${11}

echo "SCRIPT:" $SCRIPT
echo "CIS_WINDOW:" $CIS_WINDOW
echo "LD_THRESHOLD_PRUNING:" $LD_THRESHOLD_PRUNING
echo "PANEL:" $PANEL
echo "ASSAY:" $ASSAY
echo "GENE:" $GENE
echo "ROOT_DIR:" $ROOT_DIR
echo "GRS_FORMATING_FILE:" $GRS_FORMATING_FILE
echo "OUTCOME:" $OUTCOME
echo "LIST_EXPOSURE_MODEL:" $LIST_EXPOSURE_MODEL

# rm -f ALREADY EXISTS

# python /genetics/PAREG/perrotn/scripts/combinations_MVMR.py --file $LIST_EXPOSURE_MODEL

export TMPDIR=${ROOT_DIR}tmp

RUN_PARALLEL () {
  if [[ $# -eq 0 ]] ; then
          echo '
  DESCRIPTION: THIS FUNCTION WILL RUN JOBS IN PARALLEL CONSIDERING REALTIME SERVER MEMORY AND CPU USAGE

  USAGE: RUN_PARALLEL <MAX_CPU_USAGE> <MAX_MEM_USAGE> <TIME_BETWEEN_CHECKS> <YOUR COMMAND OR SCRIPT HERE>
  *** IMPORTANT *** DO NOT ADD AN AMPERSAND "&" TO THE END OF YOUR COMMAND - THIS SCRIPT DOES IT AUTOMATICALLY

  DEFINITIONS:

  <MAX_CPU_USAGE> - The maximum fraction (0.00 - 1.00) of realtime cpu usage permitted.
  i.e. the script will proceed only when realtime cpu_usage < MAX_CPU_USAGE.
  WARNING: THIS SETTING SHOULD NEVER BE SET ABOVE 0.95.

  <MAX_MEM_USAGE> - The maximum fraction (0.00 - 1.00) of realtime memory usage permitted.
  i.e. the script will proceed only when realtime cpu_usage < MAX_CPU_USAGE.
  WARNING: THIS SETTING SHOULD NEVER BE SET ABOVE 0.95.

  <MAX_PARALLEL> - Maximum number of scripts to run at one time

  <TIME_BETWEEN_CHECKS> - How often the script should check to run the next job (in seconds).
  WARNING: THIS SETTING SHOULD NEVER BE SET BELOW 5 SECONDS.

  <NAME OF SCRIPT LOG>

  <YOUR COMMAND OR SCRIPT HERE>

  Example:
  script_log=$(mktemp script_log.XXXXXX)
  for {i in 1..100} ; do
  RUN_PARALLEL 0.95 0.95 2 5 $script_log DEPTH_OF_COVERAGE INPUT_BAM.bam ALL_GENES.bed
  done'
  return 1
  fi
  max_mem=$1
  max_cpu=$2
  max_parallel=$3
  seconds_bn_checks=$4
  script_log=$5
  shift 5
  COMMAND_TO_RUN=$@
    #### MINIMUM INTERVAL, MEMORY, AND CPU USAGE
    if [[ $(echo "" | awk ''$max_mem'>0.95{print "true"}') ]] ; then
  echo "MAXIMUM MEMORY LIMIT IS SET ABOVE 0.95. PLEASE SET TO A LOWER THRESHOLD TO BE SAFE."
  return 1
  fi
  if [[ $(echo "" | awk ''$max_cpu'>0.95{print "true"}') ]] ; then
  echo "MAXIMUM CPU USAGE IS SET ABOVE 0.95. PLEASE SET TO A HIGHER THRESHOLD TO BE SAFE."
  return 1
  fi
  if [[ $(echo "" | awk ''$seconds_bn_checks'<5{print "true"}') ]] ; then
  echo "TIME INTERVAL TO CHECK CPU/MEM USAGE IS SET BELOW 5 SECONDS. PLEASE SET TO A HIGHER THRESHOLD"
  return 1
  fi
  # MAKE SURE THAT MINIMUM THRESHOLDS ARE NOT SET TOO LOW AS TO CRASH THE SERVER
  while sleep ${seconds_bn_checks} ; do # Every X seconds check the fraction of free memory and cpu. If <= 0.10, don't do anything
  # Calculate memory and cpu free
  local mem_temp=$(mktemp mem.XXXXXX)
  local cpu_temp=$(mktemp cpu.XXXXXX)
  # Take an average of 50 snapshots - should take about 10 seconds to run
  for count in {1..50} ; do
  top -b -n 1 | awk 'NR>7 && NR==FNR{mem+=$10;cpu+=$9;PIDS[$1]=1;next} END{print mem/100 >> "'$mem_temp'" ; print cpu/('$(nproc)'*100) >> "'$cpu_temp'"}'
  done
  proportion_mem_used=$(awk '{mem+=$1}END{print mem/NR}' $mem_temp) # Server Memory adds up to 100 so all we have to do is divide by 100
  proportion_cpu_used=$(awk '{cpu+=$1}END{print cpu/NR}' $cpu_temp) # CPU usage is given per core so we have to divide by num_cores * 100
  rm -f $mem_temp
  rm -f $cpu_temp

  if [[ -s "$script_log" ]] ; then 
  top -b -n 1 | awk 'NR>7 && NR==FNR{PIDS[$1]=1;next} PIDS[$1]==1{print $1} ' - $script_log > $script_log.temp
  num_scripts_running=$(wc -l $script_log | awk '{print $1}') # Number of scripts running 
  mv $script_log.temp $script_log
  else 
    num_scripts_running=0
  fi


  # Output current memory and CPU usage
  echo "% MEMORY USED: $proportion_mem_used
  % CPU USED: $proportion_cpu_used
  # SCRIPTS RUNNING: $num_scripts_running
  "
  # Run one iteration if there is enough free memory
  if [[ $(echo $proportion_cpu_used $proportion_mem_used $num_scripts_running | awk '$1<'${max_cpu}' && $2 <'${max_mem}' && $3 <'${max_parallel}'{print "true"}') ]] ; then
  $COMMAND_TO_RUN & # ./CLINICAL_EXOME_COVERAGE.sh SAMPLE_NAME INPUT_BAM CANDIDATE_GENE INCIDENTAL_GENE INVESTIGATOR_GENE CLINVAR_VCF & # RUN COMMAND OR SCRIPT IN BACKGROUND HERE &
    PID=$!
    echo $PID >> $script_log
  break # Need to break out of the loop or else this command will just keep repeating every X seconds.
  fi
  done
} # END RUN_PARALLEL FUNCTION

cd ${ROOT_DIR}PURE/MVMR_consortia/

script_log=$(mktemp script_log.XXXXXX)

N_CPU=$(nproc)

for SUFFIX in NO_MHC_MISSENSE_SPLICING # NO_MHC NO_MHC_MISSENSE NO_MHC_SPLICING
do

  ROOT_OUTPUT_DIR=${ROOT_DIR}PURE/MVMR_consortia/${SUFFIX}/
  mkdir -p $ROOT_OUTPUT_DIR
  OUTPUT_DIR=${ROOT_OUTPUT_DIR}${OUTCOME}/
  mkdir -p $OUTPUT_DIR

  for PVALUE_THRESHOLD in 0.01 # 0.000005  0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.00000005
  do
    for POPULATION in EUROPEAN # METAL_LATIN_EUROPEAN_PERSIAN
    do
          rm -f ${OUTPUT_DIR}OUTCOME_FORMATTING_FOR_MRBASE/OUTCOME_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}* \
          ${OUTPUT_DIR}EXPOSURES_FORMATTING_FOR_MRBASE/EXPOSURES_MODEL_MRBASE.FORMAT_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}* \
          ${OUTPUT_DIR}1_MATCH_TO_CONSORTIA/EXPOSURE.MATCHED_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_LD_PRUNED_${LD_THRESHOLD_PRUNING}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}* \
          ${OUTPUT_DIR}2_LD_PRUNE/EXPOSURE.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}* \
          ${OUTPUT_DIR}3_EXTRACT_OUTCOME/OUTCOME.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}* \
          ${OUTPUT_DIR}4_EXTRACT_EXPOSURES/EXPOSURES.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}*
          
          RUN_PARALLEL 0.95 0.95 ${N_CPU} 8 "$script_log" "$SCRIPT" $SUFFIX $PANEL $ASSAY $GENE $LIST_EXPOSURE_MODEL $PVALUE_THRESHOLD $CIS_WINDOW $LD_THRESHOLD_PRUNING $OUTCOME $ROOT_DIR $GRS_FORMATING_FILE $MODEL $POPULATION
    done # END POPULATION LOOP
  done # END PVALUE_THRESHOLD LOOP
done # END SUFFIX LOOP

wait 

echo "DONE"

########################################################################################################################################
ROOT_DIR="/storage/genetics_work2/perrotn/"

mkdir -p ${ROOT_DIR}tmp
export TMPDIR=${ROOT_DIR}tmp

DIR="${ROOT_DIR}PURE/MVMR_consortia/"
mkdir -p $DIR
cd $DIR || exit 1

# SNAPSHOT GRS FORMATTING FILE
GRS_FORMATING_FILE=$(mktemp GRS_FORMATTING_FILE.XXXXXX)
Rscript /genetics/PAREG/perrotn/scripts/snapshot_GRS_Formatting.r "$GRS_FORMATING_FILE"
awk 'BEGIN{FS=OFS="\t"}{print}' "$GRS_FORMATING_FILE" | tr -d '\r' > "${GRS_FORMATING_FILE}"temp
mv "${GRS_FORMATING_FILE}"temp "$GRS_FORMATING_FILE"


SCRIPT="/genetics/PAREG/perrotn/scripts/IVs_MVMR_PURE.sh"
chmod u+x $SCRIPT

CIS_WINDOW=200000
LD_THRESHOLD_PRUNING=0.1
PANEL="CMET"
ASSAY="F11"
GENE="F11"
ROOT_DIR="/storage/genetics_work2/perrotn/"
# OUTCOME="ischemic_stroke_european_GIGASTROKE_2022"
OUTCOME="all_stroke_european_GIGASTROKE_2022"
# LIST_EXPOSURE_MODEL="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/FXI_CES_VS_ISCHEMIC_STROKE.txt" # HEADER: "exposures"
# MODEL="FXI_CES_VS_ISCHEMIC_STROKE"

LIST_EXPOSURE_MODEL="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/FXI_ISCHEMIC_STROKE_VS_ALL_STROKE.txt" # HEADER: "exposures"
MODEL="FXI_ISCHEMIC_STROKE_VS_ALL_STROKE"

# LIST_EXPOSURE_MODEL="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/test_FXI_ISCHEMIC_STROKE_VS_ISCHEMIC_STROKE.txt" # HEADER: "exposures"
# MODEL="test_FXI_ISCHEMIC_STROKE_VS_ISCHEMIC_STROKE"

awk 'BEGIN{FS=OFS="\t"}{print}' "$LIST_EXPOSURE_MODEL" | tr -d '\r' > "${LIST_EXPOSURE_MODEL}"temp
mv "${LIST_EXPOSURE_MODEL}"temp "$LIST_EXPOSURE_MODEL"

# rm $(echo $(echo $LIST_EXPOSURE_MODEL | sed 's/'.txt'/_/')[0-9]*)

nohup bash /genetics/PAREG/perrotn/scripts/IVs_MVMR_PURE_parallel.sh \
$SCRIPT \
$CIS_WINDOW \
$LD_THRESHOLD_PRUNING \
$PANEL \
$ASSAY \
$GENE \
$ROOT_DIR \
${ROOT_DIR}PURE/MVMR_consortia/${GRS_FORMATING_FILE} \
$OUTCOME \
$LIST_EXPOSURE_MODEL \
$MODEL \
> /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out &
########################################################################################################################################
# in progress xxxxx

# kill $(ps -aux | grep "IVs_MVMR_PURE_parallel.sh" | awk '{print $2}' -)
# kill $(ps -aux | grep "IVs_MVMR_PURE.sh" | awk '{print $2}' -)

/storage/genetics_work2/perrotn/PURE/MVMR_consortia/NO_MHC_MISSENSE_SPLICING/all_stroke_european_GIGASTROKE_2022/

grep -ci "error" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out
grep -i "fail" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out |  grep -cv "after larger attempt(s) failed" - 
grep -ci "cannot" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out
grep -ci "syntax error" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out

grep -in "error" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out | head
grep -i "fail" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out |  grep -nv "after larger attempt(s) failed" - 
grep -ni "cannot" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out | head
grep -in "syntax error" /genetics/PAREG/perrotn/scripts/output_nohup/IVs_MVMR_PURE_parallel.out | head


