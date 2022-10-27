#!/bin/bash

# /genetics/PAREG/perrotn/scripts/Bind_MVMR_BM_F11_PURE.sh

##########################################
##### MVMR_BM WITH IVs PURE x CONSORTIA #####
##########################################

echo "Started MVMR_BM at: "$(date)
SECONDS=0

# export TMPDIR=/genetics_work/perrotn/temp
export TMPDIR=/torrent_vol/perrotn/tmp

SUFFIX=$1
PANEL=$2
ASSAY=$3
GENE=$4
PVALUE_THRESHOLD=$5
CIS_WINDOW=$6
LD_THRESHOLD_PRUNING=$7
ROOT_DIR=$8
POPULATION=$9
LIST_EXPOSURES=${10}
LIST_OUTCOMES=${11}

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
# OUTCOME="all_stroke_european_GIGASTROKE_2022"

echo "SUFFIX:" $SUFFIX
echo "PANEL:" $PANEL
echo "ASSAY:" $ASSAY
echo "GENE:" $GENE
echo "PVALUE_THRESHOLD:" $PVALUE_THRESHOLD
echo "CIS_WINDOW:" $CIS_WINDOW
echo "LD_THRESHOLD_PRUNING:" $LD_THRESHOLD_PRUNING
echo "ROOT_DIR:" $ROOT_DIR
echo "POPULATION:" $POPULATION

ROOT_OUTPUT_DIR=${ROOT_DIR}PURE/MVMR_BM_consortia/${SUFFIX}/


OUTPUT_DIR=${ROOT_OUTPUT_DIR}/ALL_MVMR_BM_F11

mkdir -p $OUTPUT_DIR

cd ${OUTPUT_DIR} || exit 1

ALL_RESULTS_FILE=${OUTPUT_DIR}/MR_ALL_RESULTS.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}.txt

rm -f $ALL_RESULTS_FILE

# TEST
# OUTCOME=$(awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $LIST_OUTCOMES | head -n 1 )
# EXPOSURE=$(awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $LIST_EXPOSURES | head -n 1)
for OUTCOME in $(awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $LIST_OUTCOMES)
do

#   echo "OUTCOME:" $OUTCOME
  DIRECTORY=${ROOT_OUTPUT_DIR}${OUTCOME}/5_MR_RESULTS/LD_${LD_THRESHOLD_PRUNING}

  for EXPOSURE in $(awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $LIST_EXPOSURES)
  do

    # echo "EXPOSURE:" $EXPOSURE
    MODEL=$(echo $ASSAY"_"$EXPOSURE"_vs_"$OUTCOME)
    # echo "MODEL:" $MODEL

    if [ ! -f ${DIRECTORY}/MR_RESULTS.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}.txt ]
    then
        echo "NO MR FILE FOR OUTCOME:" $OUTCOME "& EXPOSURE:" $EXPOSURE
        continue
    fi

    # ADD TO MR RESULTS:
    # KEEP F-stat >=10
    # PRIO MR METHOD (IVW OR EGGER)

    Rscript /genetics/PAREG/perrotn/automated_scripts/ADD_FEATURES_MVMR_RESULTS_FILE.r --file ${DIRECTORY}/MR_RESULTS.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}.txt


    if [ ! -f $ALL_RESULTS_FILE ]
    then
        head -n 1 ${DIRECTORY}/MR_RESULTS.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_PRIO_MR_METHOD.txt > $ALL_RESULTS_FILE
    fi
    
    awk 'BEGIN{FS=OFS="\t"}NR>1{print}' ${DIRECTORY}/MR_RESULTS.LD_PRUNED_${LD_THRESHOLD_PRUNING}_${POPULATION}_${PANEL}_${ASSAY}_${GENE}_CIS_${CIS_WINDOW}_PVALUE_${PVALUE_THRESHOLD}_${MODEL}_PRIO_MR_METHOD.txt >> $ALL_RESULTS_FILE

  done # END EXPOSURE LOOP
done # END OUTCOME LOOP


duration=$SECONDS
echo "Total elapsed time: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) minutes $((duration % 60)) seconds."


######################################################################################################################################################################################
SCRIPT=$1
CIS_WINDOW=$2
LD_THRESHOLD_PRUNING=$3
PANEL=$4
ASSAY=$5
GENE=$6
ROOT_DIR=$7
LIST_EXPOSURES=$8
LIST_OUTCOMES=$9

echo "SCRIPT:" $SCRIPT
echo "CIS_WINDOW:" $CIS_WINDOW
echo "LD_THRESHOLD_PRUNING:" $LD_THRESHOLD_PRUNING
echo "PANEL:" $PANEL
echo "ASSAY:" $ASSAY
echo "GENE:" $GENE
echo "ROOT_DIR:" $ROOT_DIR

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


cd ${ROOT_DIR}PURE/MVMR_BM_consortia/

script_log=$(mktemp script_log.XXXXXX)

N_CPU=$(nproc)

for SUFFIX in NO_MHC_MISSENSE_SPLICING # NO_MHC NO_MHC_MISSENSE NO_MHC_SPLICING
do
    for PVALUE_THRESHOLD in 0.01 # 0.000005  0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.00000005
    do
        for POPULATION in EUROPEAN # METAL_LATIN_EUROPEAN_PERSIAN
        do
            RUN_PARALLEL 0.95 0.95 ${N_CPU} 6 "$script_log" "$SCRIPT" $SUFFIX $PANEL $ASSAY $GENE $PVALUE_THRESHOLD $CIS_WINDOW $LD_THRESHOLD_PRUNING $ROOT_DIR $POPULATION $LIST_EXPOSURES $LIST_OUTCOMES
        done # END POPULATION LOOP
    done # END PVALUE_THRESHOLD LOOP
done # END SUFFIX LOOP

wait 

echo "DONE"

########################################################################################################################################
ROOT_DIR="/storage/genetics_work2/perrotn/"

mkdir -p ${ROOT_DIR}tmp
export TMPDIR=${ROOT_DIR}tmp

DIR="${ROOT_DIR}PURE/MVMR_BM_consortia/"
mkdir -p $DIR
cd $DIR || exit 1

SCRIPT="/genetics/PAREG/perrotn/scripts/Bind_MVMR_BM_F11_PURE.sh"
chmod u+x $SCRIPT

CIS_WINDOW=200000
LD_THRESHOLD_PRUNING=0.1
PANEL="CMET"
ASSAY="F11"
GENE="F11"
ROOT_DIR="/storage/genetics_work2/perrotn/"

LIST_EXPOSURES=/storage/genetics_work2/perrotn/PURE/BM_BM_MR_V3/LIST_EXPOSURES_BM_BM.txt
LIST_OUTCOMES=/genetics_work/perrotn/PURE/pPheWMR_consortia_V5/OUTCOMES_LIST_F11_BM_BM_MR.txt

nohup bash /genetics/PAREG/perrotn/scripts/Bind_MVMR_BM_F11_PURE_parallel.sh \
$SCRIPT \
$CIS_WINDOW \
$LD_THRESHOLD_PRUNING \
$PANEL \
$ASSAY \
$GENE \
$ROOT_DIR \
$LIST_EXPOSURES \
$LIST_OUTCOMES \
> /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out &
########################################################################################################################################
# in progress xxxxxx

# kill $(ps -aux | grep "Bind_MVMR_BM_F11_PURE_parallel.sh" | awk '{print $2}' -)
# kill $(ps -aux | grep "Bind_MVMR_BM_F11_PURE.sh" | awk '{print $2}' -)


/storage/genetics_work2/perrotn/PURE/MVMR_BM_consortia/NO_MHC_MISSENSE_SPLICING/ischemic_stroke_european_GIGASTROKE_2022/5_MR_RESULTS/LD_0.1
/storage/genetics_work2/perrotn/PURE/MVMR_BM_consortia/NO_MHC_MISSENSE_SPLICING/all_stroke_european_GIGASTROKE_2022/5_MR_RESULTS/LD_0.1

cd /storage/genetics_work2/perrotn/PURE/MVMR_BM_consortia/NO_MHC_MISSENSE_SPLICING/ALL_MVMR_BM_F11

grep -c "NO MR FILE FOR OUTCOME:" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out

grep -ci "error" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out
grep -i "fail" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out |  grep -cv "after larger attempt(s) failed" - 
grep -ci "cannot" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out
grep -ci "syntax error" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out

grep -in "error" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out | head
grep -i "fail" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out |  grep -nv "after larger attempt(s) failed" - 
grep -ni "cannot" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out | head
grep -in "syntax error" /genetics/PAREG/perrotn/scripts/output_nohup/Bind_MVMR_BM_F11_PURE_parallel.out | head






