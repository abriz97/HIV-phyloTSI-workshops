#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI.sh
# If unset default to "sim"
if [[ "$1" == "-h" || "$1" == "-help" || "$1" == "--help" ]]; then
    echo "Usage: qsub STEP=step RES=res runall_TSI.sh"
    echo "OPTIONs:"
    echo "    STEP : one of net, ali, btr, ctr, atr, tsi, dti "
    echo "    RES: determines resources fr pbs jobs (0 to 3)  [default: 0]"
    exit 0
fi


if [ -z "$STEP" ]
then
        echo "Intended use:\n"
        echo 'qsub -v STEP="xxx" runall_TSI.sh'
        exit 1
fi

${RES:=1}
${REDO:=0}
echo "running '${STEP:=sim}' analysis"

# This includes all code necessary to run PHSC pipeline to produce TSI estimates
DEEPDATA="/rds/general/project/ratmann_pangea_deepsequencedata/live"
DEEPANALYSES="/rds/general/project/ratmann_deepseq_analyses/live"
HOME="/rds/general/user/ab1820/home"

software_path="$HOME/git/HIV-phyloTSI-workshops/src/pipeline"
phyloscanner_path="$HOME/git/phyloscanner"
hivtsipath="$HOME/git/HIV-phyloTSI"

# analysis specific paths & args
out_dir_base="$DEEPANALYSES/PANGEA2_KenyaWorkshop"
out_dir_rel="$out_dir_base/TODO"

controller="$software_path/runall_TSI.sh"
inputsamples="$out_dir_base/phsc_input_samples_with_bf.rds"
CLUSIZE='50'
DATE='2023-08-22'

echo "Check that DATE, CLUSIZE, out_dir_rel and inputsamples are correctly specified"

cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

        net)
        echo "---- initialise analysis ----"
        Rscript $software_path/01_TSI_initialise.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --cluster_size $CLUSIZE \
        --controller $controller
        ;;

        ali)
        echo "---- compute alignments ----"
        if [ "$REDO" = "0" ]; then
            Rscript $software_path/02_make_deep_sequence_alignments.R \
            --out_dir_base $out_dir_base \
            --pkg_dir $software_path \
            --prog_dir $phyloscanner_path \
            --windows_start 550 \
            --windows_end 9500 \
            --sliding_width 25 \
            --n_control 0 \
            --cluster_size $CLUSIZE \
            --reference 220419_reference_set_for_PARTNERS_mafft.fasta \
            --mafft " --globalpair --maxiterate 1000 " \
            --rm_vloops FALSE \
            --controller $controller \
            --walltime_idx $RES \
            --tsi_analysis FALSE
        else 
            Rscript $software_path/02_make_deep_sequence_alignments.R \
            --out_dir_base $out_dir_base \
            --pkg_dir $software_path \
            --prog_dir $phyloscanner_path \
            --windows_start 550 \
            --windows_end 9500 \
            --sliding_width 25 \
            --n_control 0 \
            --cluster_size $CLUSIZE \
            --reference 220419_reference_set_for_PARTNERS_mafft.fasta \
            --mafft " --globalpair --maxiterate 1000 " \
            --rm_vloops FALSE \
            --controller $controller \
            --walltime_idx $RES \
            --tsi_analysis FALSE \
            --date $DATE
        fi
        ;;
        
        btr)
        echo "----- build trees ----"
        Rscript $software_path/03_make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --controller $controller \
        --walltime_idx $RES
        ;;

        atr)
        conda activate phylostan
        Rscript $software_path/04_analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --normRefFileName "$phyloscanner_path/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv" \
        --outgroupName "A1.UGANDA.2007.p191845.JX236671"  \
        --ratioBlacklistThreshold 0.01 \
        --distanceThreshold '0.02 0.05' \
        --minReadsPerHost 30 \
        --maxReadsPerHost 100 \
        --multinomial \
        --noProgressBars \
        --postHoccountBlacklisting \
        --relaxedAncestry \
        --zerolengthAdjustment \
        --date $DATE \
        --controller $controller \
        --env_name "phyloscanner" \
        --verbose TRUE
        ;;

        tsi)
        echo "----- Run HIV-TSI -----"
        Rscript $software_path/05_TSI_run_predictions.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --relationship_dir $out_dir_rel \
        --TSI_dir $hivtsipath \
        --date $DATE \
        --env_name 'hivphylotsi' \
        --controller $controller \
        --input_samples $input_samples
        ;;

        dti)
        echo "----- get dates of infection -----"
        Rscript $software_path/06_TSI_estimate_dates.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --date $DATE \
        --input_samples $input_samples \
        --controller $controller
        ;;

        *)
        echo "no R script run. STEP does not match any task.\n" 
        ;;
esac

