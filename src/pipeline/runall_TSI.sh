#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI_seroconv2.sh
# If unset default to "sim"
if [ -z "$STEP" ]
then
        echo "Intended use:\n"
        echo 'qsub -v STEP="xxx" runall_TSI_seroconv2.sh'
        exit 1
fi
echo "running '${STEP:=sim}' analysis"

# This includes all code necessary to run PHSC pipeline to produce TSI estimates
DEEPDATA="/rds/general/project/ratmann_pangea_deepsequencedata/live"
DEEPANALYSES="/rds/general/project/ratmann_deepseq_analyses/live"
HOME="/rds/general/user/ab1820/home"

software_path="$HOME/git/HIV-phyloTSI-workshops/src/pipeline"
phyloscanner_path="$HOME/git/phyloscanner"
hivtsipath="$HOME/git/HIV-phyloTSI"
out_dir_base="$DEEPANALYSES/PANGEA2_KenyaWorkshop"
out_dir_rel="$out_dir_base/TODO"
controller="$software_path/runall_TSI.sh"
CLUSIZE='50'
DATE='2022-04-25'


cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments


case $STEP in

        # In this analysis we avoid the first step of computing similarities, as there exist already
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
        Rscript $software_path/02_make_deep_sequence_alignments.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --sliding_width 10 \
        --n_control 0 \
        --cluster_size $CLUSIZE \
        --tsi_analysis TRUE \
        --controller $controller
        ;;
        
        btr)
        echo "----- build trees ----"
        Rscript $software_path/03_make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --controller $controller
        ;;

        ctr)
        echo "----- check trees ----"
        Rscript $software_path/04_check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --controller $controller
        ;;

        atr)
        conda activate phylostan
        Rscript $software_path/05_analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --blacklistReport \
        --outputNexusTree \
        --skipSummaryGraph \
        --normRefFileName "$phyloscanner_path/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv" \
        --outgroupName "A1.UGANDA.2007.p191845.JX236671"  \
        --ratioBlacklistThreshold 0.005 \
        --date $DATE \
        --env_name "phylostan" \
        --controller $controller
        ;;

        tsi)
        echo "----- Run HIV-TSI -----"
        Rscript $software_path/06_TSI_run_predictions.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --input_samples "$out_dir_base/220419_phscinput_samples.rds" \
        --TSI_dir $hivtsipath \
        --date $DATE \
        --env_name 'hivphylotsi' \
        --controller $controller
        ;;

        dti)
        echo "----- get dates of infection -----"
        Rscript $software_path/07_TSI_estimate_dates.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --date $DATE \
        --input_samples "$out_dir_base/220419_phscinput_samples.rds" \
        --controller $controller
        ;;

        pst)
        echo "----- post processing pngs -----"
        conda activate phylostan # required for gridExtra package
        Rscript $software_path/08_TSI_postprocessing_comparison.R \
        --relationship_dir $out_dir_rel \
        --TSI_dir $hivtsipath \
        --input_samples "$out_dir_base/220419_phscinput_samples.rds" \
        --controller $controller \
        ;;

        *)
        echo "no R script run. STEP does not match any task.\n" 
        ;;
esac

