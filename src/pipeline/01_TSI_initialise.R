self_relative_path <- 'pipeline/01_TSI_initialise.R'
cat( sprintf("=== Start of %s===\n", self_relative_path) )

library(data.table)
library(ggplot2)
library(optparse)
library(magrittr)

option_list <- list(
  make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'pkg.dir'
  ),
  make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  make_option(
    "--cluster_size",
    type = "integer",
    default = 50L,
    help = "Minimum cluster size for host ids groupings[default %default]",
    dest = "cluster_size"
  ),
  make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )
)

args <- parse_args(OptionParser(option_list = option_list))
print(args)

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

make.clusters <- function(DT) {
    idx <- unique(DT$UNIT_ID)
    NCLU <- DT[, floor(length(idx) / args$cluster_size) ]
    if( NCLU > 0){
        dclus <- DT[, list(ID=sample(idx), IDCLU=1:NCLU)]
    }else{
        dclus <- DT[, list(ID=sample(idx), IDCLU=1)]
    }
    dclus[, CLU_SIZE:=.N, by='IDCLU']
    setkey(dclus, IDCLU)

    # if(args$include.twosample.individuals.only)
    # {
    #     cat('Samples from same individual run in different clusters\n',
    #         '\tCheck whether this done on purpose.\n')
    #     dclus1 <- dclus[, list(
    #         ID = paste0(ID, '+'),
    #         IDCLU = IDCLU + max(IDCLU),
    #         CLU_SIZE = CLU_SIZE
    #     )]
    #     dclus <- rbind(dclus, dclus1)
    # }
    dclus
}



usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live' 
    tanya.rakai.dir <- '~/git/HIV-phyloTSI-main/RakExample_Tanya'
}else{
    indir.deepsequence_analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
    indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
    tanya.rakai.dir <- '~/git/HIV-phyloTSI/RakExample_Tanya'
}

# make input samples

mapped_samples <- file.path( indir.deepsequencedata, 'PANGEA2_KenyaWorkshop') %>%
    list.files( pattern='mapped_samples', full.names=TRUE) %>%
    readRDS()

phsc_samples <- mapped_samples[, j={
    list(
        PANGEA_ID = pangea_id,
        UNIT_ID = pangea_id,
        RENAME_ID = paste0(pangea_id, '-fq1'),
        SAMPLE_ID = file.path('PANGEA2_KenyaWorkshop/shiver_output',sequence_id),
        BF = file.path("PANGEA2_KenyaWorkshop/shiver_output", paste0(sequence_id,"_BaseFreq_WithHXB2.csv" ))
    )
}]

saveRDS(phsc_samples, file.path(args$out.dir, 'phsc_input_samples_with_bf.rds'))

# Make clusters.rds
# ______________________________

dclus <- make.clusters(phsc_samples[! UNIT_ID %like% '\\+', ])
filename=file.path(args$out.dir, 'potential_network', 'clusters.rds')
dir.create(dirname(filename))
sprintf("Saving %s\n", filename) %>% cat()
saveRDS(dclus, filename)


# make phscinput_runs_clusize...
# ______________________________
suffix <- phsc_samples[, .(UNIT_ID,PANGEA_ID, RENAME_ID, SAMPLE_ID)]
dclus[, `:=` (PTY_RUN=IDCLU,  PTY_SIZE=CLU_SIZE) ]
dclus <- merge(dclus, suffix, by.x='ID', by.y='UNIT_ID')
setnames(dclus, 'ID', 'UNIT_ID')
setkey(dclus, IDCLU)
filename=file.path(args$out.dir,
                   paste0('phscinput_runs_clusize_', args$cluster_size,'_ncontrol_0.rds'))
sprintf("Saving %s\n", filename) %>% cat()
saveRDS(dclus, filename)

qsub.next.step(file=args$controller,
               next_step='ali', 
               res=1, 
               redo=0
)
