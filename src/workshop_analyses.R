# instal required packages
# ________________________

# install.packages("here")
# install.packages("data.table")
# install.packages("ggplot2")

# Load packages
library(here)
library(data.table)
library(ggplot2) 

# Set working directory to root of repository
git_root <- here()
stopifnot(basename(git_root) == "HIV-phyloTSI-workshops")
git_data <- file.path(git_root, "data")
source(file.path(git_root, "src/R/workshop_R_helpers.R"))
git_figures <-file.path(git_root, "slides/figures")

# paths to input filesp
path_patstats   <- file.path(git_data, "ptyr1_patStats.csv")
path_maf        <- file.path(git_data, "phsc_input_samples_maf.csv")
path_hxb2       <- file.path(git_data, "HXB2_refdata.csv")
path_dates      <- file.path(git_data,  "2023-08-17_pangea2_popart_training_dates.csv")
path_tsi        <- file.path(git_data, "ptyr1_tsi.csv")

# load and preprocess inputs
patstats    <- data.table::fread(path_patstats, header=TRUE)
maf_raw     <- data.table::fread(path_maf, header=TRUE)
setnames(maf_raw, "pos", "host.id")
ddates <- preprocess_dates(path_dates)
stopifnot(all(ddates$host.id %in% maf_raw$host.id))
palette_hostid <- make_palette_based_on_dates(ddates)

ddates[!is.na(first_pos_dt), table(visit_dt == first_pos_dt)]

# plotting predictors
codon_positions <- get_codon_positions(path_hxb2)
maf             <- preprocess_maf(maf_raw)
maf_genes       <- maf_averages_by_gene(maf)
p_maf           <- plot_maf(maf, codon="12", window=50)
p_maf_avg       <- plot_average_maf_gene(maf_genes)
patstats_sub    <- preprocess_patstats(patstats) 
p_patstats      <- plot_patstats_predictors(patstats_sub)

# plotting tsi estimates
dtsi <- fread(path_tsi)
p_preds <- plot_predictors_from_tsi_output(dtsi, exclude="dual")

# merge tsi estimates with dates 
dall <- merge(dtsi, ddates)

# final plots
p_cross <- plot_cross_interval_tsisero(dall, sqroot=TRUE)
p_hist <- ggplot(dall, aes(x=RF_pred_linear, fill=host.id)) +
    geom_histogram() +
    scale_color_manual(values=palette_hostid, na.value = 'grey50') +
    scale_fill_manual(values=palette_hostid, na.value = 'grey50') +
    scale_y_continuous(expand = expansion(mult=c(0,0.2)) ) +
    theme_bw() +
    theme(legend.position='none') +
    labs( x="Estimated TSIs", y="Count", title="Shortest predictions are from seroconverters")

# save plots
ggsave(p_maf,   filename=file.path(git_figures, "maf.png"), w=10, h=8)
ggsave(p_hist,  filename=file.path(git_figures, "histogram_phyloTSI.png"), w=10, h=8)
ggsave(p_cross, filename=file.path(git_figures, "crosscheck.png"), w=10, h=8)
ggsave(p_preds, filename=file.path(git_figures, "predictors.png"), w=10, h=8)
