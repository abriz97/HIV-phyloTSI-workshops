coord_to_gene <- function(coord) {
    out <- rep(NA_character_, length(coord))
    out[coord %between% c(790, 2289)] <- "gag"
    out[coord %between% c(2253, 5093)] <- "pol"
    out[coord %between% c(6225, 7757)] <- "gp120"
    out[coord %between% c(7758, 8792)] <- "gp41"
    out
}

gene_to_midpoint <- function(gene){
    out <- rep(NA_integer_, length(gene))
    out[gene == "genome"] <- 0
    out[gene == "gag"] <- 1540
    out[gene == "pol"] <- 3673
    out[gene == "gp120"] <- 6991
    out[gene == "gp41"] <- 8275
    return(out)
}

get_codon_positions <- function(path = path_hxb2) {
    hxb2 <- data.table::fread(path_hxb2, header = TRUE)
    codon_positions <- hxb2[, list(
        pos = `HXB2 base position`,
        RF1 = `RF1 aa position`,
        RF2 = `RF2 aa position`,
        RF3 = `RF3 aa position`,
        PROTEIN_RF1 = `RF1 protein`,
        PROTEIN_RF2 = `RF2 protein`,
        PROTEIN_RF3 = `RF3 protein`
    ), ] |> subset(
        PROTEIN_RF1 == "gag" | PROTEIN_RF3 %in% c("pol", "gp120", "gp41")
    )
    setkey(codon_positions, pos)
    codon_positions[!is.na(RF1), codonRF1 := 1:.N, by = c("PROTEIN_RF1", "RF1")]
    codon_positions[!is.na(RF2), codonRF2 := 1:.N, by = c("PROTEIN_RF2", "RF2")]
    codon_positions[!is.na(RF3), codonRF3 := 1:.N, by = c("PROTEIN_RF3", "RF3")]
    codon_positions[!is.na(codonRF1) & !is.na(codonRF2)]
    # codon_positions[, table(codonRF1)]
    codon_positions[!is.na(codonRF1) & !is.na(codonRF3)]

    codon_positions[, codon_type := fifelse(codonRF1 %in% c(1, 2) | codonRF2 %in% c(1, 2) | codonRF3 %in% c(1, 2), yes = "12", no = NA_character_)]
    codon_positions[is.na(codon_type), codon_type := fifelse(codonRF1 == 3 | codonRF2 == 3 | codonRF3 == 3, yes = "3", no = NA_character_)]
    return(codon_positions)
}

preprocess_dates <- function(path=path_dates){
    ddates <- fread(path)
    ddates[, host.id := gsub("(-.*$)|(_.*$)", "-fq1", sequence_id)]
    ddates[, midpoint := as.IDate(first_pos_dt - (first_pos_dt - last_neg_dt)/2)]
    ddates[, `:=` (
        tsi_to_midpoint = (visit_dt - midpoint ) / 365.25,
        true_tsi_max = (visit_dt - last_neg_dt ) / 365.25,
        true_tsi_min = (visit_dt - first_pos_dt) / 365.25
    )]
    return(ddates)
}


preprocess_maf <- function(DT = maf_raw, DCODON = codon_positions) {
    maf <- melt(DT, id.vars = "host.id", variable.name = "pos", value.name = "maf")
    maf[, pos := as.integer(pos)]
    maf <- merge(maf, DCODON, by = "pos")
    return(maf)
}

make_palette_based_on_dates <- function(DT=ddates){
    with(DT[!is.na(tsi_to_midpoint)], {
        out <- rgb(tsi_to_midpoint, 0,tsi_to_midpoint, alpha=.5, maxColorValue=max(tsi_to_midpoint))
        names(out) <- host.id
        out
    })-> palette
    return(palette)
}

maf_averages_by_gene <- function(DT = maf) {
    rbind(
        DT[PROTEIN_RF1 == "gag", list(gene = "gag", midpoint = 1540, maf = mean(maf, na.rm = TRUE)), by = c("host.id", "codon_type")],
        DT[PROTEIN_RF3 == "pol", list(gene = "pol", midpoint = 3673, maf = mean(maf, na.rm = TRUE)), by = c("host.id", "codon_type")],
        DT[PROTEIN_RF3 == "gp120", list(gene = "gp120", midpoint = 6991, maf = mean(maf, na.rm = TRUE)), by = c("host.id", "codon_type")],
        DT[PROTEIN_RF3 == "gp41", list(gene = "gp41", midpoint = 8275, maf = mean(maf, na.rm = TRUE)), by = c("host.id", "codon_type")]
    )
}

plot_average_maf_gene <- function(DT = maf_genes) {
    ggplot(DT, aes(x = midpoint, y = maf, color = host.id)) +
        geom_line(alpha = .2) +
        geom_point() +
        facet_grid(codon_type ~ .) +
        scale_x_continuous(breaks = c(DT[, unique(midpoint)]), labels = c("gag", "pol", "gp120", "gp41")) + {
            if(exists('palette_hostid')){
                scale_color_manual(values=palette_hostid, na.value = 'grey50') 
            }else{ NULL}
        } +
        labs(x = "Gene", y = "Average Minor Allele Frequency in gene") +
        theme_bw() +
        theme(legend.position = "none")
}

plot_maf <- function(DT, window = 10, codon = "none") {
    match.arg(codon, c("12", "3", "none"))

    # subset to codon type of interest
    if (codon == "none") {
        dplot <- copy(DT)
    } else {
        dplot <- subset(DT, codon_type == codon)
    }

    # takee rolling mean
    dplot[, maf_smooth := frollmean(maf, window, align = "center"), by = .(host.id)]

    ggplot(dplot, aes(x = pos, y = maf_smooth, color = host.id)) +
        geom_line(size = .1) +
        theme_bw() +
        labs(
            x = "Genomic position with respect to HXB2",
            y = "Minor allele frequency(Rolling average)"
        ) + {
            if(exists('palette_hostid')){
                scale_color_manual(values=palette_hostid, na.value = 'grey50') 
            }else{ NULL}
        } +
        theme(legend.position = "none")
}

preprocess_patstats <- function(DT=patstats){
    predictors_from_patstats <- c('normalised.largest.rtt', 'tips', 'solo.dual.count')
    patstats_sub <- subset(DT, select=c('host.id', 'xcoord', predictors_from_patstats))
    patstats_sub[, gene := coord_to_gene(xcoord)]
    patstats_sub <- patstats_sub[
        j=lapply(.SD, mean, na.rm=TRUE),
        by=c('host.id', 'gene'), 
        .SDcols=predictors_from_patstats]
    patstats_sub[, midpoint := gene_to_midpoint(gene)]
    patstats_sub <- melt(patstats_sub, id.vars=c('host.id', 'gene','midpoint'))
}

plot_patstats_predictors <- function(DT=patstats_sub){
    ggplot(DT, aes(x=midpoint, y=value, color=host.id)) +
        geom_point() + 
        geom_line(alpha=.2) + 
        facet_grid(variable ~ ., scales="free_y") + {
            if(exists('palette_hostid')){
                scale_color_manual(values=palette_hostid, na.value = 'grey50') 
            }else{ NULL}
        } +
        theme_bw() +
        theme(legend.position="none")
}

plot_predictors_from_tsi_output <- function(DT=dtsi, exclude=character(0)){

    selected_predictors <- c(
        "gag_lrtt",
        "gag_tips",
        "gag_maf3c",
        "pol_lrtt",
        "gp120_lrtt",
        "gp120_tips",
        "gp41_tips",
        "gp41_maf12c",
        "gp41_maf3c"
    )
    tmp <- dtsi |>
        subset(select= !(names(dtsi) %like% 'RF_')) |>
        melt(id.vars = c("host.id"), variable.name = "predictor", value.name = "value") 
    tmp[, `:=` (
        selected=fifelse(predictor %in% selected_predictors, yes=TRUE, no=FALSE),
        gene=gsub("^(.*)_(.*)", "\\1", predictor),
        variable=gsub("^(.*)_(.*)", "\\2", predictor)
    )]
    tmp[, midpoint := gene_to_midpoint(gene)]
    tmp <- subset(tmp, ! variable %in% exclude)
    ggplot(tmp, aes(x=midpoint, y=value, color=host.id)) +
        geom_tile(data=tmp[selected==TRUE], aes(x=midpoint),
        fill='gold', color=NA, width=500,height=Inf, alpha=.01) +
        geom_point() +
        geom_line(data=tmp[gene != 'genome'], alpha=.1) +
        facet_wrap(.~variable, scales='free_y') +
        scale_x_continuous(breaks=tmp[, unique(midpoint)], labels=tmp[, unique(gene)]) +
        scale_y_continuous(expand = expansion(mult=c(0,0.2)) ) +
        scale_color_manual(values=palette_hostid, na.value = 'grey50') +
        theme_bw() +
        labs(x="gene", y="") +
        theme(legend.position = "none")
}

plot_cross_interval_tsisero <- function(DT=da, sqroot=FALSE){

    .f <- function(x) if(sqroot) { sqrt(x) } else {identity(x)}
    .labs <- c(x="Estimated TSI", y="Known TSI")
    if(sqroot){.labs <- paste(.labs, "(sqrt years)") }

    cols <- c(
        "host.id",
        "true_tsi_min",
        "true_tsi_max",
        "tsi_to_midpoint",
        "RF_pred_linear",
        "RF_pred_min_linear",
        "RF_pred_max_linear"
    )

    dplot <- subset(DT, ! is.na(RF_pred_linear), select=cols) 
    # dplot <- cbind(dplot, daterange2tsirange(dplot)) 
    cols2 <- setdiff(cols, "host.id")
    dplot[, (cols2) := lapply(.SD, .f), .SDcols=cols2]

    ALPHA=.1
    ggplot(dplot, aes(x=RF_pred_linear, y=tsi_to_midpoint, color=host.id)) +
        geom_abline(slope=1, intercept=0, linetype='dashed') +
        geom_point() +
        geom_errorbar(aes(ymin=true_tsi_min , ymax=true_tsi_max),alpha=ALPHA) +
        geom_errorbarh(aes(xmin=RF_pred_min_linear, xmax=RF_pred_max_linear),alpha=ALPHA) + {
            if(exists('palette_hostid')){
                scale_color_manual(values=palette_hostid, na.value = 'grey50') 
            }else{ NULL}
        } +
        # stat_smooth(method = "lm", se=FALSE) + 
        theme_bw()  +
        theme(legend.position='none') +
        labs(x=.labs[1], y=.labs[2])
}

plot_histogram_tsi <- function(DT){
    p <- ggplot(dall, aes(x=RF_pred_linear, fill=host.id)) +
        geom_histogram() +
        scale_color_manual(values=palette_hostid, na.value = 'grey50') +
        scale_fill_manual(values=palette_hostid, na.value = 'grey50') +
        scale_y_continuous(expand = expansion(mult=c(0,0.2)) ) +
        theme_bw() +
        theme(legend.position='none') +
        labs( x="Estimated TSIs", y="Count", title="Shortest predictions are from seroconverters")
    return(p)
}
