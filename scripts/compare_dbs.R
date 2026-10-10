# Compare ISFinder and ISOSDB calls for one sample (from plot_figure.R):
# per-site table with both families and whether they agree, plus a family Venn.
# Run by rules/plot.smk through Snakemake `script:`.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggvenn)
})

fam_tab <- fread(snakemake@input[["families"]], colClasses = "character")

# same as plot_IS_sites.R
is_family <- function(is) {
  is <- sub("_v[0-9]+$", "", is)
  out <- rep("unidentified", length(is))
  known <- !is.na(is) & !(is %in% c("", ".", "No identity"))
  novel <- known & grepl("^(IS_novel|novelIS)_[0-9]+_", is)
  isos  <- known & grepl("^ISOSDB[0-9]+$", is)
  other <- known & !novel & !isos
  out[novel] <- sub("^(IS_novel|novelIS)_[0-9]+_", "", is[novel])
  out[isos]  <- fam_tab$IS_fam[match(is[isos], fam_tab$IS)]
  out[other] <- vapply(strsplit(is[other], "_", fixed = TRUE),
                       function(p) if (length(p) >= 2) p[2] else p[1], "")
  missing <- known & is.na(out)
  out[missing] <- is[missing]
  out
}

keys <- c("Sample", "Chromosome", "Start_Position", "Stop_Position", "Potential_sequence")
read_hits <- function(f, db) {
  d <- fread(f, sep = "\t", quote = "", colClasses = "character")
  d <- d[, c(keys, "Potential_IS", "Alignment"), with = FALSE]
  setnames(d, c("Potential_IS", "Alignment"), c(db, paste0(db, "_alignment")))
  d
}

all_db <- merge(read_hits(snakemake@input[["isfinder"]], "ISFinder"),
                read_hits(snakemake@input[["isosdb"]], "ISOSDB"),
                by = keys, all = TRUE)
all_db[, `:=`(ISFinder_fam = is_family(ISFinder),
              ISOSDB_fam   = is_family(ISOSDB),
              site         = pmax(as.numeric(Start_Position), as.numeric(Stop_Position)))]
all_db[, agree := ISFinder_fam == ISOSDB_fam & ISFinder_fam != "unidentified"]
setorder(all_db, Chromosome, site)
fwrite(all_db, snakemake@output[["table"]], sep = "\t")

fam_sets <- list(ISFinder = setdiff(unique(all_db$ISFinder_fam), "unidentified"),
                 ISOSDB   = setdiff(unique(all_db$ISOSDB_fam), "unidentified"))

if (all(lengths(fam_sets) == 0)) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = "No IS families identified") + theme_void()
} else {
  p <- ggvenn(fam_sets, show_elements = TRUE, label_sep = "\n",
              fill_color = c("#74a9cf", "#fdae6b"), text_size = 5)
}
ggsave(snakemake@output[["venn"]], p, width = 7, height = 6)
