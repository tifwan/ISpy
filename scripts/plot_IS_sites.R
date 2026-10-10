# Plot IS insertion sites with their neighbouring genes.
# Refactored from plot_figure.R; run by rules/plot.smk through Snakemake `script:`.
#   neighbourhood: genes within `flank` bp of each insertion, one panel per site
#   overview:      insertion sites along each contig
#   compressed:    all windows on one axis, gaps cut out

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gggenes)
  library(ggrepel)
})

# ---- settings -------------------------------------------------------------
sites_file <- snakemake@input[["sites"]]      # make_sites_bed.py
genes_file <- snakemake@input[["genes"]]      # gff_to_bed.py
fai_file   <- snakemake@input[["fai"]]        # contig lengths
fam_file   <- snakemake@input[["families"]]   # ISOSDB id -> family (columns IS, IS_fam)
flank      <- as.numeric(snakemake@params[["flank"]])   # bp shown either side of each insertion
gap        <- 1000                            # width of the cut-out between segments, in plot bp

status_cols <- c("disrupted by IS" = "#d7301f",
                 "flanking gene"   = "#74a9cf",
                 "hypothetical"    = "grey80")

# ---- helpers ----------------------------------------------------------------
# gi|661922017|gb|CP008827.1| -> CP008827.1 (same as the python scripts)
norm_seqid <- function(x) vapply(strsplit(x, "|", fixed = TRUE),
                                 function(p) if (length(p)) p[length(p)] else "", "")

read_bed <- function(f, cols) {
  if (file.size(f) == 0) {
    d <- as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols))
  } else {
    d <- fread(f, sep = "\t", header = FALSE, quote = "", col.names = cols, colClasses = "character")
  }
  d[, `:=`(start0 = as.numeric(start0), end = as.numeric(end))]
  d
}

fam_tab <- fread(fam_file, colClasses = "character")

# IS family from the hit name: ISFinder IS26_IS6_unknown -> IS6,
# IS_novel_0001_IS3 -> IS3, ISOSDB123 -> family table, No identity -> unidentified
is_family <- function(is) {
  is <- sub("_v[0-9]+$", "", is)              # update_IS_db.py suffix for repeated IDs
  out <- rep("unidentified", length(is))
  known <- !is.na(is) & !(is %in% c("", ".", "No identity"))
  novel <- known & grepl("^(IS_novel|novelIS)_[0-9]+_", is)
  isos  <- known & grepl("^ISOSDB[0-9]+$", is)
  other <- known & !novel & !isos
  out[novel] <- sub("^(IS_novel|novelIS)_[0-9]+_", "", is[novel])
  out[isos]  <- fam_tab$IS_fam[match(is[isos], fam_tab$IS)]
  out[other] <- vapply(strsplit(is[other], "_", fixed = TRUE),
                       function(p) if (length(p) >= 2) p[2] else p[1], "")
  missing <- known & is.na(out)               # ISOSDB id not in the family table
  out[missing] <- is[missing]
  out
}

save_empty <- function(msg) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg) + theme_void()
  for (f in c("neighbourhood", "overview", "compressed")) ggsave(snakemake@output[[f]], p, width = 6, height = 2)
}

# ---- read inputs ----------------------------------------------------------
sites <- read_bed(sites_file, c("seqid", "start0", "end", "name", "score", "strand",
                                "sample", "IS", "alignment"))
genes <- read_bed(genes_file, c("seqid", "start0", "end", "locus", "score", "strand",
                                "gene", "product", "type"))
fai <- fread(fai_file, sep = "\t", header = FALSE, select = 1:2, col.names = c("seqid", "len"),
             colClasses = list(character = 1))
fai[, seqid := norm_seqid(seqid)]

if (nrow(sites) == 0) {
  save_empty("No insertion sites")
  quit(save = "no")
}

sites[, pos := end]
sites[, IS_fam := is_family(IS)]
sites[, is_label := sprintf("%s%s\n%s kb", IS_fam,
                            fifelse(alignment %in% "one side", " (1 end)", ""),
                            format(round(pos / 1e3, 1), big.mark = ","))]

genes[, start := start0 + 1]
genes[product == ".", product := NA_character_]
# label: gene name if there is one, otherwise the locus tag number
genes[, label := fifelse(gene != ".", gene, sub("^[^_]+_", "", locus))]

fams <- sort(setdiff(unique(sites$IS_fam), "unidentified"))
fam_cols <- c(setNames(hcl.colors(max(length(fams), 1), "Dark 3")[seq_along(fams)], fams),
              unidentified = "grey50")

# ---- genes around each insertion -----------------------------------------
nb <- rbindlist(lapply(seq_len(nrow(sites)), function(i) {
  s <- sites[i]
  g <- genes[seqid == s$seqid & end >= s$pos - flank & start <= s$pos + flank]
  if (nrow(g) == 0) return(NULL)
  g[, `:=`(name  = s$name,
           pos   = s$pos,
           site  = sprintf("%s  %s  (%s)", s$name, s$IS_fam, format(s$pos, big.mark = ",")),
           x1    = pmax(start - s$pos, -flank),   # clip to the window
           x2    = pmin(end   - s$pos,  flank),
           hit   = start <= s$pos & end >= s$pos)]
  g
}))

if (nrow(nb) == 0) {
  nb <- data.table(seqid = character(), locus = character(), start = numeric(), end = numeric(),
                   strand = character(), label = character(), product = character(),
                   name = character(), pos = numeric(), hit = logical())
} else {
  nb[, status := fcase(hit,                                  "disrupted by IS",
                       product %in% "hypothetical protein",  "hypothetical",
                       default =                             "flanking gene")]
  nb[, site := factor(site, levels = unique(site[order(pos)]))]
}

# ---- plot 1: gene neighbourhood per insertion ----------------------------
if (nrow(nb) == 0) {
  p1 <- ggplot() + annotate("text", x = 0, y = 0,
                            label = sprintf("No genes within %s bp of any insertion", flank)) +
    theme_void()
  ggsave(snakemake@output[["neighbourhood"]], p1, width = 6, height = 2)
} else {
  p1 <- ggplot(nb, aes(xmin = x1, xmax = x2, y = site,
                       fill = status, forward = strand == "+")) +
    geom_gene_arrow(arrowhead_height = unit(4, "mm"),
                    arrowhead_width  = unit(1.5, "mm")) +
    geom_gene_label(aes(label = label), align = "centre",
                    min.size = 3, grow = FALSE) +
    geom_point(data = unique(nb[, .(site)]), aes(x = 0, y = site),
               inherit.aes = FALSE, shape = 25, size = 3.5,
               fill = "black", position = position_nudge(y = 0.35)) +
    facet_wrap(~ site, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = status_cols) +
    scale_x_continuous(limits = c(-flank, flank),
                       labels = function(x) paste0(x / 1000, " kb")) +
    labs(x = "Distance from insertion site", y = NULL, fill = NULL) +
    theme_genes() +
    theme(strip.text   = element_text(hjust = 0, face = "bold"),
          axis.text.y  = element_blank(),
          legend.position = "top")
  ggsave(snakemake@output[["neighbourhood"]], p1,
         width = 9, height = 1.3 * nlevels(nb$site) + 1, limitsize = FALSE)
}

# ---- plot 2: genome-wide overview ----------------------------------------
hits <- nb[hit == TRUE, .(genes = paste(unique(label), collapse = ", ")), by = .(name, pos)]
ov <- merge(sites[, .(seqid, name, pos, IS_fam)], hits, by = c("name", "pos"), all.x = TRUE)
ov[is.na(genes), genes := "intergenic"]
ov[, lab := paste0(IS_fam, "\n", genes)]
glen <- fai[seqid %in% ov$seqid]

p2 <- ggplot(ov, aes(x = pos / 1e6)) +
  geom_segment(data = glen, aes(x = 0, xend = len / 1e6, y = 0, yend = 0),
               linewidth = 4, colour = "grey75", lineend = "round", inherit.aes = FALSE) +
  geom_segment(aes(xend = pos / 1e6, y = 0, yend = 1), colour = "grey30") +
  geom_point(aes(y = 1, colour = IS_fam), size = 3) +
  geom_text_repel(aes(y = 1, label = lab), size = 3, nudge_y = 0.5,
                  direction = "x", segment.colour = "grey60") +
  facet_wrap(~ seqid, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = fam_cols) +
  scale_y_continuous(limits = c(-0.2, 2.2)) +
  scale_x_continuous(name = "Position (Mb)") +
  expand_limits(x = 0) +
  labs(colour = "IS family") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top")
ggsave(snakemake@output[["overview"]], p2, width = 10, height = 1 + 2.5 * nrow(glen), limitsize = FALSE)

# ---- plot 3: compressed view ---------------------------------------------
# windows around each site, merged where they overlap (per contig)
win <- sites[order(seqid, pos), .(seqid, name, pos, IS_fam, is_label,
                                  ws = pmax(pos - flank, 1), we = pos + flank)]
win[, seg_local := cumsum(c(TRUE, ws[-1] > cummax(we)[-.N])), by = seqid]
win[, seg := .GRP, by = .(seqid, seg_local)]

segs <- win[, .(seqid = seqid[1], s_start = min(ws), s_end = max(we)), by = seg]
setorder(segs, seg)
segs[, width  := s_end - s_start]
segs[, offset := c(0, cumsum(head(width, -1) + gap))]
segs[, `:=`(mid = offset + width / 2,
            coord = sprintf("%s %s-%s kb", seqid, format(round(s_start / 1e3, 1), big.mark = ","),
                            format(round(s_end / 1e3, 1), big.mark = ",")))]

win <- merge(win, segs[, .(seg, s_start, offset)], by = "seg")
win[, x := offset + pos - s_start]            # site position on the compressed axis

# genes, one row each, placed on the compressed axis
genes_nb <- nb[, .(hit = any(hit)), by = .(seqid, locus, start, end, strand, label, product)]
gg <- rbindlist(lapply(seq_len(nrow(segs)), function(k) {
  s <- segs[k]
  g <- genes_nb[seqid == s$seqid & end >= s$s_start & start <= s$s_end]
  if (nrow(g) == 0) return(NULL)
  g[, `:=`(seg = s$seg,
           x1  = s$offset + pmax(start, s$s_start) - s$s_start,
           x2  = s$offset + pmin(end,   s$s_end)   - s$s_start)]
  g
}))
breaks <- segs[-1, .(bx = offset - gap / 2)]

p3 <- ggplot() +
  geom_segment(data = segs, aes(x = offset, xend = offset + width, y = 1, yend = 1),
               colour = "grey60", linewidth = 0.6) +
  geom_text(data = segs, aes(x = mid, y = 0.9, label = coord), size = 2.5, colour = "grey40")
if (nrow(gg) > 0) {
  gg[, status := fcase(hit,                                 "disrupted by IS",
                       product %in% "hypothetical protein", "hypothetical",
                       default =                            "flanking gene")]
  p3 <- p3 +
    geom_gene_arrow(data = gg, aes(xmin = x1, xmax = x2, y = 1,
                                   fill = status, forward = strand == "+"),
                    arrow_body_height = unit(9, "mm"),
                    arrowhead_height  = unit(12, "mm"),
                    arrowhead_width   = unit(4, "mm")) +
    geom_gene_label(data = gg, aes(xmin = x1, xmax = x2, y = 1, label = label),
                    height = unit(9, "mm"), align = "centre",
                    size = 12, min.size = 6, grow = FALSE,
                    padding.x = unit(1, "mm")) +
    scale_fill_manual(values = status_cols)
}
if (nrow(breaks) > 0) {
  p3 <- p3 + geom_text(data = breaks, aes(x = bx, y = 1, label = "//"), size = 5, fontface = "bold")
}
p3 <- p3 +
  geom_point(data = win, aes(x = x, y = 1.02), shape = 25, size = 3, fill = "black") +
  geom_text_repel(data = win, aes(x = x, y = 1.03, label = is_label, colour = IS_fam),
                  size = 3, angle = 45, hjust = 0, lineheight = 0.9, fontface = "bold",
                  direction = "x", nudge_y = 0.02, ylim = c(1.02, NA),
                  segment.colour = "grey60", segment.size = 0.3,
                  min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = fam_cols) +
  scale_y_continuous(limits = c(0.85, 1.4)) +
  coord_cartesian(clip = "off") +
  labs(fill = NULL) +
  theme_void() +
  theme(legend.position = "top", plot.margin = margin(30, 20, 10, 20))

ggsave(snakemake@output[["compressed"]], p3,
       width = max(10, sum(segs$width + gap) / 2500), height = 4.5, limitsize = FALSE)
