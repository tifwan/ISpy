IS_gene_annotation <- fread("/Users/tifwan/Downloads/clostet_gene.tsv")



# Plot IS insertion sites with their neighbouring genes
#install.packages(c("data.table", "ggplot2", "gggenes", "ggrepel"))

library(data.table)
library(ggplot2)
library(gggenes)
library(ggrepel)
library(ggvenn)

# ---- settings -------------------------------------------------------------
gff_file   <- "/Users/tifwan/Downloads/genes.gff"   # the cleaned GFF (no FASTA block)
sites_file <- "/Users/tifwan/Downloads/PCMP_H482_sites.bed"          # BED made from the panISa table
genome_len <- 5394056                 # CP008827.1 length (##sequence-region)
flank      <- 3000                    # bp shown either side of each insertion
ISFinder <- "/Users/tifwan/University of Michigan Dropbox/Tiffany Wan/PhD/Project_ASM_hackathon2026/output/panisa/PCMP_H482_ISfinder.tsv"
ISOSDB <- "/Users/tifwan/University of Michigan Dropbox/Tiffany Wan/PhD/Project_ASM_hackathon2026/output/panisa/PCMP_H482_ISOSDB.tsv"
isosdb_annot <- "/Users/tifwan/University of Michigan Dropbox/Tiffany Wan/PhD/Project_ASM_hackathon2026/output/IS_fam_annot.txt"
# ---- read inputs ----------------------------------------------------------

ISFinder_df <- fread(ISFinder)
ISFinder_df <- ISFinder_df %>%
  rename(ISFinder = Potential_IS)
ISOSDB_df <- fread(ISOSDB)
isosdb_annot <- fread(isosdb_annot)
ISOSDB_df <- ISOSDB_df %>%
  merge(isosdb_annot, by.x = "Potential_IS",by.y = "IS",all = TRUE)
ISOSDB_df <- ISOSDB_df %>%
  rename(ISOSDB=IS_fam) %>%
  select(-Potential_IS)

all_db <- all_db %>%
  mutate(ISFinder_fam = str_split_i(Potential_IS, "_", 2),
         site = format(pmax(Start_Position, Stop_Position), big.mark = ","),
         agree = ISFinder_fam == ISOSDB)

all_db %>% select(site, ISFinder_fam, ISOSDB, agree)

fam_sets <- list(ISFinder = unique(na.omit(all_db$ISFinder_fam)),
                 ISOSDB   = unique(na.omit(all_db$ISOSDB)))

ggvenn(fam_sets, show_elements = TRUE, label_sep = "\n",
       fill_color = c("#74a9cf", "#fdae6b","pink"), text_size =8)
all_db <- merge(ISFinder_df, ISOSDB_df, by = c("Sample","Chromosome","Start_Position","Stop_Position","Potential_sequence"))
gff <- fread(gff_file, sep = "\t", header = FALSE, quote = "",
             col.names = c("seqid", "source", "type", "start", "end",
                           "score", "strand", "phase", "attr"))


sites <- fread(sites_file, sep = "\t", header = FALSE,
               col.names = c("seqid", "start0", "end", "name"))
sites[, pos := end]

# pull a key=value field out of the GFF attribute column
get_attr <- function(x, key) {
  has <- grepl(paste0("(^|;)", key, "="), x)
  out <- rep(NA_character_, length(x))
  m <- regmatches(x, regexpr(paste0("(^|;)", key, "=[^;]+"), x))
  out[has] <- sub(paste0("^;?", key, "="), "", m)
  out
}

gff[, `:=`(gene    = get_attr(attr, "gene"),
           locus   = get_attr(attr, "locus_tag"),
           product = get_attr(attr, "product"))]
gff[!is.na(product), product := vapply(product, URLdecode, "")]  # Prokka encodes commas etc.
# label: gene name if there is one, otherwise the locus tag number
gff[, label := fifelse(!is.na(gene), gene, sub("^[^_]+_", "", locus))]

sites[, seqid := vapply(strsplit(seqid, "|", fixed = TRUE), `[`, "", 4)]
unique(sites$seqid)   # should now be "CP008827.1"
# ---- genes around each insertion -----------------------------------------
nb <- rbindlist(lapply(seq_len(nrow(sites)), function(i) {
  s <- sites[i]
  g <- gff[seqid == s$seqid & end >= s$pos - flank & start <= s$pos + flank]
  if (nrow(g) == 0) return(NULL)
  g[, `:=`(name  = s$name,
           pos   = s$pos,
           site  = sprintf("%s  (%s)", s$name, format(s$pos, big.mark = ",")),
           x1    = pmax(start - s$pos, -flank),   # clip to the window
           x2    = pmin(end   - s$pos,  flank),
           hit   = start <= s$pos & end >= s$pos)]
  g
}))

nb[, status := fcase(hit,                                  "disrupted by IS",
                     product %in% "hypothetical protein",  "hypothetical",
                     default =                             "flanking gene")]
nb[, site := factor(site, levels = unique(site[order(pos)]))]

# ---- plot 1: gene neighbourhood per insertion ----------------------------
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
  scale_fill_manual(values = c("disrupted by IS" = "#d7301f",
                               "flanking gene"   = "#74a9cf",
                               "hypothetical"    = "grey80")) +
  scale_x_continuous(limits = c(-flank, flank),
                     labels = function(x) paste0(x / 1000, " kb")) +
  labs(x = "Distance from insertion site", y = NULL, fill = NULL) +
  theme_genes() +
  theme(strip.text   = element_text(hjust = 0, face = "bold"),
        axis.text.y  = element_blank(),
        legend.position = "top")

ggsave("is_sites_neighbourhood.pdf", p1,
       width = 9, height = 1.3 * nlevels(nb$site) + 1)

# ---- plot 2: genome-wide overview ----------------------------------------
hits <- nb[hit == TRUE, .(genes = paste(unique(label), collapse = ", ")),
           by = .(name, pos)]
ov <- merge(sites[, .(name, pos)], hits, by = c("name", "pos"), all.x = TRUE)
ov[is.na(genes), genes := "intergenic"]
ov[, lab := paste0(name, "\n", genes)]

p2 <- ggplot(ov, aes(x = pos / 1e6)) +
  annotate("segment", x = 0, xend = genome_len / 1e6, y = 0, yend = 0,
           linewidth = 4, colour = "grey75", lineend = "round") +
  geom_segment(aes(xend = pos / 1e6, y = 0, yend = 1), colour = "grey30") +
  geom_point(aes(y = 1), colour = "#d7301f", size = 3) +
  geom_text_repel(aes(y = 1, label = lab), size = 3, nudge_y = 0.5,
                  direction = "x", segment.colour = "grey60") +
  scale_y_continuous(limits = c(-0.2, 2.2)) +
  scale_x_continuous(limits = c(0, genome_len / 1e6), name = "Position (Mb)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
p2
ggsave("is_sites_overview.pdf", p2, width = 10, height = 3.5)

library(data.table)
library(ggplot2)
library(gggenes)

gap <- 1000   # width of the cut-out between segments, in plot bp

# ---- 1. windows around each site, merged where they overlap --------------
win <- sites[order(pos), .(name, pos, ws = pmax(pos - flank, 1), we = pos + flank)]
win[, seg := cumsum(c(TRUE, ws[-1] > cummax(we)[-.N]))]

segs <- win[, .(s_start = min(ws), s_end = max(we)), by = seg]
segs[, width  := s_end - s_start]
segs[, offset := c(0, cumsum(head(width, -1) + gap))]

win <- merge(win, segs, by = "seg")
win[, x := offset + pos - s_start]            # site position on the compressed axis

# ---- 2. genes, one row each, placed on the compressed axis ---------------
genes <- nb[, .(hit = any(hit)), by = .(locus, start, end, strand, label, product)]

gg <- rbindlist(lapply(seq_len(nrow(segs)), function(k) {
  s <- segs[k]
  g <- genes[end >= s$s_start & start <= s$s_end]
  g[, `:=`(seg = s$seg,
           x1  = s$offset + pmax(start, s$s_start) - s$s_start,
           x2  = s$offset + pmin(end,   s$s_end)   - s$s_start)]
  g
}))
gg[, status := fcase(hit,                                 "disrupted by IS",
                     product %in% "hypothetical protein", "hypothetical",
                     default =                            "flanking gene")]

# segment labels (genomic coordinates) and break marks
segs[, `:=`(mid = offset + width / 2,
            coord = sprintf("%s–%s kb", format(round(s_start / 1e3, 1), big.mark = ","),
                            format(round(s_end   / 1e3, 1), big.mark = ",")))]
breaks <- segs[-1, .(bx = offset - gap / 2)]

# ---- 3. plot --------------------------------------------------------------

lab <- as.data.table(all_db)[, .(pos       = pmax(Start_Position, Stop_Position),
                                 IS_fam    = ISOSDB,
                                 isos_side = Alignment.y)]

win <- merge(win, lab, by = "pos", all.x = TRUE)

win[, IS_fam := fifelse(is.na(IS_fam), "unidentified", IS_fam)]
win[, is_label := sprintf("%s%s\n%s kb",
                          IS_fam,
                          fifelse(isos_side %in% "one side", " (1 end)", ""),
                          format(round(pos / 1e3, 1), big.mark = ","))]
p1 <- ggplot() +
  geom_segment(data = segs, aes(x = offset, xend = offset + width, y = 1, yend = 1),
               colour = "grey60", linewidth = 0.6) +
  geom_gene_arrow(data = gg, aes(xmin = x1, xmax = x2, y = 1,
                                 fill = status, forward = strand == "+"),
                  arrow_body_height = unit(9, "mm"),
                  arrowhead_height  = unit(12, "mm"),
                  arrowhead_width   = unit(4, "mm")) +
  geom_gene_label(data = gg, aes(xmin = x1, xmax = x2, y = 1, label = label),
                  height = unit(9, "mm"), align = "centre",
                  size = 12, min.size = 6, grow = FALSE,
                  padding.x = unit(1, "mm")) +
  geom_point(data = win, aes(x = x, y = 1.02), shape = 25, size = 3, fill = "black") +
  geom_text_repel(data = win, aes(x = x, y = 1.03, label = is_label, colour = IS_fam),
                  size = 3, angle = 45, hjust = 0, lineheight = 0.9, fontface = "bold",
                  direction = "x", nudge_y = 0.02, ylim = c(1.02, NA),
                  segment.colour = "grey60", segment.size = 0.3,
                  min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = c("IS5" = "#1b9e77", "IS1" = "#d95f02",
                                 "IS110" = "#7570b3", "IS481" = "#e7298a",
                                 "unidentified" = "grey50")) +
  geom_text(data = breaks, aes(x = bx, y = 1, label = "//"), size = 5, fontface = "bold") +
  scale_fill_manual(values = c("disrupted by IS" = "#d7301f",
                               "flanking gene"   = "#74a9cf",
                               "hypothetical"    = "grey80")) +
  scale_y_continuous(limits = c(0.85, 1.4)) +
  coord_cartesian(clip = "off") +
  labs(fill = NULL) +
  theme_void() +
  theme(legend.position = "top", plot.margin = margin(30, 20, 10, 20))
p1
ggsave("is_sites_compressed.pdf", p1,
       width = max(10, sum(segs$width + gap) / 2500), height = 2.8)


ggsave("Dropbox (University of Michigan)/PhD/Project_ASM_hackathon2026/is_sites_compressed.jpeg", p1,
       width = max(10, sum(segs$width + gap) / 2500), height = 4.5)
