
suppressPackageStartupMessages({
  library(tidyverse); library(ComplexHeatmap); library(circlize)
  library(grid); library(cluster); library(patchwork); library(ggdendro)
})


suppressWarnings(Sys.setlocale("LC_ALL", "en_US.UTF-8"))

# parameters
INPUT       <- "Silvia_project_human_proteome.txt"
P_SIG       <- 0.1            # omnibus + read-out FDR threshold (q-value)
K           <- 7             # number of trend clusters                            
K_SCAN      <- 2:12          # range scanned in the diagnostic
FONT_FAMILY <- "DejaVu Sans"  # renders the Delta glyph under cairo
set.seed(1)

# palette
strain_cols <- c(Control="#000000", AP1="#CA2222", "AP1\u0394SpnA"="#305E9D",
                 BMJ71="#3D8E42", "BMJ71\u0394SpnA"="#7D7E7E")
tp_cols   <- c(`3h`="#B0B0B0", `24h`="#3C3C3C")   
dir_cols  <- c(Up="#B2182B", Down="#2166AC", `n.s.`="#F0F0F0") 
divergent <- c("#2166AC","#F7F7F7","#B2182B")   
accent    <- "#CA2222"                          

# cairo PNG helper 
png_cairo <- function(file, w, h) png(file, width=w, height=h, units="in",
                                      res=300, type="cairo", family=FONT_FAMILY)
theme_paper <- function(base = 11)
  theme_classic(base_size = base) +
  theme(text          = element_text(family = FONT_FAMILY, colour = "black"),
        axis.line     = element_line(colour = "black", linewidth = 0.4),
        axis.ticks    = element_line(colour = "black", linewidth = 0.4),
        axis.title    = element_text(size = 11, face = "bold"),
        axis.text     = element_text(size = 9,  colour = "black"),
        legend.title  = element_text(size = 9,  face = "bold"),
        legend.text   = element_text(size = 9),
        strip.background = element_blank(),
        strip.text    = element_text(size = 9, face = "bold"))


# 1. Parsing the data file

raw <- read_tsv(INPUT, show_col_types = FALSE)
samp_cols <- grep("LJH", names(raw), value = TRUE)                 # 30 samples
expr <- as.matrix(as.data.frame(lapply(raw[, samp_cols], function(x)
          as.numeric(str_replace(as.character(x), ",", ".")))))    # comma decimals
rownames(expr) <- str_split_i(raw[["T: PG.ProteinGroups"]], ";", 1)
colnames(expr) <- samp_cols

core   <- str_match(samp_cols, "LJH_T2504_THP1_AA_(.*?)_S\\d")[, 2]
parts  <- str_match(core, "^(?:(AP1mut|BMJ71mut|AP1|BMJ71)_)?(0h|24h)_(\\d)$")
strain <- parts[, 2]; strain[is.na(strain)] <- "Control"
strain <- recode(strain, "AP1mut" = "AP1\u0394SpnA", "BMJ71mut" = "BMJ71\u0394SpnA")
strain_lv <- names(strain_cols)
meta <- tibble(col = samp_cols,
               strain = factor(strain, levels = strain_lv),
               tp = factor(recode(parts[,3], "0h"="3h"), levels = c("3h","24h")),
               grp = paste(strain, recode(parts[,3], "0h"="3h"), sep = "_"))
stopifnot(!any(is.na(meta$strain)), !any(is.na(meta$tp)))

sym_lookup <- tibble(uniprot = rownames(expr), protein = raw[["T: PG.ProteinNames"]]) |>
  mutate(symbol = str_remove(protein, "_HUMAN.*")) |>
  distinct(uniprot, .keep_all = TRUE)

col_order <- meta |> arrange(strain, tp) |> pull(col)             # strain-major
meta_ord  <- meta[match(col_order, meta$col), ]


# 2. omnibus one-way ANOVA (F-test) across groups to keep proteins that change in ANY group

gf    <- factor(meta$grp[match(colnames(expr), meta$col)])
aov_p <- apply(expr, 1, function(y)
           tryCatch(summary(aov(y ~ gf))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_))
aov_q <- p.adjust(aov_p, "BH")
sel   <- names(aov_q)[which(aov_q < P_SIG)]
message(sprintf("omnibus ANOVA q < %.2g : %d / %d proteins selected", P_SIG, length(sel), nrow(expr)))

zscore <- function(m){ z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
z <- zscore(expr[sel, col_order, drop = FALSE])

# 3. natural ward.D2 clustering of the z-scored rows

d   <- dist(z, method = "euclidean")
hc  <- hclust(d, method = "ward.D2")
cl  <- cutree(hc, k = K)
ord <- names(sort(table(cl), decreasing = TRUE))                  # size-order labels
cl2 <- factor(setNames(paste0("C", seq_along(ord)), ord)[as.character(cl)],
              levels = paste0("C", 1:K))
message(sprintf("clusters (size-ordered): %s",
                paste(levels(cl2), table(cl2), sep=":", collapse="  ")))

# 4. per-protein Welch t-test + BH

tt_dir <- function(a_cols, b_cols){
  A <- expr[sel, a_cols, drop=FALSE]; B <- expr[sel, b_cols, drop=FALSE]
  fc <- rowMeans(A) - rowMeans(B)
  p  <- vapply(seq_along(sel), function(i)
          tryCatch(t.test(A[i,], B[i,])$p.value, error=function(e) NA_real_), numeric(1))
  q  <- p.adjust(p, "BH")
  ifelse(q < P_SIG & fc > 0, "Up", ifelse(q < P_SIG & fc < 0, "Down", "n.s."))
}
cf <- function(s, t) meta$col[meta$strain == s & meta$tp == t]

# infection course: 3h -> 24h per strain (includes Control as the reference)
time_bars <- sapply(strain_lv, function(s) tt_dir(cf(s,"24h"), cf(s,"3h")))
colnames(time_bars) <- strain_lv
# SpnA effect: KO - WT per background per timepoint
spna_bars <- cbind(
  `AP1 3h`   = tt_dir(cf("AP1\u0394SpnA","3h"),  cf("AP1","3h")),
  `AP1 24h`  = tt_dir(cf("AP1\u0394SpnA","24h"), cf("AP1","24h")),
  `BMJ71 3h` = tt_dir(cf("BMJ71\u0394SpnA","3h"),  cf("BMJ71","3h")),
  `BMJ71 24h`= tt_dir(cf("BMJ71\u0394SpnA","24h"), cf("BMJ71","24h")))

# AP1 SpnA-dependence set (KO vs WT, full proteome, BH within timepoint) 
spna_test <- function(tpt){                      # KO - WT at one timepoint, full proteome
  A <- expr[, cf("AP1\u0394SpnA", tpt), drop=FALSE]; B <- expr[, cf("AP1", tpt), drop=FALSE]
  fc <- rowMeans(A) - rowMeans(B)
  p  <- vapply(seq_len(nrow(expr)), function(i)
          tryCatch(t.test(A[i,], B[i,])$p.value, error=function(e) NA_real_), numeric(1))
  tibble(uniprot = rownames(expr), log2FC = fc, p = p, q = p.adjust(p, "BH"), tp = tpt)
}
spna     <- bind_rows(spna_test("3h"), spna_test("24h"))
spna_sig <- spna |> filter(q < P_SIG) |> pull(uniprot) |> unique()
message(sprintf("AP1 \u0394SpnA vs AP1 significant (q<%.2g): 3h=%d  24h=%d  union=%d",
                P_SIG,
                sum(spna$tp=="3h"  & spna$q<P_SIG, na.rm=TRUE),
                sum(spna$tp=="24h" & spna$q<P_SIG, na.rm=TRUE),
                length(spna_sig)))

# 4b. heatmap

zl <- quantile(abs(z), 0.99); col_fun <- colorRamp2(c(-zl,0,zl), divergent)

top_anno <- HeatmapAnnotation(
  Strain = meta_ord$strain, Timepoint = meta_ord$tp,
  col = list(Strain = strain_cols, Timepoint = tp_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"), simple_anno_size = unit(3.5,"mm"),
  annotation_legend_param = list(
    Strain    = list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8)),
    Timepoint = list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8))))

tb <- as.data.frame(time_bars); names(tb) <- paste0(strain_lv, " (t)")
left_anno <- rowAnnotation(df = tb,
  col = setNames(rep(list(dir_cols), ncol(tb)), names(tb)),
  annotation_name_gp = gpar(fontsize=6.5), annotation_name_rot = 45,
  simple_anno_size = unit(2.6,"mm"),
  show_legend = c(TRUE, rep(FALSE, ncol(tb)-1)),
  annotation_legend_param = setNames(list(list(title = "Direction (3h\u219224h / \u0394SpnA)",
    title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8))), names(tb)[1]))
# Right side: instead of a significance track, MARK where the AP1 dSpnA-significant
# proteins fall within the naturally-clustered rows (labelled with UniProt accessions).
mark_idx <- which(rownames(z) %in% spna_sig)
mark_idx <- mark_idx[order(match(rownames(z)[mark_idx], spna_sig))]  # stable order
message(sprintf("Marking %d of %d AP1 \u0394SpnA-significant proteins present in the omnibus-selected rows",
                length(mark_idx), length(spna_sig)))
mark_title <- "AP1 \u0394SpnA sig."          # Delta only in a plain string, not a backtick name
right_anno <- rowAnnotation(
  spna_mark = anno_mark(at = mark_idx, labels = rownames(z)[mark_idx],
    labels_gp = gpar(fontsize = 6.5, fontface = "italic"),
    link_gp = gpar(lwd = 0.5), padding = unit(0.6, "mm")),
  annotation_label = mark_title,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"), annotation_name_rot = 90)

ht <- Heatmap(z, name = "z-score", col = col_fun,
  top_annotation = top_anno, left_annotation = left_anno, right_annotation = right_anno,
  cluster_columns = FALSE, column_split = meta_ord$strain, column_gap = unit(1,"mm"),
  show_column_names = FALSE,
  cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
  row_split = cl2, cluster_row_slices = FALSE,          # keep C1..CK order
  row_title_gp = gpar(fontsize=10, fontface="bold"), row_title_rot = 0,  # labels only, no band
  show_row_names = FALSE, use_raster = TRUE, raster_quality = 3,
  column_title_gp = gpar(fontsize=9, fontface="bold"),
  heatmap_legend_param = list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8)))

draw_ht <- function() draw(ht, heatmap_legend_side="right",
                           annotation_legend_side="right", merge_legend=TRUE)
png_cairo("heatmap_natural_omnibus.png", 10.5, 12); draw_ht(); dev.off()
cairo_pdf("heatmap_natural_omnibus.pdf", width=10.5, height=12, family=FONT_FAMILY); draw_ht(); dev.off()

# 5. mean z-score trajectory per cluster

zt   <- as.data.frame(z); zt$cluster <- cl2; zt$uniprot <- rownames(z)
long <- zt |> pivot_longer(-c(cluster, uniprot), names_to="col", values_to="zval") |>
  left_join(tibble(col = col_order, strain = meta_ord$strain, tp = meta_ord$tp), by="col")
prof <- long |> group_by(cluster, strain, tp) |>
  summarise(mean_z = mean(zval), se = sd(zval)/sqrt(n()), .groups="drop")
sizes <- as.data.frame(table(cl2)) |> setNames(c("cluster","n"))
prof  <- prof |> left_join(sizes, by="cluster") |> mutate(facet = paste0(cluster," (n=",n,")"))

# One panel per cluster, stacked in a SINGLE vertical column (patchwork) so the
# strip sits flush beside the tall omnibus heatmap. Shared y-range; x-axis and
# legend drawn once at the bottom.
yr <- range(c(prof$mean_z - prof$se, prof$mean_z + prof$se))
facet_lv <- levels(prof$cluster)
one_panel <- function(cl_id, is_last){
  d <- filter(prof, cluster == cl_id)
  ttl <- d$facet[1]
  g <- ggplot(d, aes(tp, mean_z, color=strain, group=strain)) +
    geom_hline(yintercept=0, color="gray70", linewidth=0.3, linetype="dotted") +
    geom_ribbon(aes(ymin=mean_z-se, ymax=mean_z+se, fill=strain), alpha=0.15, color=NA) +
    geom_line(linewidth=0.9) + geom_point(size=1.8) +
    scale_color_manual(values=strain_cols, name="Strain") +
    scale_fill_manual(values=strain_cols, guide="none") +
    scale_y_continuous(limits=yr) +
    labs(x=NULL, y=ttl) +                        # cluster label as the y-axis title
    theme_paper() +
    theme(axis.title.y = element_text(size=9, face="bold", angle=90),
          plot.margin  = margin(2,4,2,4))
  if (!is_last) g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                               axis.line.x=element_blank())
  else          g <- g + labs(x="Infection course")
  g
}
panels  <- lapply(seq_along(facet_lv), function(i) one_panel(facet_lv[i], i==length(facet_lv)))
p_prof  <- wrap_plots(panels, ncol=1) +
  plot_layout(guides="collect") &
  theme(legend.position="bottom", legend.text=element_text(size=8),
        legend.title=element_text(size=9,face="bold"))
ph <- max(6, 1.5*length(facet_lv))               # ~1.5 in per cluster, tall to match heatmap
png_cairo("cluster_trend_profiles.png", 3.2, ph); print(p_prof); dev.off()
cairo_pdf("cluster_trend_profiles.pdf", width=3.2, height=ph, family=FONT_FAMILY); print(p_prof); dev.off()

# 6. clustering diagnostic 

zmat <- as.matrix(z)
wss <- sapply(K_SCAN, function(k){
  g <- cutree(hc, k = k)
  sum(sapply(unique(g), function(gi){
    rows <- zmat[g == gi, , drop=FALSE]
    if (nrow(rows) < 2) return(0)
    sum(sweep(rows, 2, colMeans(rows))^2)
  }))
})
sil <- sapply(K_SCAN, function(k) mean(silhouette(cutree(hc, k = k), d)[, "sil_width"]))
diag_df <- tibble(k = K_SCAN, wss = wss, sil = sil)

p_elbow <- ggplot(diag_df, aes(k, wss)) +
  geom_line(color="black", linewidth=0.6) +
  geom_point(size=2, color="black") +
  geom_point(data=filter(diag_df, k==K), size=3.5, shape=21, fill=accent, color="black", stroke=0.6) +
  scale_x_continuous(breaks=K_SCAN) +
  labs(x="Number of clusters (k)", y="Total within-cluster SS") +
  theme_paper()
p_sil <- ggplot(diag_df, aes(k, sil)) +
  geom_line(color="black", linewidth=0.6) +
  geom_point(size=2, color="black") +
  geom_point(data=filter(diag_df, k==K), size=3.5, shape=21, fill=accent, color="black", stroke=0.6) +
  scale_x_continuous(breaks=K_SCAN) +
  labs(x="Number of clusters (k)", y="Average silhouette width") +
  theme_paper()

# dendrogram (ggdendro) with the K-cut height marked
h_cut  <- mean(rev(hc$height)[c(K-1, K)])           # height between K-1 and K merges
ddata  <- ggdendro::dendro_data(as.dendrogram(hc), type = "rectangle")
p_dend <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.2, color="gray30") +
  geom_hline(yintercept=h_cut, color=accent, linewidth=0.7, linetype="dashed") +
  annotate("text", x=Inf, y=h_cut, label=paste0("k = ",K," cut"),
           hjust=1.05, vjust=-0.5, size=3, color=accent, family=FONT_FAMILY) +
  labs(x=paste0("Proteins (n=", length(sel), ")"), y="Ward.D2 merge height") +
  theme_paper() +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())

diag_fig <- (p_elbow | p_sil) / p_dend +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")
png_cairo("clustering_diagnostic.png", 10, 8); print(diag_fig); dev.off()
cairo_pdf("clustering_diagnostic.pdf", width=10, height=8, family=FONT_FAMILY); print(diag_fig); dev.off()

# 7. per-cluster protein lists 

clusters_tbl <- tibble(uniprot = rownames(z), cluster = as.character(cl2)) |>
  left_join(sym_lookup |> select(uniprot, symbol, protein), by = "uniprot") |>
  left_join(tibble(uniprot = names(aov_q), omnibus_q = aov_q), by = "uniprot") |>
  arrange(cluster, omnibus_q)
write_csv(clusters_tbl, "clusters_all_proteins.csv")

dir.create("cluster_lists", showWarnings = FALSE)
for (cc in levels(cl2)){
  ids <- clusters_tbl$uniprot[clusters_tbl$cluster == cc]
  writeLines(ids, file.path("cluster_lists", paste0("cluster_", cc, ".txt")))
}
message(sprintf("exported %d proteins across %d cluster lists into cluster_lists/",
                nrow(clusters_tbl), nlevels(cl2)))

# 8. AP1 SpnA dependence 
# spna / spna_sig were computed in Section 4 (reused here for the focus heatmap).
if (length(spna_sig) >= 2) {
  z_sp <- zscore(expr[spna_sig, col_order, drop=FALSE])

  # right-side read-out: direction of the SpnA effect (KO-WT) at each timepoint
  spna_dir <- spna |> filter(uniprot %in% spna_sig) |>
    mutate(dir = ifelse(q < P_SIG & log2FC > 0, "Up",
                 ifelse(q < P_SIG & log2FC < 0, "Down", "n.s."))) |>
    select(uniprot, tp, dir) |> pivot_wider(names_from = tp, values_from = dir)
  spna_dir <- spna_dir[match(spna_sig, spna_dir$uniprot), ]
  # unicode Delta only in plain strings (not backtick names — Rscript rejects \u there)
  lab3  <- "\u0394SpnA 3h"; lab24 <- "\u0394SpnA 24h"
  rt <- data.frame(spna_dir$`3h`, spna_dir$`24h`, check.names = FALSE)
  names(rt) <- c(lab3, lab24)

  right_sp <- rowAnnotation(df = rt,
    col = setNames(list(dir_cols, dir_cols), c(lab3, lab24)),
    annotation_name_gp = gpar(fontsize=8), annotation_name_rot = 45,
    simple_anno_size = unit(3,"mm"), show_legend = c(TRUE, FALSE),
    annotation_legend_param = setNames(list(list(title = "\u0394SpnA effect (KO\u2212WT)",
      title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8))), lab3))

  zl_sp <- quantile(abs(z_sp), 0.99); col_sp <- colorRamp2(c(-zl_sp,0,zl_sp), divergent)
  ht_sp <- Heatmap(z_sp, name = "z-score", col = col_sp,
    top_annotation = top_anno, right_annotation = right_sp,
    cluster_columns = FALSE, column_split = meta_ord$strain, column_gap = unit(1,"mm"),
    show_column_names = FALSE,
    cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
    row_labels = spna_sig, row_names_gp = gpar(fontsize=7.5, fontface="italic"),
    column_title_gp = gpar(fontsize=9, fontface="bold"),
    heatmap_legend_param = list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8)))

  draw_sp <- function() draw(ht_sp, heatmap_legend_side="right",
                             annotation_legend_side="right", merge_legend=TRUE)
  png_cairo("heatmap_AP1-dSpnA_significant.png", 8.5, 8); draw_sp(); dev.off()
  cairo_pdf("heatmap_AP1-dSpnA_significant.pdf", width=8.5, height=8, family=FONT_FAMILY); draw_sp(); dev.off()

  spna |> filter(uniprot %in% spna_sig) |>
    left_join(sym_lookup |> select(uniprot, symbol, protein), by = "uniprot") |>
    select(uniprot, symbol, tp, log2FC, q) |> arrange(tp, q) |>
    write_csv("proteins_AP1-dSpnA_significant.csv")
  message(sprintf("AP1 \u0394SpnA focus heatmap written (%d proteins).", length(spna_sig)))
} else {
  message("AP1 \u0394SpnA focus heatmap skipped: fewer than 2 significant proteins.")
}

# 9. SpnA-responsive proteins vs the uninfected Control

if (length(spna_sig) >= 2) {

  gm <- function(s, t) rowMeans(expr[spna_sig, cf(s, t), drop = FALSE])   # group mean, log2
  build_db <- function(tpt)
    tibble(uniprot     = spna_sig,
           AP1_vs_ctrl = gm("AP1", tpt)             - gm("Control", tpt),
           KO_vs_ctrl  = gm("AP1\u0394SpnA", tpt)   - gm("Control", tpt))

  db3  <- build_db("3h")
  db24 <- build_db("24h")
  lev  <- db3$uniprot[order(db3$AP1_vs_ctrl)]        # lock row order to panel A
  db3  <- db3  |> mutate(uniprot = factor(uniprot, levels = lev))
  db24 <- db24 |> mutate(uniprot = factor(uniprot, levels = lev))

  db_cols <- c("AP1 (WT)"      = unname(strain_cols["AP1"]),           # red
               "AP1dSpnA (KO)" = unname(strain_cols["AP1\u0394SpnA"])) # light pink

  db_panel <- function(df, xlab, show_y) {
    long <- df |>
      pivot_longer(c(AP1_vs_ctrl, KO_vs_ctrl), names_to = "cond", values_to = "d") |>
      mutate(cond = factor(ifelse(cond == "AP1_vs_ctrl", "AP1 (WT)", "AP1dSpnA (KO)"),
                           levels = c("AP1 (WT)", "AP1dSpnA (KO)")))
    p <- ggplot(df, aes(y = uniprot)) +
      geom_vline(xintercept = 0, color = "gray55", linewidth = 0.5, linetype = "dotted") +
      geom_segment(aes(x = AP1_vs_ctrl, xend = KO_vs_ctrl, yend = uniprot),
                   color = "gray78", linewidth = 0.5) +
      geom_point(data = long, aes(x = d, color = cond), size = 2.4) +
      scale_color_manual(values = db_cols, name = NULL,
                         labels = c("AP1 (WT)", "AP1\u0394SpnA (KO)")) +
      scale_x_continuous(breaks = seq(-6, 6, 2)) +
      labs(x = xlab, y = NULL) +
      theme_paper() +
      theme(legend.position = "top")
    if (show_y) p + theme(axis.text.y = element_text(size = 7.5, face = "italic"))
    else        p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                          axis.line.y = element_blank())
  }

  pA <- db_panel(db3,  "log2 abundance vs Control (3 h)",  TRUE)
  pB <- db_panel(db24, "log2 abundance vs Control (24 h)", FALSE)
  db_comp <- (pA | pB) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "top",
          plot.tag = element_text(size = 13, face = "bold", family = FONT_FAMILY))

  ggsave("spna_vs_control_3h_24h.png", db_comp, width = 9.6, height = 6.8,
         dpi = 300, type = "cairo")
  cairo_pdf("spna_vs_control_3h_24h.pdf", width = 9.6, height = 6.8, family = FONT_FAMILY)
  print(db_comp); dev.off()
  message(sprintf("SpnA polarised dumbbell written (%d proteins; 3h WT-KO gap %.2f, 24h %.2f).",
                  length(spna_sig),
                  mean(abs(db3$AP1_vs_ctrl  - db3$KO_vs_ctrl)),
                  mean(abs(db24$AP1_vs_ctrl - db24$KO_vs_ctrl))))
} else {
  message("SpnA polarised dumbbell skipped: fewer than 2 significant proteins.")
}

