options(stringsAsFactors = F)

# Dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
library(devtools)
install_github("lldelisle/usefulLDfunctions", upgrade = "never")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("reticulate")
use_condaenv("baredSC", required = T)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# One of the dependency of ggpubr is only available for R=4.
if ( ! require("ggpubr") & ! rversionAbove(4)){
  packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  safelyLoadAPackageInCRANorBioconductor("ggpubr")
}

# Personnal functions
weight_gauss <- function(x, new_x, sigma){
  x <- unlist(x)
  dist <- (new_x %*% t(rep(1, length(x))) - rep(1, length(new_x)) %*% t(x)) ^ 2
  weight <- exp( - dist / (2 * sigma ^ 2))
  return(weight / rowSums(weight))
}
gauss_smooth <- function(x, y, new_x, sigma){
  return(weight_gauss(x, new_x, sigma) %*% y)
}
get_raw_count_matrix <- function(ad, h5ad_fn){
  data_ad <- ad$read_h5ad(h5ad_fn)
  # I just extract the raw counts
  raw_counts <- as.matrix(py_to_r(data_ad$T$X))
  rownames(raw_counts) <- py_to_r(data_ad$var_names$values)
  colnames(raw_counts) <- py_to_r(data_ad$obs_names$values)
  return(raw_counts)
}
# Change working directory:
setwd(commandArgs(TRUE)[1])
# Get the output directory
output.directory <- commandArgs(TRUE)[2]

# Store all raw values in list per experiment:
all.data <- list()

# I assume the input data are in the current directory

# For the control experiments, the data were already preprocessed
# I need to convert the format:
ad <- import("anndata", convert = FALSE)
# I just extract the raw counts
raw_counts <- get_raw_count_matrix(ad, "svensson_chromium_control.h5ad")
# I put separately each experiment
all.data[["Svensson2"]] <- as.matrix(raw_counts[, grepl(pattern = "20312_", x = colnames(raw_counts))])
all.data[["Svensson1"]] <- as.matrix(raw_counts[, grepl(pattern = "20311_", x = colnames(raw_counts))])
all.data[["Klein"]] <- get_raw_count_matrix(ad, "klein_indrops_control_GSM1599501.h5ad")
all.data[["Macosko"]] <- get_raw_count_matrix(ad, "macosko_dropseq_control.h5ad")
all.data[["Zheng"]] <- get_raw_count_matrix(ad, "zheng_gemcode_control.h5ad")


# Homogeneous cells:
# # I followed what is written in the supplementary methods of
# # Svensson, V. Droplet scRNA-seq is not zero-inflated. Nat Biotechnol 38, 147–150 (2020). https://doi.org/10.1038/s41587-019-0379-5
# # To separate NIH3T3 cells and HEK293T cells were first filtered to have at least 2,000 UMIs
# # and then assigned to HEK293T if they had 20 times more UMI’s
# # from human than from mouse, and vice versa for NIH3T3.
# # Other cells were discarded.
cells <- as.matrix(Read10X("filtered_feature_bc_matrix/"))
dim(cells)
# [1] 112137   5366
cells <- cells[, colSums(cells > 0) > 2000]
dim(cells)
# [1] 112137   4652
human_umis <- colSums(cells[grepl("hg19", rownames(cells)),])
mouse_umis <- colSums(cells[!grepl("hg19", rownames(cells)),])
all.data[["NIH3T3"]] <- cells[!grepl("hg19", rownames(cells)), mouse_umis > 20 * human_umis]
all.data[["HEK293T"]] <- cells[grepl("hg19", rownames(cells)), human_umis > 20 * mouse_umis]

# I remove all genes with no count:
all.data <- lapply(all.data, function(m){
  m[rowSums(m) > 0, ]
})

# We duplicate each data where we removed the ERCC:
all.data.noERCC <- lapply(all.data, function(mat){mat[!grepl("ERCC", rownames(mat)), ]})
names(all.data.noERCC) <- paste0(names(all.data), "_noERCC")

all.data.noERCC <- subsetByNamesOrIndices(all.data.noERCC, which(sapply(all.data.noERCC, nrow) != 0))

all.data <- c(all.data, all.data.noERCC)

# I measure the Xi = ki / Ni (for each cell the normalized expression)
norm.counts <- lapply(all.data, function(m){apply(m, 2, function(v){v / sum(v)})})

# I estimate the lambda  = mean gene expression by Mg = sum(ki) / sum(Ni)
gene_norm_means <- lapply(all.data, function(m){rowSums(m) / sum(m)})

# I get the Nis
tot.counts <- lapply(all.data, colSums)
# I check Nis:
counts.df <- data.frame(tot.counts = unlist(tot.counts),
                        expe = rep(names(tot.counts), sapply(tot.counts, length)))
g <- ggplot(counts.df, aes(x = "", y = tot.counts)) +
  geom_violin() +
  facet_grid(expe ~ ., scales = "free")
ggsave(file.path(output.directory, "Nis.pdf"), width = 5, height = 10)

# I evaluate the phi_poisson to compare to Svensson:
ddply(counts.df, .(expe), summarize,
      phi.poisson = round(var(tot.counts)/mean(tot.counts)^2, 3))
#                expe phi.poisson
# 1           HEK293T      0.2298
# 2    HEK293T_noERCC      0.2298
# 3             Klein      0.0451
# 4      Klein_noERCC      0.0454
# 5           Macosko      0.1151
# 6    Macosko_noERCC      0.1345
# 7            NIH3T3      0.1092
# 8     NIH3T3_noERCC      0.1092
# 9         Svensson1      0.0969
# 10 Svensson1_noERCC      0.1042
# 11        Svensson2      0.3787
# 12 Svensson2_noERCC      0.4247
# 13            Zheng      0.0426

# Select like in Ahlmann-Eltze et al. and check the phi_poisson
all.data.subset <- lapply(names(all.data), function(expe){
  threshold <- median(tot.counts[[expe]]) * c(1, 1.3)
  all.data[[expe]][, which(tot.counts[[expe]] < threshold[2] & tot.counts[[expe]] > threshold[1])]
})
names(all.data.subset) <- paste0(names(all.data), "_subset")
tot.counts.subset <- lapply(tot.counts, function(v){
  threshold <- median(v) * c(1, 1.3)
  return(v[which(v < threshold[2] & v > threshold[1])])
})
names(tot.counts.subset) <- paste0(names(tot.counts), "_subset")
counts.df.subset <- data.frame(tot.counts = unlist(tot.counts.subset),
                               expe = rep(names(tot.counts.subset), sapply(tot.counts.subset, length)))
ddply(counts.df.subset, .(expe), summarize,
      phi.poisson = round(var(tot.counts)/mean(tot.counts)^2, 3))
#                       expe phi.poisson
# 1    HEK293T_noERCC_subset       0.005
# 2           HEK293T_subset       0.005
# 3      Klein_noERCC_subset       0.005
# 4             Klein_subset       0.005
# 5    Macosko_noERCC_subset       0.005
# 6           Macosko_subset       0.005
# 7     NIH3T3_noERCC_subset       0.005
# 8            NIH3T3_subset       0.005
# 9  Svensson1_noERCC_subset       0.004
# 10        Svensson1_subset       0.004
# 11 Svensson2_noERCC_subset       0.005
# 12        Svensson2_subset       0.005
# 13            Zheng_subset       0.003


# The estimate of normalized variance of Xi
norm_var_est <- lapply(names(all.data), function(expe){
  rowSums(rep(1, nrow(all.data[[expe]])) %*% t(tot.counts[[expe]]) * 
            (norm.counts[[expe]] - gene_norm_means[[expe]]) ^ 2) /
    (length(tot.counts[[expe]]) - 1)
})
names(norm_var_est) <- names(all.data)

# Fit NB on norm_var_est:
slopeAVarNB <- lapply(names(all.data), function(expe){
  lm(norm_var_est[[expe]] ~  1 * gene_norm_means[[expe]] + I(gene_norm_means[[expe]]^2) + 0)$coefficients
})
names(slopeAVarNB) <- names(all.data)
phiNB <- lapply(names(all.data), function(expe){
  slopeAVarNB[[expe]] / 
    ((sum(tot.counts[[expe]])^2 - sum(tot.counts[[expe]]^2)) / 
       ((length(tot.counts[[expe]]) - 1) * sum(tot.counts[[expe]])))
})
names(phiNB) <- names(all.data)
slopeAVarNBlog <- lapply(names(all.data), function(expe){
  coef(nls(log(y) ~ log(x + a * x^2),
           data = data.frame(x = gene_norm_means[[expe]], y = norm_var_est[[expe]]),
           start=list(a = 0)))
})
names(slopeAVarNBlog) <- names(all.data)
phiNBlog <- lapply(names(all.data), function(expe){
  slopeAVarNBlog[[expe]] / 
    ((sum(tot.counts[[expe]])^2 - sum(tot.counts[[expe]]^2)) / 
       ((length(tot.counts[[expe]]) - 1) * sum(tot.counts[[expe]])))
})
names(phiNBlog) <- names(all.data)

# Prop of zeros
# Observed
prop.zeros <- lapply(all.data, function(m){rowMeans(m==0)})
# Predicted by Poisson
prop.zeros.poisson <- lapply(names(all.data), function(expe){
  rowMeans(exp(- gene_norm_means[[expe]] %*% t(tot.counts[[expe]])))
})
# Predicted by NB
prop.zeros.NB <- lapply(names(all.data), function(expe){
  rowMeans((1 / (1 + gene_norm_means[[expe]] %*% t(tot.counts[[expe]]) * phiNB[[expe]])) ^ (1 / phiNB[[expe]]))
})
# Predicted by NB fitted on log
prop.zeros.NBlog <- lapply(names(all.data), function(expe){
  rowMeans((1 / (1 + gene_norm_means[[expe]] %*% t(tot.counts[[expe]]) * phiNBlog[[expe]])) ^ (1 / phiNBlog[[expe]]))
})


# Like Svensson:
k_means <- lapply(all.data, rowMeans)
k_var <- lapply(all.data, function(m){apply(m, 1, var)})

# Fit NB on k_var:
slopeKVarNB <- lapply(names(all.data), function(expe){
  lm(k_var[[expe]] ~  1 * k_means[[expe]] + I(k_means[[expe]]^2) + 0)$coefficients
})
names(slopeKVarNB) <- names(all.data)

# I go gene centrist:
big.df <- data.frame(gene_means = unlist(gene_norm_means),
                     expe = rep(names(all.data), sapply(all.data, nrow)),
                     gene_name = unlist(lapply(all.data, rownames)),
                     norm_var_est = unlist(norm_var_est),
                     prop_zeros = unlist(prop.zeros),
                     prop_zeros_poisson = unlist(prop.zeros.poisson),
                     prop_zeros_NB = unlist(prop.zeros.NB),
                     prop_zeros_NBlog = unlist(prop.zeros.NBlog),
                     slopeAVarNB = rep(unlist(slopeAVarNB), sapply(all.data, nrow)),
                     phiNB = rep(unlist(phiNB), sapply(all.data, nrow)),
                     slopeAVarNBlog = rep(unlist(slopeAVarNBlog), sapply(all.data, nrow)),
                     phiNBlog = rep(unlist(phiNBlog), sapply(all.data, nrow)),
                     k_means = unlist(k_means),
                     k_var = unlist(k_var),
                     slopeKVarNB = rep(unlist(slopeKVarNB), sapply(all.data, nrow))
)
# The norm_var_est when it follows a NB:
big.df$norm_var_est_NB <-  big.df$gene_means + big.df$slopeAVarNB * big.df$gene_means ^ 2
# The norm_var_est when it follows a NB with NB fitted in log:
big.df$norm_var_est_NBlog <-  big.df$gene_means + big.df$slopeAVarNBlog * big.df$gene_means ^ 2
# The k_var when it follows a NB (like in Svensson):
big.df$k_var_NB <-  big.df$k_means + big.df$slopeKVarNB * big.df$k_means ^ 2


big.df$expe.name <- gsub("_noERCC", "", big.df$expe)
big.df$expe_type <- "Control"
big.df$expe_type[big.df$expe.name %in% c("HEK293T", "NIH3T3")] <- "Cells"

write.table(big.df, file.path(output.directory, "big.df.txt"), row.names = F,
            quote = F, sep = "\t")

big.df$expe_type <- factor(big.df$expe_type, levels = c("Control", "Cells"))
big.df <- ddply(big.df, .(expe), mutate,
                ngenes = length(gene_name))
# There are a lot of genes. It can be misleading on the graph
# I prefer to display a trend line with error bars
# I apply a kde
sigma <- 0.1
n.points <- 1000

smooth <- do.call(rbind, lapply(unique(big.df$expe), function(my.expe){
  my.df <- subset(big.df, expe == my.expe)
  df <- data.frame(logx = unique(quantile(log10(my.df$gene_means), seq(0, 1, length.out = n.points))))
  weight <- with(subset(my.df, gene_means > 0),
                 weight_gauss(log10(gene_means), df$logx, sigma))
  df$logy <- with(subset(my.df, gene_means > 0),
               weight %*% log10(norm_var_est))
  df$sigmalog <- with(subset(my.df, gene_means > 0),
                   sqrt(weight %*% (log10(norm_var_est)^2) - df$logy^2))
  df$min.logx <- min(log10(my.df$gene_means))
  df$expe = my.expe
  return(df)
})
)
smooth$x <- 10^smooth$logx
smooth$y <- 10^smooth$logy
smooth$lowy <- 10^(smooth$logy - smooth$sigmalog)
smooth$highy <- 10^(smooth$logy + smooth$sigmalog)
smooth <- merge(smooth, unique(big.df[, c("expe", "expe_type", "expe.name")]))

# Labeller:
my.labeller <- c("Svensson et al. (1)", "Svensson et al. (2)", "Klein et al.", "Macosko et al.", "Zheng et al.",
                 "NIH3T3", "HEK293T")
names(my.labeller) <- c("Svensson1", "Svensson2", "Klein", "Macosko", "Zheng",
                        "NIH3T3", "HEK293T")

# Reproduce Svensson analysis:
g1 <- ggplot(subset(big.df, !grepl("noERCC", expe)), aes(x = k_means)) +
  geom_point(aes(y = k_var, alpha = 300 / ngenes), show.legend = F) +
  geom_abline(aes(intercept = 0, slope = 1, color = "Poisson"), show.legend = F) +
  geom_line(aes(y = k_var_NB, colour = "NB")) +
  facet_wrap(expe_type ~ expe, scales = 'free', nrow = 2, ncol = 5,
             labeller = labeller(expe = my.labeller)) +
  scale_x_continuous(trans='log10')  +
  scale_y_continuous(trans='log10') +
  xlab("Means of k") +
  ylab("Variances of k") +
  theme_classic() +
  theme(legend.position="bottom",
        text = element_text(size = 15)) +
  labs(color = "") +
  scale_alpha_identity()

ggsave(file.path(output.directory, "figSvensson.pdf"), width = 12, height = 6)
ggsave(file.path(output.directory, "figSvensson.png"), width = 12, height = 6,
       dpi = 500)


# I plot the variance estimate as function of the mean expression
# With ERCC included
ggplot.list <- list(list(), list())
names(ggplot.list) <- levels(big.df$expe_type)

for (my.expe in names(my.labeller)){
  sub.df <- subset(big.df, expe == my.expe)
  ggplot.list[[sub.df$expe_type[1]]][[my.expe]] <- ggplot(sub.df, aes(x = gene_means)) +
    geom_ribbon(data = subset(smooth, expe == my.expe),
                aes(ymin = lowy, ymax = highy, x = x, fill = "data"),
                alpha = 0.2) +
    geom_point(aes(y = norm_var_est), alpha = 0.1) +
    geom_abline(aes(intercept = 0, slope = 1, color = "Poisson"), show.legend = F) +
    geom_line(aes(y = norm_var_est_NBlog, colour = "NB")) +
    # geom_line(aes(y = norm_var_est_NB, colour = "NB")) +
    geom_line(data = subset(smooth, expe == my.expe),
              aes(x, y, colour = "data")) +
    facet_wrap(. ~ expe.name,
               labeller = labeller(expe.name = my.labeller)) +
    scale_x_continuous(trans='log10')  +
    scale_y_continuous(trans='log10') +
    xlab(expression(italic(M[g]))) +
    ylab(expression(italic(tilde(V)[g]))) +
    theme_classic() +
    theme(text = element_text(size = 15)) +
    geom_blank(data = data.frame(gene_means = rep(1e-6, 3), y = 1e-6, 
                                 groupCol = c("data", "Poisson", "NB")), 
               aes(y = y, color = groupCol, fill = groupCol)) +
    scale_fill_manual(name = "",
                      values = c('data' = "yellow",
                                 "Poisson" = "white",
                                 "NB" = "white")) +
    scale_color_discrete("",
                         labels = c('data', "Poisson", "NB"))
}
g <- ggarrange(ggarrange(plotlist = ggplot.list[["Control"]], legend = "none",
                         nrow = 1, labels = paste(c("A", rep("", length(ggplot.list[["Control"]]) - 1)),
                                                  letters[1:length(ggplot.list[["Control"]])])),
          ggarrange(
            ggarrange(plotlist = ggplot.list[["Cells"]], legend = "none",
                    nrow = 1, labels = paste(c("B", rep("", length(ggplot.list[["Cells"]]) - 1)),
                                             letters[1:length(ggplot.list[["Cells"]])])),
            get_legend(ggplot.list[["Cells"]][[1]]),
            nrow = 1 , widths = c(2, length(ggplot.list[["Control"]]) - 2)
            ),
          nrow = 2)
ggsave(file.path(output.directory, "fig1_withERCC.png"), width = 15, height = 6, dpi = 500, bg = "white")
ggsave(file.path(output.directory, "fig1_withERCC.pdf"), width = 15, height = 6)

# Same without ERCC
ggplot.list <- list(list(), list())
names(ggplot.list) <- levels(big.df$expe_type)

for (my.expe in intersect(paste0(names(my.labeller), "_noERCC"), unique(big.df$expe))){
  sub.df <- subset(big.df, expe == my.expe)
  ggplot.list[[sub.df$expe_type[1]]][[my.expe]] <- ggplot(sub.df, aes(x = gene_means)) +
    geom_ribbon(data = subset(smooth, expe == my.expe),
                aes(ymin = lowy, ymax = highy, x = x, fill = "data"),
                alpha = 0.2) +
    geom_point(aes(y = norm_var_est), alpha = 0.1) +
    geom_abline(aes(intercept = 0, slope = 1, color = "Poisson"), show.legend = F) +
    geom_line(aes(y = norm_var_est_NBlog, colour = "NB")) +
    # geom_line(aes(y = norm_var_est_NB, colour = "NB")) +
    geom_line(data = subset(smooth, expe == my.expe),
              aes(x, y, colour = "data")) +
    facet_wrap(. ~ expe.name,
               labeller = labeller(expe.name = my.labeller)) +
    scale_x_continuous(trans='log10')  +
    scale_y_continuous(trans='log10') +
    xlab(expression(italic(M[g]))) +
    ylab(expression(italic(tilde(V)[g]))) +
    theme_classic() +
    theme(text = element_text(size = 13)) +
    geom_blank(data = data.frame(gene_means = rep(1e-6, 3), y = 1e-6, 
                                 groupCol = c("data", "Poisson", "NB")), 
               aes(y = y, color = groupCol, fill = groupCol)) +
    scale_fill_manual(name = "",
                      values = c('data' = "yellow",
                                 "Poisson" = "white",
                                 "NB" = "white")) +
    scale_color_discrete("",
                         labels = c('data', "Poisson", "NB"))
}
g <- ggarrange(ggarrange(plotlist = ggplot.list[["Control"]], legend = "none",
                         nrow = 1, labels = paste(c("A", rep("", length(ggplot.list[["Control"]]) - 1)),
                                                  letters[1:length(ggplot.list[["Control"]])])),
               ggarrange(
                 ggarrange(plotlist = ggplot.list[["Cells"]], legend = "none",
                           nrow = 1, labels = paste(c("B", rep("", length(ggplot.list[["Cells"]]) - 1)),
                                                    letters[1:length(ggplot.list[["Cells"]])])),
                 get_legend(ggplot.list[["Cells"]][[1]]),
                 nrow = 1 , widths = c(2, length(ggplot.list[["Control"]]) - 2)
               ),
               nrow = 2)
ggsave(file.path(output.directory, "fig1_noERCC.png"), width = 12, height = 6, dpi = 500, bg = "white")
ggsave(file.path(output.directory, "fig1_noERCC.pdf"), width = 12, height = 6)

# We plot the figure with colors for different gene type
ggplot.list <- list(list(), list())
names(ggplot.list) <- levels(big.df$expe_type)

for (my.expe in names(my.labeller)){
  sub.df <- subset(big.df, expe == my.expe)
  ggplot.list[[sub.df$expe_type[1]]][[my.expe]] <- ggplot(sub.df, aes(x = gene_means)) +
    geom_point(aes(y = norm_var_est, color = grepl("ERCC", gene_name)), alpha = 0.1) +
    facet_wrap(. ~ expe.name,
               labeller = labeller(expe.name = my.labeller)) +
    scale_x_continuous(trans='log10')  +
    scale_y_continuous(trans='log10') +
    xlab(expression(italic(M[g]))) +
    ylab(expression(italic(tilde(V)[g]))) +
    scale_color_manual(name = "gene type", labels = c('TRUE' = 'ERCC', 'FALSE' = 'gene'),
                       values = c('TRUE' = "blue", 'FALSE' = "red")) +
    guides(colour = guide_legend(override.aes = list(alpha = 0.5))) +
    theme_classic() +
    theme(text = element_text(size = 15))
}
g <- ggarrange(ggarrange(plotlist = ggplot.list[["Control"]], legend = "none",
                         nrow = 1, labels = paste(c("A", rep("", length(ggplot.list[["Control"]]) - 1)),
                                                  letters[1:length(ggplot.list[["Control"]])])),
               ggarrange(
                 ggarrange(plotlist = ggplot.list[["Cells"]], legend = "none",
                           nrow = 1, labels = paste(c("B", rep("", length(ggplot.list[["Cells"]]) - 1)),
                                                    letters[1:length(ggplot.list[["Cells"]])])),
                 get_legend(ggplot.list[["Control"]][[1]]),
                 nrow = 1 , widths = c(2, length(ggplot.list[["Control"]]) - 2)
               ),
               nrow = 2)
ggsave(file.path(output.directory, "figS1_geneColor.png"), width = 15, height = 6, dpi = 500, bg = "white")
ggsave(file.path(output.directory, "figS1_geneColor.pdf"), width = 15, height = 6)

# In order to generate data with realistic Nis, I take nih3t3 experiment:
nih3t3 <- CreateSeuratObject(counts = all.data[["NIH3T3"]])
set.seed(1)
# I create 1 subset of 300, another one of 800, another one with all
# I create 1 subset of 300, another one of 800, another one with all
df <- data.frame(nCount_RNA = nih3t3$nCount_RNA,
                 group = sample(rep(c('group1', 'group2', 'group3'),
                                    c(300, 500, length(nih3t3$nCount_RNA) - 800)),
                                length(nih3t3$nCount_RNA))
)
g <- ggplot(df, aes(x = group, y = nCount_RNA)) +
  geom_violin() +
  geom_jitter(size = .1)
ggsave(file.path(output.directory, "figS1D.pdf"), width = 7, height = 5)
write.table(df,
            file = file.path(output.directory, "nih3t3_nRNA.txt"),
            quote = F, sep = "\t", row.names = F)
