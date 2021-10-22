options(stringsAsFactors = F)

# Dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# One of the dependency of ggpubr is only available for R=4.
if ( ! require("ggpubr") & ! rversionAbove(4)){
  packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  safelyLoadAPackageInCRANorBioconductor("ggpubr")
}

wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgene\txmax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]

setwd(wd)

# Read the input table
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "gene", "xmax", "group"))
my.table$i <- 1:nrow(my.table)
# Extract from input the number of cells
my.table$ncells <- as.numeric(sapply(my.table$full.path.input, function(s){tail(unlist(strsplit(gsub(".txt$", "", s), "_")), 1)}))

# This part is specific to data generated:
# Process in detail the generation
meta.generation <- data.frame(gene = unique(my.table$gene))
meta.generation$label <- gsub("_gauss", "\n+ gauss", gsub("_uniform", "\n+ uniform", meta.generation$gene))
meta.generation$nb <- 1 + lengths(regmatches(meta.generation$label, gregexpr("\n", meta.generation$label)))
meta.generation$split <- strsplit(meta.generation$gene, "_")
meta.generation$first.model <- sapply(meta.generation$split, function(v){v[1]})
meta.generation$locs <- sapply(meta.generation$split, function(v){v[seq(3, length(v), 4)]})
meta.generation$min.loc <- sapply(meta.generation$locs, min)
meta.generation$scales <- sapply(meta.generation$split, function(v){v[seq(4, length(v), 4)]})
meta.generation$min.scales <- sapply(meta.generation$scales, min)
meta.generation$prop.zero <- sapply(meta.generation$split, function(v){1 - sum(as.numeric(v[seq(2, length(v), 4)]))})
meta.generation$pretty <- apply(meta.generation, 1, function(v){
  split.string <- v['split'][[1]]
  pretty.str <- ""
  for (i in 1:(length(split.string) / 4)){
    if(split.string[4 * (i-1) + 1] == "gauss"){
      first.letter <- "N"
    } else {
      first.letter <- "U"
    }
    if (i > 1){
      # pretty.str <- paste0(pretty.str, " + ")
      pretty.str <- paste0(pretty.str, "\n")
    }
    pretty.str <- paste0(pretty.str, first.letter, "(", split.string[4 * (i-1) + 3], ",", split.string[4 * (i-1) + 4], ")")
  }
  if (v['prop.zero'] > 0){
    # pretty.str <- paste0(pretty.str, " + ", as.numeric(v['prop.zero']) * 100, "% of 0")
    pretty.str <- paste0(pretty.str, "\n", as.numeric(v['prop.zero']) * 100, "% of 0")
  }
  return(pretty.str)
})
meta.generation <- meta.generation[order(meta.generation$nb, meta.generation$first.model, meta.generation$min.loc, meta.generation$min.scale), ]
my.table$gene <- factor(my.table$gene, levels = meta.generation$gene)

label.generation <- meta.generation$pretty
names(label.generation) <- meta.generation$gene

label.run <- c("Only the good\nvalue of m", "m=1, 2, 3 and 4\nin parallel\nThen, merge of 4")
names(label.run) <- c("goodm", "default")

# Collect info on run
full.df <- NULL
files <- list.files(directory, pattern = ".log", full.names = T)
for (file in files){
  temp.df <- read.delim(file, sep = "|")
  temp.df.cleaned <- subset(temp.df, grepl("batch", JobID))
  if (nrow(temp.df.cleaned) > 0){
    temp.df.cleaned$i <- as.numeric(sapply(temp.df.cleaned$JobID, function(s){strsplit(gsub(".batch", "", s), "_")[[1]][2]}))
    temp.df.cleaned$cond <- strsplit(basename(file), "__")[[1]][1]
    full.df <- rbind(full.df, temp.df.cleaned)
  }
}
full.df <- merge(full.df, my.table, all.x = T)
# Get time in minutes from (d-)HH:MM:SS
full.df$time.min <- sapply(full.df$Elapsed, function(s){
  v1 <- strsplit(s, "-")[[1]]
  if (length(v1) == 1){
    d <- 0
  } else {
    d <- as.numeric(v1[1])
    v1 <- v1[2]
  }
  v <- as.numeric(strsplit(v1, ":")[[1]])
  return(d * 24  * 60 + v[1] * 60 + v[2] + v[3] / 60)
})

# Attribute a unique ID
full.df$id <- apply(full.df[, c("cond", "gene", "ncells")], 1, paste, collapse = "_")

# Do not keep run which went out of time
full.df <- subset(full.df, State != "CANCELLED")

full.df <- ddply(full.df, .(id), mutate, nb.rep = length(id))
full.df$nb.rep <- factor(full.df$nb.rep, levels = sort(unique(full.df$nb.rep)))

meta.cond <- data.frame(cond = unique(full.df$cond))
meta.cond$split <- strsplit(meta.cond$cond, "_")
meta.cond$core.per.task <- sapply(meta.cond$split, function(v){as.numeric(strsplit(v[3], "core")[[1]][1])})
meta.cond$core.per.task <- factor(meta.cond$core.per.task, levels = sort(unique(meta.cond$core.per.task)))
meta.cond$run <- sapply(meta.cond$split, "[[", 2)

full.df <- merge(full.df, meta.cond)

# Get the slope if it is linear based on the 125000 cells:
slopes <- ddply(full.df, .(core.per.task, run, gene), summarise,
                slope = mean(time.min[ncells == 125000]) / 125000)

# Choose breaks
my.breaks <- c(1, 3, 10, 30, 60, 180, 5*60, 10* 60, 20 * 60)
my.breaks.labels <- sapply(my.breaks, function(t){paste(formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
                                                        ,formatC(t %% 60, width = 2, format = "d", flag = "0")
                                                        ,sep = ":")
})

sub.full.df <- subset(full.df, nb.rep == 3 & gene != "gauss_0.25_0.75_0.25")
sub.slopes <- subset(slopes, run == "goodm" & gene != "gauss_0.25_0.75_0.25")

sub.full.df$run <- factor(sub.full.df$run, levels = names(label.run))
sub.slopes$run <- factor(sub.slopes$run, levels = names(label.run))

# Statistics at 25'000 cells:
summary(sub.full.df$time.min[sub.full.df$ncells == 25000 & sub.full.df$run == "goodm" & sub.full.df$core.per.task == 1])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.033   2.725   2.892   2.858   3.237   3.400 
mean(sub.full.df$time.min[sub.full.df$ncells == 25000 & sub.full.df$run == "goodm" & sub.full.df$core.per.task == 4])
# [1] 1.547222
mean(sub.full.df$MaxRSS[sub.full.df$ncells == 25000 & sub.full.df$run == "goodm"])
# [1] 335973205
mean(sub.full.df$time.min[sub.full.df$ncells == 25000 & sub.full.df$run == "default"])
# [1] 27.84583
mean(sub.full.df$MaxRSS[sub.full.df$ncells == 25000 & sub.full.df$run == "default"])
# [1] 1474195797

# Plot the time in log/log scale
g1 <- ggplot(sub.full.df, aes(x = ncells, y = time.min, color = core.per.task)) +
  stat_summary(data = sub.full.df, aes(size = "average of 3 runs"), fun = "mean", geom = "point") +
  stat_summary(data = sub.full.df, fun = "mean", geom = "line") +
  geom_point(aes(size = "individual run")) +
  facet_grid(run ~ gene,
             labeller = labeller(gene = label.generation,
                                 run = label.run)) +
  scale_x_log10("Number of cells",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous("Time (HH:MM)",
                     trans='log10',
                     breaks = my.breaks,
                     labels = my.breaks.labels) +
  scale_size_manual("Value",
                    values = c('individual run' = 1,
                               'average of 3 runs' = 3)) +
  scale_color_discrete("Core per task") +
  geom_abline(data = sub.slopes, aes(slope = 1, intercept = log10(slope)),
              lty = 2) +
  theme_classic() + 
  theme(
    panel.grid.major = element_line(colour='grey90', size = 0.2)
  )

# Plot the memory in log/log scale
g2 <- ggplot(sub.full.df, aes(x = ncells, y = MaxRSS / 2^30, color = core.per.task)) +
  stat_summary(data = sub.full.df, aes(size = "average of 3 runs"), fun = "mean", geom = "point") +
  stat_summary(data = sub.full.df, fun = "mean", geom = "line") +
  geom_point(aes(size = "individual run")) +
  facet_grid(run ~ gene,
             labeller = labeller(gene = label.generation,
                                 run = label.run)) +
  scale_x_log10("Number of cells",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous("Memory used (GB)",
                     trans='log10') +
  scale_size_manual("Value",
                    values = c('individual run' = 1,
                               'average of 3 runs' = 3)) +
  scale_color_discrete("Core per task") +
  theme_classic() + 
  theme(
    panel.grid.major = element_line(colour='grey90', size = 0.2)
  )

g <- ggarrange(g1, g2, labels = c("A", "B"), common.legend = T, legend = "bottom")
ggsave(paste0(output.prefix, "_main.pdf"), height = 5, width = 10)

# Get all log files to know the number of samples
logs.files <- list.files(directory, pattern = "^[1-4]gauss.*.log$", full.names = T, recursive = T)
my.df <- do.call(rbind, lapply(logs.files, function(fn){
  i <- as.numeric(tail(strsplit(gsub(".log$", "", basename(fn)), "_")[[1]], 1))
  nsampconsidered <- as.numeric(strsplit(tail(grep("Considering the last", readLines(fn), value = T), 1), " ")[[1]][4])
  nsamp <- (nsampconsidered - 1) * 4 / 3
  ngauss <- as.numeric(substr(basename(fn), 1, 1))
  directory <- tail(strsplit(fn, "/")[[1]], 2)[1]
  directory.s <- strsplit(directory, "core(s)?_")[[1]]
  seed <- ifelse(length(directory.s) == 2, yes = directory.s[2], no = "1")
  return(c(i, nsamp, ngauss, seed, fn))
}))
my.df <- as.data.frame(my.df)
colnames(my.df) <- c("i", "nsamp", "ngauss", "seed", "fn")
my.df$ngauss <- factor(my.df$ngauss)
my.df <- merge(my.df, my.table)
# We remove when the same seed was used twice
# For example, with 1 core and 4 cores
# Or with goodm and default
my.df.u <- unique(my.df[, setdiff(colnames(my.df), "fn")])

my.df.u <- ddply(my.df.u, .(i, ngauss),
                 mutate, nb.rep = length(i))

expected <- data.frame(gene = unique(my.df.u$gene),
                       expected = sapply(strsplit(label.generation[unique(my.df.u$gene)], "\n"), 
                                         length))

ggplot(subset(my.df.u, nb.rep >= 3 & gene != "gauss_0.25_0.75_0.25"), aes(x = ngauss, y = nsamp, fill = ngauss, color = ngauss)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 3) +
  facet_grid(ncells ~ gene,
             labeller = labeller(gene = label.generation)) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(colour='grey90', size = 0.2)
  ) +
  xlab("Number of Gaussians in the model") +
  ylab("Number of samples in MCMC at convergence") +
  geom_vline(data = subset(expected, gene != "gauss_0.25_0.75_0.25"),
             aes(xintercept = expected,
                 lty = "expected"), alpha = 0.2) +
  scale_color_discrete("m\n(Number of\nGaussians\nin the model)") +
  scale_fill_discrete("m\n(Number of\nGaussians\nin the model)") +
  scale_linetype_manual("",
                        values = c('expected' = 2),
                        label = c('Number of\nGaussians\nin the real\ndistribution'))

ggsave(paste0(output.prefix, "_sup.pdf"), height = 5.5, width = 6)
