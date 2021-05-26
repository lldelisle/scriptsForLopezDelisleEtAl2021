options(stringsAsFactors = F)

# Dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("stringi")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# One of the dependency of ggpubr is only available for R=4.
if ( ! require("ggpubr") & ! rversionAbove(4)){
  packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  safelyLoadAPackageInCRANorBioconductor("ggpubr")
}

# Personal function
expected.values <- function(name, x){
  if(length(unique(round((x[-1] - x[-length(x)])/(x[2] - x[1]), 5))) > 1){
    stop("I did not expect non regular x")
  }
  delta.x <- x[2] - x[1]
  p <- strsplit(name, "_")[[1]]
  val <- rep(0, length(x))
  for (i in 1:(length(p) / 4)){
    distrib <- p[4 * (i - 1) + 1]
    amp <- as.numeric(p[4 * (i - 1) + 2])
    loc <- as.numeric(p[4 * (i - 1) + 3])
    scale <- as.numeric(p[4 * (i - 1) + 4])
    if (distrib == "uniform"){
      cur.val <- dunif(x, min = loc, max = scale + loc)
      val <- val + cur.val * amp / sum(cur.val)
    } else if (distrib == "gauss"){
      cur.val <- dnorm(x, mean = loc, sd = scale)
      val <- val + cur.val * amp / sum(cur.val)
    }
  }
  val[x==0] <- val[x==0] + (1 - sum(val))
  return(val / delta.x)
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

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "gene", "xmax", "group"))
my.table$i <- rownames(my.table)

# I assume the input is in the current directory
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$input), function(fn){
  df <- read.delim(list.files(pattern = fn, recursive = T, full.names = T)[1])
})
names(data) <- unique(my.table$input)
nb.per.group <- NULL
# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "gene", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "gene")]), unique(meta.data[, c("i", "gene")]), by = "i")
if(!all(temp$gene.x == temp$gene.y)){
  # The gene had "_" in its name:
  meta.data$gene <- my.table$gene[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("logevid.file", "gene")], 1, function(v){
    return(paste(strsplit(gsub("_pdf.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)
}
# Find the group name
meta.data$group <- apply(meta.data[, c("id", "i")], 1, function(v){
  group.name <- my.table$group[my.table$i == v[2]]
  if (! is.na(group.name) & group.name != ""){
    return(paste0(group.name, gsub(paste0("_", v[2], "$"), "", strsplit(v[1], paste0("_", group.name))[[1]][2])))
  } else {
    return("all")
  }
})
if (! is.null(nb.per.group)){
  meta.data$group <- factor(meta.data$group, levels = names(nb.per.group))
}

# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Process the input data:
norm.values <- cbind(log(1 + 10^4 * subset(data[[1]], select = as.character(unique(pdfs$gene))) / data[[1]]$nCount_RNA),
                     subset(data[[1]], select = unique(my.table$group)))
# Reshape
norm.values <- melt(norm.values, 
                    measure.vars = as.character(unique(pdfs$gene)),
                    variable_name = "gene")

# rbind for each group name + update the group name
full.norm.values <- NULL
for(group.name in unique(my.table$group)){
  temp.df <- norm.values[, c("gene", "value")]
  temp.df$group <- paste0(group.name, norm.values[, group.name])
  full.norm.values <- rbind(full.norm.values, temp.df)
}
# Add the values when all cells were used
comb.norm.values <- norm.values[, c("gene", "value")]
comb.norm.values$group <- "all"
full.norm.values <- rbind(full.norm.values, comb.norm.values)

# Put the levels in good order:
if (! is.null(nb.per.group)){
  full.norm.values$group <- factor(full.norm.values$group, levels = names(nb.per.group))
}

# Specify labels for panels
label.group <- c(expression(paste("FL ", italic("Pitx1"^"+/+"))), expression(paste("HL ", italic("Pitx1"^"+/+"))), expression(paste("HL ", italic("Pitx1"^"Pen-/Pen-"))))
names(label.group) <- c("limbtypeFLWT", "limbtypeHLWT", "limbtypeHLPendel")
color.group <- c("black", scales::hue_pal()(2))

label_as_notation <- function(x, ...) {
  list(lapply(unlist(x), function(X){
    switch(as.character(X),
           "limbtypeHLWT" = parse(text = "paste(\"HL \", italic(\"Pitx1\"^\"+/+\"))"),
           "limbtypeHLPendel" = parse(text = "paste(\"HL \", italic(\"Pitx1\"^\"Pen-/Pen-\"))"),
           "limbtypeFLWT" = parse(text = "paste(\"FL \", italic(\"Pitx1\"^\"+/+\"))"),
           X
    )
  }
  ))
}
# First plot all
my.df <- pdfs
g <- ggplot(my.df) +
  geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
  facet_grid(group ~ gene, scales='free',
             labeller = label_as_notation) +
  geom_density(data = subset(full.norm.values, gene %in% my.df$gene & group %in% my.df$group),
               aes(x = value, color = "Density from data"), fill = NA,
               trim = T, key_glyph = "path")  +
  geom_line(aes(x = x, y = mean, color = "baredSC")) +
  theme_classic()  +
  geom_blank(data = data.frame(x = rep(1, 2), y = 0.5, 
                               groupCol = c("Density from data", "baredSC")), 
             aes(y = y, color = groupCol, fill = groupCol)) +
  scale_fill_manual(name = "",
                    values = c('Density from data' = NA, 'baredSC' = "darkgreen")) +
  scale_color_manual(name = "",
                     values = c('Density from data' = "red", 'baredSC'="darkgreen")) +
  xlab(parse(text = paste0("paste(\"log(1 + \",10^4, italic(\"", unique(meta.data$gene), "\"),\")\")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 3)) +
  expand_limits(x=0)

ggsave(paste0(output.prefix, "_combinedModels_all.pdf"),
       width = 3 + 2 * length(unique(my.df$gene)), height = 6, limitsize = F)

g1 <- ggplot(subset(full.norm.values, gene %in% my.df$gene & group %in% my.df$group)) +
  geom_density(aes(x = value, color = group), fill = NA,
               trim = T, key_glyph = "path") +
  theme_classic() +
  xlab(parse(text = paste0("paste(\"log(1 + \",10^4, italic(\"", unique(meta.data$gene), "\"),\")\")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1)) +
  expand_limits(x=0) + 
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  ggtitle("Normalized expression\nfrom scRNA-seq")
g2 <- ggplot(my.df) +
  geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = group), alpha = 0.3, color = NA) +
  geom_line(aes(x = x, y = mean, color = group), key_glyph = "path") +
  theme_classic() +
  xlab(parse(text = paste0("paste(\"log(1 + \",10^4, italic(\"", unique(meta.data$gene), "\"),\")\")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1.5)) +
  expand_limits(x=0) + 
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  scale_fill_manual(name = "Genotype",
                    values = color.group,
                    breaks = names(label.group),
                    labels = label.group) +
  ggtitle("Normalized expression\nafter baredSC")

g.all <- ggarrange(g1, g2, labels = c("C", "D"),
                   legend = "bottom",
                   legend.grob = get_legend(g2, "bottom"))
ggsave(paste0(output.prefix, "_combinedModels_fig_Rouco_cd.pdf"), g.all,
       width = 6, height = 3.5)
