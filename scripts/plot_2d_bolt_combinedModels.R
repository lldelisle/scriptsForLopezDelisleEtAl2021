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
corr.from.pdf <- function(df){
  delta.x <- (max(df$x) - min(df$x)) / (length(unique(df$x)) - 1)
  delta.y <- (max(df$y) - min(df$y)) / (length(unique(df$y)) - 1)
  # Evaluate pdf1d:
  pdfx <- ddply(df, .(x), summarize,
                val = sum(value))
  pdfx$val <- pdfx$val * delta.y
  pdfy <- ddply(df, .(y), summarize,
                val = sum(value))
  pdfy$val <- pdfy$val * delta.x
  
  # Evaluate the correlation:
  mux <- sum(pdfx$x * pdfx$val) * delta.x
  muy <- sum(pdfy$y * pdfy$val) * delta.y
  pdfx$cx <- pdfx$x - mux
  pdfy$cy <- pdfy$y - muy
  varx <- sum(pdfx$cx ^ 2 * pdfx$val) * delta.x
  vary <- sum(pdfy$cy ^ 2 * pdfy$val) * delta.y
  df$cx <- pdfx$cx[match(df$x, pdfx$x)]
  df$cy <- pdfy$cy[match(df$y, pdfy$y)]
  covxy <- with(df, sum(cx * cy * value)) * delta.x * delta.y
  return(covxy / sqrt(varx * vary))
}

wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgenex\tgeney\txmax\tymax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "genex", "geney", "xmax", "ymax", "group"))
my.table$i <- rownames(my.table)

# I assume the input is in the current directory
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$input), function(fn){
  df <- read.delim(list.files(pattern = fn, recursive = T, full.names = T)[1])
})
names(data) <- unique(my.table$input)

my.table$genes <- paste0(my.table$genex, "VS", my.table$geney)

# Get the size of each group
nb.per.group <- NULL


# Get the input data:
input <- list()
for (i in unique(my.table$i)){
  genex <- gsub("-", ".", my.table$genex[my.table$i == i])
  if (!substr(genex, 1, 1) %in% c(letters, LETTERS)){
    genex <- paste0("X", genex)
  }
  geney <- gsub("-", ".", my.table$geney[my.table$i == i])
  if (!substr(geney, 1, 1) %in% c(letters, LETTERS)){
    geney <- paste0("X", geney)
  }
  group <- my.table$group[my.table$i == i]
  if (! is.na(group) && group != ""){
    input.data <- data[[my.table$input[my.table$i == i]]][, c(genex, geney, "nCount_RNA", group)]
    colnames(input.data) <- c("x", "y", "nCount_RNA", "group")
  } else {
    input.data <- data[[my.table$input[my.table$i == i]]][, c(genex, geney, "nCount_RNA")]
    colnames(input.data) <- c("x", "y", "nCount_RNA")
    input.data$group <- "all"
  }
  input[[i]] <- input.data
}

# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf2d_flat.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf2d_flat.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "genes", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "genes")]), unique(meta.data[, c("i", "genes")]), by = "i")
if(!all(temp$genes.x == temp$genes.y)){
  # The gene had "_" in its name:
  meta.data$genes <- my.table$genes[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("file", "genes")], 1, function(v){
    return(paste(strsplit(gsub("_pdf2d_flat.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$genes, "__", meta.data$info)
}
# Process genes
meta.data <- cbind(meta.data, matrix(unlist(strsplit(meta.data$genes, "VS")), ncol = 2, byrow = T))
colnames(meta.data)[-1:0 + ncol(meta.data)] <- c("genex", "geney")
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
  df$value <- df$mean
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)
pdfs$condition <- "baredSC"

# Compute the same pdfs from data post Poisson:

pdfs.input <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  i <- unique(meta.data$i[meta.data$id == my.id])
  group.name <- my.table$group[my.table$i == i]
  my.group <- unique(meta.data$group[meta.data$id == my.id])
  if (! is.na(group.name) & group.name != ""){
    data.per.cell <- subset(input[[i]], group == gsub(paste0("^", group.name), "", my.group))
  } else {
    data.per.cell <- input[[i]]
  }
  labels.x <- unique(sort(pdfs$x[pdfs$id == my.id]))
  if (length(labels.x) == 0){
    labels.x <- unique(sort(pdfs$x[pdfs$i == 1]))
  }
  xmax.val <- my.table[my.table$i == i, "xmax"]
  delta.x <- xmax.val / length(labels.x)
  data.per.cell$cut_x <- cut(log(1 + 10 ^ 4 * data.per.cell$x / data.per.cell$nCount_RNA), 
                             breaks = seq(from =  0, 
                                          to = xmax.val, 
                                          length.out = length(labels.x) + 1),
                             right = FALSE,
                             labels = labels.x)
  labels.y <- unique(sort(pdfs$y[pdfs$id == my.id]))
  if (length(labels.y) == 0){
    labels.y <- unique(sort(pdfs$y[pdfs$i == 1]))
  }
  ymax.val <- my.table[my.table$i == i, "ymax"]
  delta.y <- ymax.val / length(labels.y)
  data.per.cell$cut_y <- cut(log(1 + 10 ^ 4 * data.per.cell$y / data.per.cell$nCount_RNA), 
                             breaks = seq(from =  0, 
                                          to = ymax.val,  
                                          length.out = length(labels.y) + 1),
                             right = FALSE,
                             labels = labels.y)
  new.df <- with(data.per.cell, as.data.frame(table(cut_x, cut_y)))
  new.df$x <- as.numeric(as.character(new.df$cut_x))
  new.df$y <- as.numeric(as.character(new.df$cut_y))
  new.df$value <- new.df$Freq / sum(new.df$Freq) / delta.x / delta.y
  new.df$id <- my.id
  return(new.df[, c("x", "y", "value", "id")])
}))
pdfs.input$condition <- "input"

# Get the corr
meta.data$file.cor <- gsub("pdf2d_flat", "corr", meta.data$file)
corr <- do.call(rbind, lapply(meta.data$file.cor, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file.cor <- fn
  return(df)
}))
# Add the meta data
corr <- merge(corr, meta.data)
# A label is formatted:
corr$label <- paste0(round(corr$mean, 2), "[-", round(corr$mean, 2) - round(corr$low, 2), "]",
                     "^{+", round(corr$high, 2) - round(corr$mean, 2), "}")
# For the p-value, a superior value is given (still only mean + sd)
corr$p.label <- sapply(corr$pval + corr$error, function(v){paste0("p<", format(v, digits = 2))})
corr$condition <- "baredSC"
corr$condition <- factor(corr$condition, levels = c("input", "baredSC"))


my.cor <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  df <- subset(pdfs.input, id == my.id)
  cor.input <- corr.from.pdf(df)
  df <- subset(pdfs, id == my.id)
  cor.mcmc <- corr.from.pdf(df)
  return(data.frame(id = my.id, corr = c(cor.input, cor.mcmc), condition = c("input", "baredSC")))
}))
my.cor$condition <- factor(my.cor$condition, levels = c("input", "baredSC"))
my.cor <- merge(my.cor, unique(meta.data[, c("i", "group", "id")]))


# Get the proportion of each model
prop.model <- do.call(rbind, lapply(unique(meta.data$id), function(my.id){
  temp.df <- as.data.frame(matrix(na.omit(
    as.numeric(unlist(
      strsplit(grep("Using", readLines(file.path(directory, gsub("_pdf2d_flat.txt", ".log", meta.data$file[meta.data$id == my.id]))),
                    value = T),
               " ")
    ))), byrow = T, ncol = 2))
  samples.per.model <- ddply(temp.df, .(V2), summarise, min(V1))
  colnames(samples.per.model) <- c("ngauss", "samples")
  samples.per.model$ngauss <- samples.per.model$ngauss + 1
  samples.per.model$prop <- samples.per.model$samples / sum(samples.per.model$samples)
  samples.per.model$id <- my.id
  return(samples.per.model)
}))

prop.model <- merge(prop.model, unique(meta.data[, c("id", "i", "group")]))
prop.model$condition <- "baredSC"
prop.model$condition <- factor(prop.model$condition, levels = c("input", "baredSC"))


# ggplot(pdfs, aes(x, y)) +
#   geom_tile(aes(fill = log(1 + value))) +
#   facet_wrap(i + group ~ .) +
#   theme_classic() +
#   scale_fill_gradient(low="white", high="black") 

all.pdfs <- rbind(pdfs[, c("x", "y", "value", "id", "condition")], pdfs.input)
all.pdfs <- merge(all.pdfs, unique(meta.data[, c("i", "id", "group")]))
all.pdfs$condition <- factor(all.pdfs$condition, levels = c("input", "baredSC"))

all.pdfs_round <- all.pdfs
all.pdfs_round[, c("x", "y")] <- round(all.pdfs_round[, c("x", "y")], 5)

label.group <- paste("cluster", gsub("seurat_clusters", "", unique(all.pdfs_round$group)))
names(label.group) <- unique(all.pdfs_round$group)

g <- ggplot(all.pdfs_round, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  facet_grid(condition ~ group,
             labeller = labeller(group = label.group)) +
  # geom_text(data = prop.model, aes(x = 1 + 0.5 * ngauss, label = ngauss, size = prop),
  #           y =  max(all.pdfs_round$y) * 0.9,
  #           show.legend = F) +
  # geom_point(data = prop.model, aes(size = prop),
  #            x = 1, y = 1, colour = NA) +
  geom_text(data = corr,
            aes(label = label),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 0.9,
            size = 2.5, hjust = 0,
            parse = T) +
  geom_text(data = corr,
            aes(label = p.label),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 0.6,
            size = 2.5, hjust = 0) +
  geom_text(data = subset(my.cor, condition != "baredSC"),
            aes(label = round(corr, 2)),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 0.9,
            size = 2.5, hjust = 0) +
  theme_classic() +
  # scale_size("proportion\nof each model",
  #            guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
  scale_fill_gradient(low="white", high="black") +
  xlab(parse(text = paste0("paste(\"log(1 + \",10^4, italic(\"", unique(meta.data$genex), "\"),\")\")"))) +
  ylab(parse(text = paste0("paste(\"log(1 + \",10^4, italic(\"", unique(meta.data$geney), "\"),\")\")")))
ggsave(paste0(output.prefix, "_combinedModels.pdf"), g, width = 1.8 + 1.3 * length(unique(all.pdfs_round$id)), height = 3.5, limitsize = F)
