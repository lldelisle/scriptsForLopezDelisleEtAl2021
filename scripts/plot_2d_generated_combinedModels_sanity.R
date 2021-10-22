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
safelyLoadAPackageInCRANorBioconductor("ggh4x")

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

pdf.2d.norm <- function(df, mux, muy, sigx, sigy, corr){
  delta.x <- (max(df$x) - min(df$x)) / (length(unique(df$x)) - 1)
  delta.y <- (max(df$y) - min(df$y)) / (length(unique(df$y)) - 1)
  df$ux <- (df$x - mux) / sigx
  df$uy <- (df$y - muy) / sigy
  df$lognorm <- -0.5 * (df$ux * df$ux + df$uy * df$uy - 2 * corr * df$ux * df$uy) / (1 - corr * corr)
  df$lognorm.rescale <- df$lognorm - max(df$lognorm) # To avoid overflow
  df$value.rel <- exp(df$lognorm.rescale)
  df$value <- df$value.rel / (sum(df$value.rel) * delta.x * delta.y)
  return(df[, c("x", "y", "value")])
}

pdf.2d.post <- function(df.xy, df.cells){
  # df.xy contains the grid x,y
  # df.cells contains mux, muy, sigx, sigy
  df.all.cells <- apply(df.cells[, c("mux", "muy", "sigx", "sigy")], 1, function(v){pdf.2d.norm(df.xy, v["mux"], v["muy"], v["sigx"], v["sigy"], 0)$value})
  return(data.frame(x = df.xy$x, y = df.xy$y, value = apply(df.all.cells, 1, mean)))
}

wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgenex\tgeney\txmin\txmax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]
# The sanity results are in the directory
directory.sanity <- commandArgs(TRUE)[5]

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "genex", "geney", "xmin", "xmax", "group"))
my.table$i <- rownames(my.table)

# I assume there are not 2 tables with same basename
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$full.path.input), function(fn){
  df <- read.delim(fn)
})
names(data) <- unique(my.table$input)

my.table$genes <- paste0(my.table$genex, "VS", my.table$geney)

# Get the size of each group
if (length(unique(my.table$input)) == 1 & length(setdiff(unique(my.table$group), "")) == 1){
  nb.per.group <- table(data[[1]][, setdiff(unique(my.table$group), "")])
  names(nb.per.group) <- paste0(setdiff(unique(my.table$group), ""), names(nb.per.group))
  nb.per.group["all"] <- sum(nb.per.group)
} else {
  nb.per.group <- NULL
}

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

# Get the input data before Poisson:
beforeP <- list()
for (i in unique(my.table$i)){
  genex <- paste0(gsub("-", ".", my.table$genex[my.table$i == i]), "_expression")
  if (!substr(genex, 1, 1) %in% c(letters, LETTERS)){
    genex <- paste0("X", genex)
  }
  geney <- paste0(gsub("-", ".", my.table$geney[my.table$i == i]), "_expression")
  if (!substr(geney, 1, 1) %in% c(letters, LETTERS)){
    geney <- paste0("X", geney)
  }
  group <- my.table$group[my.table$i == i]
  if (! is.na(group) && group != ""){
    input.data <- data[[my.table$input[my.table$i == i]]][, c(genex, geney, group)]
    colnames(input.data) <- c("x", "y", "group")
  } else {
    input.data <- data[[my.table$input[my.table$i == i]]][, c(genex, geney)]
    colnames(input.data) <- c("x", "y")
    input.data$group <- "all"
  }
  beforeP[[i]] <- input.data
}


# Get Sanity results
m.sanity <- list()
m.sanity.e <- list()
for (cur.input in unique(my.table$input)){
  sub.directory.sanity <- gsub(".txt$", "", cur.input)
  m.sanity[[cur.input]] <- as.matrix(read.delim(file.path(directory.sanity, sub.directory.sanity, "log_transcription_quotients.txt"), row.names = 1))
  m.sanity.e[[cur.input]] <- as.matrix(read.delim(file.path(directory.sanity, sub.directory.sanity, "ltq_error_bars.txt"), row.names = 1))
}
sanity.res <- list()
for (i in unique(my.table$i)){
  genex <- my.table$genex[my.table$i == i]
  geney <- my.table$geney[my.table$i == i]
  cur.input <- my.table$input[my.table$i == i]
  group <- my.table$group[my.table$i == i]
  if (! is.na(group) && group != ""){
    return(NULL)
  } else {
    input.data <- as.data.frame(cbind(t(m.sanity[[cur.input]])[, c(genex, geney)],
                                      t(m.sanity.e[[cur.input]])[, c(genex, geney)]))
    colnames(input.data) <- c("mux", "muy", "sigx", "sigy")
    input.data$group <- "all"
  }
  sanity.res[[i]] <- input.data
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

# Compute the same pdfs from data pre and post Poisson:
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
  xmin.val <- my.table[my.table$i == i, "xmin"]
  xmax.val <- my.table[my.table$i == i, "xmax"]
  delta.x <- (xmax.val - xmin.val) / length(labels.x)
  x.vals <- log(data.per.cell$x / data.per.cell$nCount_RNA)
  x.vals[is.infinite(x.vals)] <- xmin.val
  data.per.cell$cut_x <- cut(x.vals, 
                             breaks = seq(from = xmin.val, 
                                          to = xmax.val, 
                                          length.out = length(labels.x) + 1),
                             right = FALSE,
                             labels = labels.x)
  labels.y <- unique(sort(pdfs$y[pdfs$id == my.id]))
  if (length(labels.y) == 0){
    labels.y <- unique(sort(pdfs$y[pdfs$i == 1]))
  }
  ymin.val <- xmin.val
  ymax.val <- xmax.val
  delta.y <- (ymax.val - ymin.val) / length(labels.y)
  y.vals <- log(data.per.cell$y / data.per.cell$nCount_RNA)
  y.vals[is.infinite(y.vals)] <- ymin.val
  data.per.cell$cut_y <- cut(y.vals, 
                             breaks = seq(from =  ymin.val, 
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
pdfs.input$condition <- "postPoisson"

pdfs.pre <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  i <- unique(meta.data$i[meta.data$id == my.id])
  group.name <- my.table$group[my.table$i == i]
  my.group <- unique(meta.data$group[meta.data$id == my.id])
  if (! is.na(group.name) & group.name != ""){
    data.per.cell <- subset(beforeP[[i]], group == gsub(paste0("^", group.name), "", my.group))
  } else {
    data.per.cell <- beforeP[[i]]
  }
  labels.x <- unique(sort(pdfs$x[pdfs$id == my.id]))
  if (length(labels.x) == 0){
    labels.x <- unique(sort(pdfs$x[pdfs$i == 1]))
  }
  
  xmin.val <- my.table[my.table$i == i, "xmin"]
  xmax.val <- my.table[my.table$i == i, "xmax"]
  delta.x <- (xmax.val - xmin.val) / length(labels.x)
  x.vals <- data.per.cell$x
  x.vals[is.infinite(x.vals)] <- xmin.val
  data.per.cell$cut_x <- cut(x.vals, 
                             breaks = seq(from = xmin.val, 
                                          to = xmax.val, 
                                          length.out = length(labels.x) + 1),
                             right = FALSE,
                             labels = labels.x)
  labels.y <- unique(sort(pdfs$y[pdfs$id == my.id]))
  if (length(labels.y) == 0){
    labels.y <- unique(sort(pdfs$y[pdfs$i == 1]))
  }
  ymin.val <- xmin.val
  ymax.val <- xmax.val
  delta.y <- (ymax.val - ymin.val) / length(labels.y)
  y.vals <- data.per.cell$y
  y.vals[is.infinite(y.vals)] <- ymin.val
  data.per.cell$cut_y <- cut(y.vals, 
                             breaks = seq(from = ymin.val, 
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
pdfs.pre$condition <- "prePoisson"

pdfs.sanity.ltq <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  i <- unique(meta.data$i[meta.data$id == my.id])
  my.group <- unique(meta.data$group[meta.data$id == my.id])
  if (my.group != "all"){
    return(NULL)
  }
  data.per.cell <- sanity.res[[i]][, c("mux", "muy")]
  colnames(data.per.cell) <- c("x", "y")
  labels.x <- unique(sort(pdfs$x[pdfs$id == my.id]))
  if (length(labels.x) == 0){
    labels.x <- unique(sort(pdfs$x[pdfs$i == 1]))
  }
  xmin.val <- my.table[my.table$i == i, "xmin"]
  xmax.val <- my.table[my.table$i == i, "xmax"]
  delta.x <- (xmax.val - xmin.val) / length(labels.x)
  x.vals <- data.per.cell$x
  x.vals[is.infinite(x.vals)] <- xmin.val
  data.per.cell$cut_x <- cut(x.vals, 
                             breaks = seq(from = xmin.val, 
                                          to = xmax.val, 
                                          length.out = length(labels.x) + 1),
                             right = FALSE,
                             labels = labels.x)
  labels.y <- unique(sort(pdfs$y[pdfs$id == my.id]))
  if (length(labels.y) == 0){
    labels.y <- unique(sort(pdfs$y[pdfs$i == 1]))
  }
  ymin.val <- xmin.val
  ymax.val <- xmax.val
  delta.y <- (ymax.val - ymin.val) / length(labels.y)
  y.vals <- data.per.cell$y
  y.vals[is.infinite(y.vals)] <- ymin.val
  data.per.cell$cut_y <- cut(y.vals, 
                             breaks = seq(from = ymin.val, 
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
pdfs.sanity.ltq$condition <- "SanityLTQ"

pdfs.sanity.post <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  i <- unique(meta.data$i[meta.data$id == my.id])
  my.group <- unique(meta.data$group[meta.data$id == my.id])
  if (my.group != "all"){
    return(NULL)
  }
  data.per.cell <- sanity.res[[i]]
  df.xy <- unique(pdfs[pdfs$id == my.id, c("x", "y")])
  new.df <- pdf.2d.post(df.xy, data.per.cell)
  new.df$id <- my.id
  return(new.df[, c("x", "y", "value", "id")])
}))
pdfs.sanity.post$condition <- "SanityPost"

pdfs.original.distrib <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  df.xy <- unique(pdfs[pdfs$id == my.id, c("x", "y")])
  if (my.id == "0_0.5_0.5_0_xVS0_0.5_0.5_0_y__5"){
    temp1 <- pdf.2d.norm(df.xy, -9, -7.5, 0.3, 0.12, 0)
    temp2 <- pdf.2d.norm(df.xy, -7.5, -9, 0.12, 0.3, 0)
    new.df <- merge(temp1, temp2, by = c("x", "y"))
    new.df$value <- (new.df$value.x + new.df$value.y) / 2
  } else if (my.id == "1_-9_-9_0.5_0.5_0.5_xVS1_-9_-9_0.5_0.5_0.5_y__1"){
    new.df <- pdf.2d.norm(df.xy, -9, -9, 0.5, 0.5, 0.5)
  } else if (my.id == "1_-9_-9_0.5_0.5_-0.5_xVS1_-9_-9_0.5_0.5_-0.5_y__2"){
    new.df <- pdf.2d.norm(df.xy, -9, -9, 0.5, 0.5, -0.5)
  } else if (my.id == "1_-9_-9_0.5_0.5_0_xVS1_-9_-9_0.5_0.5_0_y__3"){
    new.df <- pdf.2d.norm(df.xy, -9, -9, 0.5, 0.5, 0)
  } else if (my.id == "0.5_0_0_0.5_xVS0.5_0_0_0.5_y__4"){
    temp1 <- pdf.2d.norm(df.xy, -9, -9, 0.3, 0.3, 0)
    temp2 <- pdf.2d.norm(df.xy, -7.5, -7.5, 0.12, 0.12, 0)
    new.df <- merge(temp1, temp2, by = c("x", "y"))
    new.df$value <- (new.df$value.x + new.df$value.y) / 2
  } else if (my.id == "0_0.5_0.5_0_xVS0_0.5_0.5_0_y__7"){
    temp1 <- pdf.2d.norm(df.xy, -9, -8, 0.3, 0.12, 0)
    temp2 <- pdf.2d.norm(df.xy, -8, -9, 0.12, 0.3, 0)
    new.df <- merge(temp1, temp2, by = c("x", "y"))
    new.df$value <- (new.df$value.x + new.df$value.y) / 2
  } else if (my.id == "0.5_0_0_0.5_xVS0.5_0_0_0.5_y__6"){
    temp1 <- pdf.2d.norm(df.xy, -9, -9, 0.3, 0.3, 0)
    temp2 <- pdf.2d.norm(df.xy, -8, -8, 0.12, 0.12, 0)
    new.df <- merge(temp1, temp2, by = c("x", "y"))
    new.df$value <- (new.df$value.x + new.df$value.y) / 2
  } else {
    return(NULL)
  }
  new.df$id <- my.id
  return(new.df[, c("x", "y", "value", "id")])
}))
pdfs.original.distrib$condition <- "OriginalDistrib"

labeller.cond <- c("original\ndistribution",
                   "generated\nexpression", "simulated\nscRNA-seq\nnormalized counts",
                   "expression\ninfered\nby baredSC", "expression\ninfered\nby Sanity",
                   "posterior\n from Sanity")
names(labeller.cond) <- c("OriginalDistrib",
                          "prePoisson", "postPoisson",
                          "baredSC", "SanityLTQ", "SanityPost")

# Get the corr from baredSC
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
corr$condition <- factor(corr$condition, levels = names(labeller.cond))

# Get the cor from Sanity ltq:
corr.sanity.ltq <- do.call(rbind, lapply(unique(pdfs$id), function(my.id){
  i <- unique(meta.data$i[meta.data$id == my.id])
  my.group <- unique(meta.data$group[meta.data$id == my.id])
  if (my.group != "all"){
    return(NULL)
  }
  data.per.cell <- sanity.res[[i]][, c("mux", "muy")]
  return(c(my.id, cor(data.per.cell[, 1], data.per.cell[, 2])))
}))
corr.sanity.ltq <- as.data.frame(corr.sanity.ltq)
colnames(corr.sanity.ltq) <- c("id", "corr")
corr.sanity.ltq$corr <- as.numeric(corr.sanity.ltq$corr)
# Add the meta data
corr.sanity.ltq <- merge(corr.sanity.ltq, meta.data)
corr.sanity.ltq$condition <- "SanityLTQ"
corr.sanity.ltq$condition <- factor(corr.sanity.ltq$condition, levels = names(labeller.cond))


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
prop.model$condition <- factor(prop.model$condition, levels = names(labeller.cond))

# Labeller experiments:
my.labeller <- c("weakly\ncorrelated", "weakly\nanti-correlated", "independent",
                 "2 populations\ncorrelated", "2 populations\nanti-correlated",
                 "2 populations\ncorrelated\ncloser", "2 populations\nanti-correlated\ncloser")
names(my.labeller) <- as.character(1:length(my.labeller))

all.pdfs <- rbind(pdfs[, c("x", "y", "value", "id", "condition")],
                  pdfs.pre, pdfs.input, pdfs.sanity.ltq, pdfs.sanity.post,
                  pdfs.original.distrib)
all.pdfs <- merge(all.pdfs, unique(meta.data[, c("i", "id", "group")]))
all.pdfs$condition <- factor(all.pdfs$condition,
                             levels = names(labeller.cond))

all.pdfs_round <- all.pdfs
all.pdfs_round[, c("x", "y")] <- round(all.pdfs_round[, c("x", "y")], 5)
all.pdfs_round$i <- factor(all.pdfs_round$i, levels = names(my.labeller))

my.cor <- apply(unique(all.pdfs[, c("id", "condition")]), 1, function(v){
  df <- subset(all.pdfs, id == v[1] & condition == v[2])
  cor.input <- corr.from.pdf(df)
  return(c(id = unname(v[1]), corr = corr.from.pdf(df), condition = unname(v[2])))
})
my.cor <- as.data.frame(t(my.cor))
my.cor$corr <- as.numeric(my.cor$corr)
my.cor$condition <- factor(my.cor$condition, levels = names(labeller.cond))
my.cor <- merge(my.cor, unique(meta.data[, c("i", "group", "id")]))

corr$i <- factor(corr$i, levels = names(my.labeller))
my.cor$i <- factor(my.cor$i, levels = names(my.labeller))
prop.model$i <- factor(prop.model$i, levels = names(my.labeller))

ggplot(all.pdfs_round, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  facet_grid(condition ~ i + group) +
  geom_text(data = corr,
            aes(label = label),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 0,
            parse = T) +
  geom_text(data = corr,
            aes(label = p.label),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.4,
            size = 2.5, hjust = 0) +
  geom_text(data = my.cor,
            aes(label = round(corr, 2)),
            x = max(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 1) +
  geom_text(data = corr.sanity.ltq,
            aes(label = round(corr, 2)),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 0) +
  theme_classic() +
  xlab(expression(paste("log(", lambda[g1], ")"))) + 
  ylab(expression(paste("log(", lambda[g2], ")"))) + 
  scale_size("proportion\nof each model",
             guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
  scale_fill_gradient(low="white", high="black") 

g <- ggplot(all.pdfs_round, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  facet_grid(condition ~ i + group, 
             labeller = labeller(i = my.labeller,
                                 condition = labeller.cond)) + # ,
                                 # group = label.group)) +
  geom_text(data = prop.model, aes(x = -8 + 0.5 * ngauss, label = ngauss, size = prop),
            y =  max(all.pdfs_round$y) * 1.1,
            show.legend = F) +
  geom_point(data = prop.model, aes(size = prop),
             x = max(all.pdfs_round$x) * 0.9, y = max(all.pdfs_round$y) * 0.9, colour = NA) +
  geom_text(data = corr,
            aes(label = label),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 0,
            parse = T) +
  geom_text(data = corr,
            aes(label = p.label),
            x = min(all.pdfs_round$x), y = min(all.pdfs_round$y) * 0.9,
            size = 2.5, hjust = 0) +
  geom_text(data = subset(my.cor, ! condition %in% c("baredSC", "SanityLTQ")),
            aes(label = round(corr, 2)),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 0) +
  geom_text(data = corr.sanity.ltq,
            aes(label = round(corr, 2)),
            x = min(all.pdfs_round$x), y = max(all.pdfs_round$y) * 1.1,
            size = 2.5, hjust = 0) +
  theme_classic() +
  xlab(expression(paste("log(", lambda[g1], ")"))) + 
  ylab(expression(paste("log(", lambda[g2], ")"))) + 
  scale_size("proportion\nof each model",
             guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
  scale_fill_gradient(low="white", high="black") 
ggsave(paste0(output.prefix, "_combinedModels_all.pdf"), g, width = 3 + 1.5 * length(unique(all.pdfs_round$id)), height = 3 + 1.5 * length(labeller.cond), limitsize = F)

# Selected experiments to plot:
figs <- list('figS6' = 1:7
)
ggplots <- list()
for(fig.name in names(figs)){
  # Main:
  my.df <- subset(all.pdfs_round, i %in% figs[[fig.name]] & group == "all" & condition != "prePoisson")
  my.corr <- subset(corr, i %in% figs[[fig.name]] & group == "all")
  my.cor.sanity <- subset(corr.sanity.ltq, i %in% figs[[fig.name]] & group == "all")
  my.cor.from.pdf <- subset(my.cor, i %in% figs[[fig.name]] & group == "all" & condition != "prePoisson")
  my.prop.model <- subset(prop.model, i %in% figs[[fig.name]] & group == "all")
  g <- ggplot(my.df, aes(x, y)) +
    geom_tile(aes(fill = log(1 + value))) +
    facet_grid(condition ~ i, labeller = labeller(i = my.labeller,
                                                  condition = labeller.cond)) +
    geom_text(data = my.cor.sanity,
              aes(label =  round(corr, 2)),
              x = min(my.df$x), y = max(my.df$y) * 1.1,
              size = 2.5, hjust = 0) +
    geom_text(data = my.corr,
              aes(label = label),
              x = min(my.df$x), y = max(my.df$y) * 1.1,
              size = 2.5, hjust = 0,
              parse = T) +
    geom_text(data = my.corr,
              aes(label = p.label),
              x = min(my.df$x), y = min(all.pdfs_round$y) * 0.9,
              size = 2.5, hjust = 0) +
    geom_text(data = subset(my.cor.from.pdf, ! condition %in% c("baredSC", "SanityLTQ")),
              aes(label = round(corr, 2)),
              x = min(my.df$x), y = max(my.df$y) * 1.1,
              size = 2.5, hjust = 0) +
    theme_classic() +
    xlab(expression(paste("log(", lambda[g1], ")"))) + 
    ylab(expression(paste("log(", lambda[g2], ")"))) + 
    scale_fill_gradient(low="white", high="black")
    
    ggplots[[fig.name]] <- g
    
  ggsave(paste0(output.prefix, "_combinedModels_", fig.name, ".pdf"), g,
         width = 3 + 1 * length(unique(my.df$i)), height = 2 + 1 * length(unique(my.df$condition)), limitsize = F)
}
