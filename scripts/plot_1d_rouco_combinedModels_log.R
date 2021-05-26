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

get.post <- function(x, mu, e){
  m <- apply(data.frame(mu = mu, e = e), 1, function(v){dnorm(x, mean = v[1], sd = v[2])})
  post <- rowMeans(m)
  return(post)
}

wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgene\txmin\txmax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "gene", "xmin", "xmax", "group"))
my.table$i <- rownames(my.table)

# I assume the input is in the current directory
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$input), function(fn){
  df <- read.delim(list.files(pattern = fn, recursive = T, full.names = T)[1])
})
names(data) <- unique(my.table$input)
nb.per.group <- NULL
# Find pdf files
pdf.files <- grep("pretty",
                  list.files(path = directory, pattern = "1-.*pdf.txt"),
                  value = T,
                  invert = T)

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
norm.values <- cbind(log(subset(data[[1]], select = as.character(unique(pdfs$gene))) / data[[1]]$nCount_RNA),
                     subset(data[[1]], select = unique(my.table$group)))
# Reshape
norm.values <- melt(norm.values, 
                    measure.vars = as.character(unique(pdfs$gene)),
                    variable_name = "gene")

# For each gene substitute -Inf by xmin:
for(my.gene in unique(norm.values$gene)){
  norm.values$value[norm.values$gene == my.gene & is.infinite(norm.values$value)] <- min(my.table$xmin[my.table$gene == my.gene])
}
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

# Run and get Sanity:
all.pdf.sanity <- NULL
all.sanity.values <- NULL
for (my.input in unique(my.table$input)){
  my.data <- data[[my.input]]
  for(group.name in unique(my.table$group[my.table$input == my.input])){
    for (group.value in unique(my.data[, group.name])){
      temp.data <- as.data.frame(my.data[my.data[, group.name] == group.value, unique(my.table$gene[my.table$input == my.input])])
      temp.data[, ncol(temp.data) + 1] <- my.data$nCount_RNA[my.data[, group.name] == group.value] - rowSums(temp.data)
      new.data <- cbind(data.frame(GeneID=c(unique(my.table$gene[my.table$input == my.input]), "complement")), t(temp.data))
      temp.dir <- tempdir()
      write.table(new.data, file.path(temp.dir, "input.txt"), quote = F, sep = "\t", row.names = F)
      system(paste0("Sanity/bin/Sanity -f ", temp.dir, "/input.txt -d ", temp.dir))
      m.sanity <- as.matrix(read.delim(file.path(temp.dir, "log_transcription_quotients.txt"), row.names = 1))
      m.sanity.e <- as.matrix(read.delim(file.path(temp.dir, "ltq_error_bars.txt"), row.names = 1))
      if (nrow(m.sanity) == 2){
        sanity.values <- data.frame(cell.id = colnames(m.sanity),
                                    gene = rownames(m.sanity)[1],
                                    value = m.sanity[1, ])
      } else {
        sanity.values <- melt(t(m.sanity[1:(nrow(m.sanity) - 1), ]))
        colnames(sanity.values)[1:2] <- c("cell.id", "gene")
      }
      sanity.values$group <- paste0(group.name, group.value)
      all.sanity.values <- rbind(all.sanity.values, sanity.values)
      pdf.sanity <- do.call(rbind, lapply(intersect(as.character(sanity.values$gene), as.character(pdfs$gene)), function(my.gene){
        x <- sort(pdfs$x[pdfs$gene == my.gene & pdfs$group == paste0(group.name, group.value)])
        post <- get.post(x,
                         m.sanity[my.gene, ],
                         m.sanity.e[my.gene, ])
        return(data.frame(x = x, mean = post, gene = my.gene))
      }))
      pdf.sanity$group <- paste0(group.name, group.value)
      all.pdf.sanity <- rbind(all.pdf.sanity, pdf.sanity)
    }
    # With all:
    if (any(is.na(my.table$group[my.table$input == my.input]))){
      temp.data <- as.data.frame(my.data[, unique(my.table$gene[my.table$input == my.input])])
      temp.data[, ncol(temp.data) + 1] <- my.data$nCount_RNA - rowSums(temp.data)
      new.data <- cbind(data.frame(GeneID=c(unique(my.table$gene[my.table$input == my.input]), "complement")), t(temp.data))
      temp.dir <- tempdir()
      write.table(new.data, file.path(temp.dir, "input.txt"), quote = F, sep = "\t", row.names = F)
      system(paste0("Sanity/bin/Sanity -f ", temp.dir, "/input.txt -d ", temp.dir))
      m.sanity <- as.matrix(read.delim(file.path(temp.dir, "log_transcription_quotients.txt"), row.names = 1))
      m.sanity.e <- as.matrix(read.delim(file.path(temp.dir, "ltq_error_bars.txt"), row.names = 1))
      if (nrow(m.sanity) == 2){
        sanity.values <- data.frame(cell.id = colnames(m.sanity),
                                    gene = rownames(m.sanity)[1],
                                    value = m.sanity[1, ])
      } else {
        sanity.values <- melt(t(m.sanity[1:(nrow(m.sanity) - 1), ]))
        colnames(sanity.values)[1:2] <- c("cell.id", "gene")
      }
      sanity.values$group <- "all"
      all.sanity.values <- rbind(all.sanity.values, sanity.values)
      pdf.sanity <- do.call(rbind, lapply(intersect(as.character(sanity.values$gene), as.character(pdfs$gene)), function(my.gene){
        x <- sort(pdfs$x[pdfs$gene == my.gene & pdfs$group == "all"])
        post <- get.post(x,
                         m.sanity[my.gene, ],
                         m.sanity.e[my.gene, ])
        return(data.frame(x = x, mean = post, gene = my.gene))
      }))
      pdf.sanity$group <- "all"
      all.pdf.sanity <- rbind(all.pdf.sanity, pdf.sanity)
    }
  }
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
  # geom_density(data = subset(all.sanity.values, gene %in% my.df$gene),
  #              aes(x = value, color = "Density from Sanity"), fill = NA,
  #              trim = T, key_glyph = "path")  +
  geom_density(data = subset(full.norm.values, gene %in% my.df$gene & group %in% my.df$group),
               aes(x = value, color = "Density from data"), fill = NA,
               trim = T, key_glyph = "path")  +
  geom_line(aes(x = x, y = mean, color = "baredSC")) +
  geom_line(data = subset(all.pdf.sanity, gene %in% my.df$gene),
            aes(x = x, y = mean, color = "Posterior distribution from Sanity")) +
  theme_classic()  +
  geom_blank(data = data.frame(x = rep(1, 4), y = 0.5, 
                               groupCol = c("Density from data",
                                            "baredSC",
                                            "Density from Sanity",
                                            "Posterior distribution from Sanity")), 
             aes(y = y, color = groupCol, fill = groupCol)) +
  scale_fill_manual(name = "",
                    values = c('Density from data' = NA,
                               'baredSC' = "darkgreen",
                               'Density from Sanity' = NA,
                               "Posterior distribution from Sanity" = NA)) +
  scale_color_manual(name = "",
                     values = c('Density from data' = "red",
                                'baredSC'="darkgreen",
                                'Density from Sanity' = "purple",
                                "Posterior distribution from Sanity" = "pink")) +
  xlab(expression(paste("log(", italic("Pitx1"), ")"))) + 
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1))

ggsave(paste0(output.prefix, "_combinedModels_all.pdf"),
       width = 3 + 2 * length(unique(my.df$gene)), height = 6, limitsize = F)

gd <- ggplot(subset(full.norm.values, gene %in% my.df$gene & group %in% my.df$group)) +
  geom_density(aes(x = value, color = group), fill = NA,
               trim = T, key_glyph = "path") +
  theme_classic() +
  xlab(expression(paste("log(", italic("Pitx1"), ")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  ggtitle("Normalized expression\nfrom scRNA-seq")
gb <- ggplot(my.df) +
  geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = group), alpha = 0.3, color = NA) +
  geom_line(aes(x = x, y = mean, color = group), key_glyph = "path") +
  theme_classic() +
  xlab(expression(paste("log(", italic("Pitx1"), ")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  scale_fill_manual(name = "Genotype",
                    values = color.group,
                    breaks = names(label.group),
                    labels = label.group) +
  ggtitle("Normalized expression\nafter baredSC")

# gsv <- ggplot(all.sanity.values) +
#   geom_density(aes(x = value, color = group), fill = NA,
#                trim = T, key_glyph = "path")  +
#   theme_classic() +
#   xlab("log norm expression") + 
#   ylab("Density") +
#   coord_cartesian(ylim=c(0, 1.5)) +
#   scale_color_manual(name = "Genotype",
#                      values = color.group,
#                      breaks = names(label.group),
#                      labels = label.group) +
#   scale_fill_manual(name = "Genotype",
#                     values = color.group,
#                     breaks = names(label.group),
#                     labels = label.group) +
#   ggtitle("Density of Sanity values")
# 
# 
gs <- ggplot(all.pdf.sanity) +
  geom_line(aes(x = x, y = mean, color = group)) +
  theme_classic() +
  xlab(expression(paste("log(", italic("Pitx1"), ")"))) +
  ylab("Density") +
  coord_cartesian(ylim=c(0, 1.5)) +
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  scale_fill_manual(name = "Genotype",
                    values = color.group,
                    breaks = names(label.group),
                    labels = label.group) +
  ggtitle("Posterior distribution\nfrom Sanity")

# Read the FACS data:
data <- do.call(rbind, lapply(gsub("limbtype", "", names(label.group)), function(limbtype){
  values <- read.csv(file.path(directory, paste0("../", limbtype, ".csv")))$FL1.AREA
  return(data.frame(limbtype = paste0("limbtype", limbtype), value = values))
}))
gf <- ggplot(data, aes(x = log(value))) +
  geom_density(aes(color = limbtype), key_glyph = "path") +
  coord_cartesian(ylim = c(0, 0.8), xlim = c(1, 8.5)) +
  theme_classic() + 
  xlab("log(GFP)") +
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  ggtitle("Fluorescence level")
g.all <- ggarrange(gf, gb, gd, gs, labels = LETTERS[1:4],
                   legend = "bottom",
                   legend.grob = get_legend(gb, "bottom"))
ggsave(paste0(output.prefix, "_combinedModels_figS_Rouco.pdf"), g.all,
       width = 6, height = 6)
