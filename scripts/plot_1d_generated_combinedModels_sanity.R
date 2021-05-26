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
# The sanity results are in the directory
directory.sanity <- commandArgs(TRUE)[5]

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "gene", "xmin", "xmax", "group"))
my.table$i <- rownames(my.table)

# I assume the input basename is unique
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$full.path.input), function(fn){
  df <- read.delim(fn)
})
names(data) <- unique(my.table$input)

# Get the size of each group
if (length(unique(my.table$input)) == 1 & length(setdiff(na.omit(unique(my.table$group)), "")) == 1){
  nb.per.group <- table(data[[1]][, setdiff(unique(my.table$group), "")])
  names(nb.per.group) <- paste0(setdiff(unique(my.table$group), ""), names(nb.per.group))
  nb.per.group["all"] <- sum(nb.per.group)
} else {
  nb.per.group <- NULL
}

# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pretty.*pdf.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pretty_pdf.txt", "", pdf.files), function(v){
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
  meta.data$info <- apply(meta.data[, c("file", "gene")], 1, function(v){
    return(paste(strsplit(gsub("_pdf.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)
}
# Find the group name
meta.data$group <- apply(meta.data[, c("id", "i")], 1, function(v){
  group.name <- my.table$group[my.table$i == v[2]]
  if (! is.na(group.name) & group.name != ""){
    return(paste0(group.name, gsub(paste0("_", v[2], "_pretty$"), "", strsplit(v[1], paste0("_", group.name))[[1]][2])))
  } else {
    return("all")
  }
})
if (! is.null(nb.per.group)){
  meta.data$group <- factor(meta.data$group, levels = names(nb.per.group))
}

# This part is specific to data generated:
# Process in detail the generation
meta.generation <- data.frame(gene = unique(meta.data$gene))
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
meta.generation <- meta.generation[order(meta.generation$nb, meta.generation$first.model, meta.generation$min.loc, meta.generation$min.scale, meta.generation$prop.zero), ]
meta.data$gene <- factor(meta.data$gene, levels = meta.generation$gene)

# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Process the input data:
# Get the values generated before and after poisson
# This is because it is generated data:
colnames(data[[1]]) <- gsub("_\\.", "_-", colnames(data[[1]]))
norm.values <- cbind(log(subset(data[[1]], select = as.character(unique(pdfs$gene))) / data[[1]]$nCount_RNA),
                     subset(data[[1]], select = paste0(as.character(unique(pdfs$gene)), "_expression")),
                     subset(data[[1]], select = "group"))
# Reshape
norm.values <- melt(norm.values, 
                    measure.vars = c(as.character(unique(pdfs$gene)),
                                     paste0(as.character(unique(pdfs$gene)), "_expression")),
                    variable_name = "gene")
# Add a column to specify before or after poisson
norm.values$poisson <- "after"
norm.values$poisson[grepl("expression", norm.values$gene)] <- "before"
norm.values$gene <- gsub("_expression", "", norm.values$gene)
# For each gene substitute -Inf by xmin:
for(my.gene in unique(as.character(norm.values$gene))){
  norm.values$value[norm.values$gene == my.gene & is.infinite(norm.values$value)] <- min(my.table$xmin[my.table$gene == my.gene])
}

# Update the group name
norm.values$group <- paste0("group", norm.values$group)
# Add the values when all cells were used
comb.norm.values <- norm.values
comb.norm.values$group <- "all"
norm.values <- rbind(norm.values, comb.norm.values)

# Put the levels in good order:
if (! is.null(nb.per.group)){
  norm.values$group <- factor(norm.values$group, levels = names(nb.per.group))
}
norm.values$gene <- factor(norm.values$gene, levels = levels(pdfs$gene))

# Plot the distribution used to generate the before poisson
expected.values.df <- do.call(rbind, lapply(unique(pdfs$gene), function(col){
  df <- data.frame(x = unique(pdfs$x[pdfs$gene == col]))
  df$delta.x <- unique(pdfs$delta.x[pdfs$gene == col])
  df$val <- expected.values(as.character(col), df$x)
  df$gene <- col
  return(df)
}))

# Duplicate expected values for each group:
if (! is.null(nb.per.group)){
  all.expected.values.df <- expected.values.df[rep(1:nrow(expected.values.df), length(nb.per.group)), ]
  all.expected.values.df$group <- factor(rep(names(nb.per.group), each = nrow(expected.values.df)), levels = names(nb.per.group))
} else {
  all.expected.values.df <- expected.values.df
  all.expected.values.df$group <- "all"
}
# Get the proportion of each model
prop.model <- do.call(rbind, lapply(unique(meta.data$id), function(my.id){
  temp.df <- as.data.frame(matrix(na.omit(
    as.numeric(unlist(
      strsplit(grep("Using", readLines(file.path(directory, gsub("_pdf.txt", ".log", meta.data$file[meta.data$id == my.id]))),
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

prop.model <- merge(prop.model, unique(meta.data[, c("id", "group", "gene")]))

# Get Sanity results
m.sanity <- as.matrix(read.delim(file.path(directory.sanity, "log_transcription_quotients.txt"), row.names = 1))
m.sanity.e <- as.matrix(read.delim(file.path(directory.sanity, "ltq_error_bars.txt"), row.names = 1))
# Reshape
sanity.values <- melt(t(m.sanity[1:(nrow(m.sanity) - 1), ]))
colnames(sanity.values)[1:2] <- c("cell.id", "gene")
sanity.values$group <- "all"
sanity.values$gene <- factor(sanity.values$gene, levels = meta.generation$gene)

pdf.sanity <- do.call(rbind, lapply(intersect(as.character(sanity.values$gene), as.character(pdfs$gene)), function(my.gene){
  x <- sort(pdfs$x[pdfs$gene == my.gene & pdfs$group == "all"])
  post <- get.post(x,
                   m.sanity[my.gene, ],
                   m.sanity.e[my.gene, ])
  return(data.frame(x = x, mean = post, group = "all", gene = my.gene))
}))
pdf.sanity$gene <- factor(pdf.sanity$gene, levels = meta.generation$gene)
# Specify labels for panels
# label.group <- paste0(nb.per.group, " cells")
# names(label.group) <- names(nb.per.group)
label.generation <- meta.generation$pretty
names(label.generation) <- meta.generation$gene

temp.df.color <- data.frame(x = rep(1, 5), y = 0.5, 
                            groupCol = c('Original distribution',
                                         #  'Density from data before poisson',
                                         'Density from data',
                                         'baredSC',
                                         'Density from Sanity',
                                         'Posterior distribution\nfrom Sanity'))
ggplot.list <- list()
for (my.gene in intersect(unique(my.table$gene), pdfs$gene)){
  my.df <- subset(pdfs, group == "all" & gene == my.gene)
  ggplot.list[[my.gene]] <- ggplot(my.df) +
    geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
    facet_wrap(. ~ gene, scales='free',
               labeller = labeller(gene = label.generation)) +
    geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "after" & group %in% my.df$group),
                 aes(x = value, color = "Density from data"), fill = NA,
                 trim = T, key_glyph = "path")  +
    geom_density(data = subset(sanity.values, gene %in% my.df$gene & group %in% my.df$group),
                 aes(x = value, color = "Density from Sanity"), fill = NA,
                 trim = T, key_glyph = "path")  +
    # geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "before" & group %in% my.df$group),
    #              aes(x = value, color = "Density from data before poisson"), fill = NA,
    #              trim = T, key_glyph = "path")  +
    geom_line(data = subset(pdf.sanity, gene %in% my.df$gene & group %in% my.df$group),
              aes(x = x, y = mean, color = "Posterior distribution\nfrom Sanity")) +
    geom_line(data = subset(all.expected.values.df, gene %in% my.df$gene & group %in% my.df$group),
              aes(x = x, y = val, color = "Original distribution")) +
    geom_line(aes(x = x, y = mean, color = "baredSC")) +
    # geom_text(data = subset(prop.model, group %in% my.df$group),
    #           aes(x = -9 + 0.8 * ngauss, label = ngauss, size = prop),
    #           y = Inf, vjust = 1) +
    geom_blank(data = temp.df.color, 
               aes(y = y, color = groupCol, fill = groupCol)) +
    theme_classic() +
    scale_fill_manual(name = "",
                      values = c('Original distribution' = NA,
                                 #  'Density from data before poisson' = NA,
                                 'Density from data' = NA,
                                 'baredSC' = "darkgreen",
                                 'Density from Sanity' = NA,
                                 'Posterior distribution\nfrom Sanity' = NA),
                      breaks = temp.df.color$groupCol) +
    scale_color_manual(name = "",
                       values = c('Original distribution' = "blue",
                                  # 'Density from data before poisson' = "black",
                                  'Density from data' = "red",
                                  'baredSC' = "darkgreen",
                                  'Density from Sanity' = "purple",
                                  'Posterior distribution\nfrom Sanity' = "plum2"),
                       breaks = temp.df.color$groupCol) +
    xlab(expression(paste("log(", lambda[g], ")"))) + 
    ylab("Density")
  if (meta.generation$prop.zero[meta.generation$gene == my.gene] > 0){
    ggplot.list[[my.gene]] <- ggplot.list[[my.gene]] +
      scale_y_continuous(limits= c(0, 3))
  }
  leg <- get_legend(ggplot.list[[my.gene]])
  ggplot.list[[my.gene]] <- ggplot.list[[my.gene]] +
    theme(legend.position = "none")
}
ggplot.list[["leg"]] <- leg
g <- ggarrange(plotlist = ggplot.list, labels = LETTERS[1:length(unique(my.table$gene))])
ggsave(paste0(output.prefix, "_combinedModels_all_cells.pdf"),
       width = 5 + 3 * length(unique(my.df$gene)), height = 5, limitsize = F)

# # First plot all
# my.df <- pdfs
# g <- ggplot(my.df) +
#   geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
#   facet_grid(group ~ gene, scales='free',
#              labeller = labeller(# group = label.group,
#                                  gene = label.generation)) +
#   geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "after"),
#                aes(x = value, color = "Density from data"), fill = NA,
#                trim = T, key_glyph = "path")  +
#   geom_density(data = subset(sanity.values, gene %in% my.df$gene),
#                aes(x = value, color = "Density from Sanity"), fill = NA,
#                trim = T, key_glyph = "path")  +
#   # geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "before"),
#   #              aes(x = value, color = "Density from data before poisson"), fill = NA,
#   #              trim = T, key_glyph = "path")  +
#   geom_line(data = subset(pdf.sanity, gene %in% my.df$gene),
#             aes(x = x, y = mean, color = "Posterior distribution from Sanity")) +
#   geom_line(data = subset(all.expected.values.df, gene %in% my.df$gene),
#             aes(x = x, y = val, color = "Original distribution")) +
#   geom_line(aes(x = x, y = mean, color = "baredSC")) +
#   geom_text(data = prop.model, aes(x = -9 + 0.8 * ngauss, label = ngauss, size = prop),
#             y = 2, show.legend = F) +
#   geom_point(data = prop.model, aes(size = prop),
#              x = 1, y = 1, colour = NA) +
#   geom_blank(data = temp.df.color, 
#              aes(y = y, color = groupCol, fill = groupCol)) +
#   theme_classic() +
#   scale_size("proportion\nof each model",
#              guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
#   scale_fill_manual(name = "",
#                     values = c('Original distribution' = NA,
#                                #  'Density from data before poisson' = NA,
#                                'Density from data' = NA,
#                                'baredSC' = "darkgreen",
#                                'Density from Sanity' = NA,
#                                'Posterior distribution from Sanity' = NA),
#                     breaks = temp.df.color$groupCol) +
#   scale_color_manual(name = "",
#                      values = c('Original distribution' = "blue",
#                                 # 'Density from data before poisson' = "black",
#                                 'Density from data' = "red",
#                                 'baredSC' = "darkgreen",
#                                 'Density from Sanity' = "purple",
#                                 'Posterior distribution from Sanity' = "plum2"),
#                      breaks = temp.df.color$groupCol) +
#   xlab("log norm expression") + 
#   ylab("Density") +
#   coord_cartesian(ylim=c(0, 3))
# 
# ggsave(paste0(output.prefix, "_combinedModels_all.pdf"),
#        width = 3 + 2 * length(unique(my.df$gene)), height = 6, limitsize = F)
