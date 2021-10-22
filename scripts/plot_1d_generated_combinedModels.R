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

# I assume there are not 2 tables with same basename
my.table$input <- basename(my.table$full.path.input)

data <- lapply(unique(my.table$full.path.input), function(fn){
  df <- read.delim(fn)
})
names(data) <- unique(my.table$input)

# Plot the FigS2:
my.gene <- "gauss_0.5_0.75_0.25_gauss_0.5_2_0.2"
my.title <- "N(0.75, 0.25)\nN(2, 0.2)\n2361 cells"
n.bins <- 80
real <- data[[1]][, paste0(my.gene, "_expression")]
from.counts <- log(1 + 10 ^ 4 * data[[1]][, my.gene] / data[[1]]$nCount_RNA)
breaks <- seq(0, 3, length.out = n.bins)
mids <- (breaks[1:(length(breaks) - 1)] + breaks[2:length(breaks)]) / 2
delta.bin <- breaks[2] - breaks[1]
original.distrib <- (dnorm(mids, mean = 0.75, sd = 0.25) + dnorm(mids, mean = 2, sd = 0.2)) / 2 * nrow(data[[1]]) * delta.bin
new.df <- data.frame(value = c(real, from.counts), condition = rep(c("Simulated", "From data"), each = nrow(data[[1]])))
my.colors <- c("Original distribution" = "blue",
               "Simulated" = "lightblue",
               "From data" = "darkred",
               "Density from data" = "red"
)
ggplot(new.df, aes(x = value, fill = condition, color = condition)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = delta.bin) +
  geom_line(data = data.frame(value = mids, y = original.distrib, condition = "Original distribution"),
            aes(y = y)) +
  geom_density(data = subset(new.df, condition == "From data", select = "value"),
               aes(x = value, color = "Density from data", fill = NULL, y = ..density.. * nrow(data[[1]]) * delta.bin)) +
  theme_classic() + 
  geom_blank(data = data.frame(a = my.colors, b = names(my.colors)),
             aes(fill = b, color = b), inherit.aes = F) +
  scale_fill_manual("",
                    values = my.colors,
                    breaks = names(my.colors)) +
  scale_color_manual("",
                     values = my.colors,
                     breaks = names(my.colors)) +
  xlab(expression(paste("log(1 + ", 10^{4}, lambda[g], ")"))) +
  ylab("Number of cells") +
  ggtitle(my.title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output.prefix, "_figS2.pdf"),
       width = 5.5, height = 3.5)


# Get the size of each group
if (length(unique(my.table$input)) == 1 & length(setdiff(unique(my.table$group), "")) == 1){
  nb.per.group <- table(data[[1]][, setdiff(unique(my.table$group), "")])
  names(nb.per.group) <- paste0(setdiff(unique(my.table$group), ""), names(nb.per.group))
  nb.per.group["all"] <- sum(nb.per.group)
} else {
  nb.per.group <- NULL
}

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
  meta.data$info <- apply(meta.data[, c("file", "gene")], 1, function(v){
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
meta.generation <- meta.generation[order(meta.generation$nb, meta.generation$first.model, meta.generation$min.loc, meta.generation$min.scale), ]
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
norm.values <- cbind(log(1 + 10^4 * subset(data[[1]], select = as.character(unique(pdfs$gene))) / data[[1]]$nCount_RNA),
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
all.expected.values.df <- expected.values.df[rep(1:nrow(expected.values.df), length(nb.per.group)), ]
all.expected.values.df$group <- factor(rep(names(nb.per.group), each = nrow(expected.values.df)), levels = names(nb.per.group))

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

# Specify labels for panels
label.group <- paste0(nb.per.group, " cells")
names(label.group) <- names(nb.per.group)
label.generation <- meta.generation$pretty
names(label.generation) <- meta.generation$gene

temp.df.color <- data.frame(x = rep(1, 3), y = 0.5, 
                            groupCol = c('Original distribution',
                                        #  'Density from data before poisson',
                                         'Density from data',
                                         'baredSC'))
# First plot all
my.df <- pdfs
g <- ggplot(my.df) +
  geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
  facet_grid(group ~ gene, scales='free',
             labeller = labeller(group = label.group,
                                 gene = label.generation)) +
  geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "after"),
               aes(x = value, color = "Density from data"), fill = NA,
               trim = T, key_glyph = "path")  +
  # geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "before"),
  #              aes(x = value, color = "Density from data before poisson"), fill = NA,
  #              trim = T, key_glyph = "path")  +
  geom_line(data = subset(all.expected.values.df, gene %in% my.df$gene),
            aes(x = x, y = val, color = "Original distribution")) +
  geom_line(aes(x = x, y = mean, color = "baredSC")) +
  geom_text(data = prop.model, aes(x = 1 + 0.5 * ngauss, label = ngauss, size = prop),
            y = 2, show.legend = F) +
  geom_point(data = prop.model, aes(size = prop),
             x = 1, y = 1, colour = NA) +
  geom_blank(data = temp.df.color, 
             aes(y = y, color = groupCol, fill = groupCol)) +
  theme_classic() +
  scale_size("proportion\nof each model",
             guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
  scale_fill_manual(name = "",
                    values = c('Original distribution' = NA,
                              #  'Density from data before poisson' = NA,
                               'Density from data' = NA,
                               'baredSC' = "darkgreen"),
                    breaks = temp.df.color$groupCol) +
  scale_color_manual(name = "",
                     values = c('Original distribution' = "blue",
                                # 'Density from data before poisson' = "black",
                                'Density from data' = "red",
                                'baredSC' = "darkgreen"),
                     breaks = temp.df.color$groupCol) +
  xlab(expression(paste("log(1 + ", 10^4, lambda[g], ")"))) + 
  ylab("Density") +
  coord_cartesian(ylim=c(0, 3)) +
  expand_limits(x=0)

ggsave(paste0(output.prefix, "_combinedModels_all.pdf"),
       width = 3 + 2 * length(unique(my.df$gene)), height = 6, limitsize = F)

# Selected experiments to plot:
figs <- list('fig2a' = c("gauss_1_0.5_0.5", "gauss_1_1_0.2", "gauss_1_0.75_0.25", "gauss_1_1.5_0.5"),
             'fig2b' = c("uniform_1_0_1", "uniform_1_0_2"),
             'fig2c' = c("gauss_0.5_0.5_0.15_gauss_0.5_1.5_0.15", "gauss_0.5_0.75_0.25_gauss_0.5_2_0.2",
                         "gauss_0.5_1_0.5_gauss_0.5_2.5_0.15",
                         "gauss_0.3_0.5_0.15_gauss_0.4_1.5_0.15_gauss_0.3_2.5_0.15",
                         "gauss_0.3_0.5_0.2_gauss_0.3_1.25_0.2_gauss_0.4_2_0.2"),
             'fig2d' = c("gauss_0.25_0.75_0.25", "gauss_0.25_1.5_0.25")
)
ggplots <- list()
for(fig.name in names(figs)){
  if (fig.name == "fig2a"){
    ymax = 2.5
  } else {
    ymax = 2
  }
  my.df <- subset(pdfs, pdfs$gene %in% figs[[fig.name]])
  g <- ggplot(my.df) +
    geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
    facet_grid(group ~ gene, scales='free',
               labeller = labeller(group = label.group,
                                   gene = label.generation)) +
    geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "after"),
                 aes(x = value, color = "Density from data"), fill = NA,
                 trim = T, key_glyph = "path")  +
    # geom_density(data = subset(norm.values, gene %in% my.df$gene & poisson == "before"),
    #              aes(x = value, color = "Density from data before poisson"), fill = NA,
    #              trim = T, key_glyph = "path")  +
    geom_line(data = subset(all.expected.values.df, gene %in% my.df$gene),
              aes(x = x, y = val, color = "Original distribution")) +
    geom_line(aes(x = x, y = mean, color = "baredSC")) +
    # geom_text(data = subset(prop.model, gene %in% my.df$gene), aes(x = 1 + 0.5 * ngauss, label = ngauss, size = prop),
    #           y = 0.75 * ymax, show.legend = F) +
    # geom_point(data = subset(prop.model, gene %in% my.df$gene), aes(size = prop),
    #            x = 1, y = 0.75 * ymax, colour = NA) +
    geom_blank(data = temp.df.color, 
               aes(y = y, color = groupCol, fill = groupCol)) +
    theme_classic() +
    # scale_size("proportion\nof each model",
    #            guide = guide_legend(override.aes = list(colour = "black", shape = utf8ToInt("1")))) +
    scale_fill_manual(name = "",
                      values = c('Original distribution' = NA,
                                #  'Density from data before poisson' = NA,
                                 'Density from data' = NA,
                                 'baredSC' = "darkgreen"),
                      breaks = temp.df.color$groupCol) +
    scale_color_manual(name = "",
                       values = c('Original distribution' = "blue",
                                  # 'Density from data before poisson' = "black",
                                  'Density from data' = "red",
                                  'baredSC' = "darkgreen"),
                       breaks = temp.df.color$groupCol) +
    xlab(expression(paste("log(1 + ", 10^4, lambda[g], ")"))) + 
    ylab("Density") +
    coord_cartesian(ylim=c(0, ymax)) +
    expand_limits(x=0)
  if (fig.name != 'fig2a'){
    ggplots[[fig.name]] <- g +
      theme(legend.position = "none")
  } else {
    ggplots[[fig.name]] <- g
  }
  ggsave(paste0(output.prefix, "_combinedModels_", fig.name, ".pdf"), g,
         width = 3 + 2 * length(unique(my.df$gene)), height = 6, limitsize = F)
}
leg <- get_legend(ggplots[['fig2a']])
g.all <- 
  ggarrange(ggarrange(ggplots[['fig2a']] + theme(legend.position = "none"),
                      ggplots[['fig2b']],
                      labels = c("A", "B"), widths = c(1.2, 1)),
            ggplots[['fig2c']],
            ggarrange(ggplots[['fig2d']], leg,
                      labels = c("D", ""), widths = c(1, 1.2)),
            ncol = 1, labels = c("", "C", ""))
ggsave(paste0(output.prefix, "_combinedModels_fig2abcd.pdf"), g.all,
       width = 10, height = 15)
