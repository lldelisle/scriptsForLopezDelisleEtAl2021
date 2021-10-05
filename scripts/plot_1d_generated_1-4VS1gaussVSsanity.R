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
pdf.files <- list.files(path = directory, pattern = "1.*gauss.*.*pretty.*pdf.txt")

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

meta.data$post <- gsub("_pdf.txt", "_posterior_per_cell.txt", meta.data$file)

label.generation <- meta.generation$pretty
names(label.generation) <- meta.generation$gene

label.i <- label.generation[sort(match(my.table$gene, meta.generation$gene))]
names(label.i) <- my.table$i[order(match(my.table$gene, meta.generation$gene))]

# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data[, c("file", "i", "post", "group", "gene", "model")])
pdfs$i <- factor(pdfs$i, levels = names(label.i))

# Get the posterior per cell from baredSC:
baredsc.post <- do.call(rbind, lapply(meta.data$post, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$post <- fn
  return(df)
}))
# Add the meta data
baredsc.post <- merge(baredsc.post, meta.data[, c("i", "post", "group", "gene", "model")])
baredsc.post$i <- factor(baredsc.post$i, levels = names(label.i))

pdf.baredsc.post <- do.call(rbind, lapply(meta.data$post, function(fn){
  x <- sort(pdfs$x[pdfs$post == fn & pdfs$group == "all"])
  post <- get.post(x,
                   baredsc.post$mu[baredsc.post$post == fn],
                   baredsc.post$sd[baredsc.post$post == fn])
  return(data.frame(x = x, mean = post, group = "all", post = fn))
}))
# Add the meta data
pdf.baredsc.post <- merge(pdf.baredsc.post, meta.data[, c("i", "post", "group", "gene", "model")])
pdf.baredsc.post$gene <- factor(pdf.baredsc.post$gene, levels = meta.generation$gene)
pdf.baredsc.post$i <- factor(pdf.baredsc.post$i, levels = names(label.i))

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

# Add the i info
norm.values <- merge(norm.values, unique(meta.data[, c("gene", "i")]))

# For each id substitute -Inf by xmin:
for(my.i in unique(as.character(norm.values$i))){
  norm.values$value[norm.values$i == my.i & is.infinite(norm.values$value)] <- my.table$xmin[my.table$i == my.i]
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
norm.values$i <- factor(norm.values$i, levels = names(label.i))

# Plot the distribution used to generate the before poisson
expected.values.df <- do.call(rbind, lapply(unique(pdfs$i), function(my.i){
  df <- data.frame(x = unique(pdfs$x[pdfs$i == my.i]))
  df$delta.x <- unique(pdfs$delta.x[pdfs$i == my.i])
  my.gene <- my.table$gene[my.table$i == my.i]
  df$val <- expected.values(as.character(my.gene), df$x)
  df$gene <- my.gene
  df$i <- my.i
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
all.expected.values.df$i <- factor(all.expected.values.df$i, levels = names(label.i))

# Get Sanity results
m.sanity <- as.matrix(read.delim(file.path(directory.sanity, "log_transcription_quotients.txt"), row.names = 1))
m.sanity.e <- as.matrix(read.delim(file.path(directory.sanity, "ltq_error_bars.txt"), row.names = 1))
# Reshape
sanity.values <- melt(t(m.sanity[1:(nrow(m.sanity) - 1), ]))
colnames(sanity.values)[1:2] <- c("cell.id", "gene")
sanity.values$group <- "all"
sanity.values$gene <- factor(sanity.values$gene, levels = meta.generation$gene)

# Add the i info
sanity.values <- merge(sanity.values, unique(meta.data[, c("gene", "i")]))
sanity.values$i <- factor(sanity.values$i, levels = names(label.i))

pdf.sanity <- do.call(rbind, lapply(unique(sanity.values$i), function(my.i){
  x <- sort(pdfs$x[pdfs$i == my.i & pdfs$group == "all"])
  my.gene <- my.table$gene[my.table$i == my.i]
  post <- get.post(x,
                   m.sanity[my.gene, ],
                   m.sanity.e[my.gene, ])
  return(data.frame(x = x, mean = post, group = "all", gene = my.gene, i = my.i))
}))
pdf.sanity$gene <- factor(pdf.sanity$gene, levels = meta.generation$gene)
pdf.sanity$i <- factor(pdf.sanity$i, levels = names(label.i))


sanity.priors <- data.frame(
  mu = read.delim(file.path(directory.sanity, "mu.txt"), h = F)$V1,
  dmu = read.delim(file.path(directory.sanity, "d_mu.txt"), h = F)$V1,
  var = read.delim(file.path(directory.sanity, "variance.txt"), h = F)$V1,
  gene = read.delim(file.path(directory.sanity, "geneID.txt"), h = F)$V1
)

# Add the id info
sanity.priors<- merge(sanity.priors, unique(meta.data[, c("gene", "i")]))

pdf.sanity.prior <- do.call(rbind, lapply(unique(sanity.priors$i), function(my.i){
  x <- sort(pdfs$x[pdfs$i == my.i & pdfs$group == "all"])
  my.gene <- my.table$gene[my.table$i == my.i]
  mu <- sanity.priors$mu[sanity.priors$i == my.i]
  dmu <- sanity.priors$dmu[sanity.priors$i == my.i]
  my.var <- sanity.priors$var[sanity.priors$i == my.i]
  my.range <- data.frame(dnorm(x, mu - dmu, sqrt(my.var)), dnorm(x, mu + dmu, sqrt(my.var)))
  return(data.frame(x = x, mean = dnorm(x, mu, sqrt(my.var)), low = apply(my.range, 1, min), high = apply(my.range, 1, max), group = "all", gene = my.gene, i = my.i))
}))
pdf.sanity.prior$gene <- factor(pdf.sanity.prior$gene, levels = meta.generation$gene)
pdf.sanity.prior$i <- factor(pdf.sanity.prior$i, levels = names(label.i))

# Plot main figure
my.fill.colors <- c('Original distribution' = "white",
                    #  'Density from data before poisson' = "white",
                    'Density from data' = "white",
                    # 'baredSC_1gauss' = "green",
                    'baredSC' = "darkgreen",
                    'Density from Sanity' = "white",
                    'Posterior distribution\nfrom Sanity' = "white")
my.colors <- c('Original distribution' = "blue",
               # 'Density from data before poisson' = "black",
               'Density from data' = "red",
               # 'baredSC_1gauss' = "green",
               'baredSC' = "darkgreen",
               'Density from Sanity' = "purple",
               'Posterior distribution\nfrom Sanity' = "plum2")
ggplot.list <- list()
for (my.i in sort(setdiff(unique(pdfs$i), "6"))){
  my.df <- subset(pdfs, group == "all" & i == my.i)
  my.gene <- unique(my.df$gene)
  ggplot.list[[my.i]] <- ggplot(subset(my.df, model == "1-4gauss")) +
    geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = "baredSC"), alpha = 0.3, color = NA) +
    facet_wrap(. ~ i, scales='free',
               labeller = labeller(i = label.i)) +
    geom_density(data = subset(norm.values, i %in% my.df$i & poisson == "after" & group %in% my.df$group),
                 aes(x = value, color = "Density from data"), fill = NA,
                 trim = T, key_glyph = "path")  +
    geom_density(data = subset(sanity.values, i %in% my.df$i & group %in% my.df$group),
                 aes(x = value, color = "Density from Sanity"), fill = NA,
                 trim = T, key_glyph = "path")  +
    geom_line(data = subset(pdf.sanity, i %in% my.df$i & group %in% my.df$group),
              aes(x = x, y = mean, color = "Posterior distribution\nfrom Sanity")) +
    geom_line(data = subset(all.expected.values.df, i %in% my.df$i & group %in% my.df$group),
              aes(x = x, y = val, color = "Original distribution")) +
    geom_line(aes(x = x, y = mean, color = "baredSC")) +
    geom_blank(data = data.frame(groupCol = names(my.colors)), 
               aes(color = groupCol, fill = groupCol)) +
    theme_classic() +
    scale_fill_manual(name = "",
                      values = my.fill.colors,
                      breaks = names(my.colors)) +
    scale_color_manual(name = "",
                       values = my.colors,
                       breaks = names(my.colors)) +
    xlab(expression(paste("log(", lambda[g], ")"))) + 
    ylab("Density")
  if (meta.generation$prop.zero[meta.generation$gene == my.gene] > 0){
    ggplot.list[[my.i]] <- ggplot.list[[my.i]] +
      scale_y_continuous(limits= c(0, 3))
  }
  leg <- get_legend(ggplot.list[[my.i]])
  ggplot.list[[my.i]] <- ggplot.list[[my.i]] +
    theme(legend.position = "none")
}
ggplot.list[["leg"]] <- leg
g <- ggarrange(plotlist = ggplot.list, labels = LETTERS[1:(length(ggplot.list) - 1)])
ggsave(paste0(output.prefix, "_combinedModels_all_cells.pdf"),
       width = 5 + 3 * length(unique(my.df$gene)), height = 5, limitsize = F)

# Plot the comparison with Sanity step by step
plotted.i <- setdiff(names(label.i), "2")
my.df <- subset(pdfs, i %in% plotted.i)
my.df$colorGroup <- "baredSC"
my.df$colorGroup[my.df$model == "1gauss"] <- "baredSC_1gauss"
my.df$panel <- "A"
my.norm.values <- subset(norm.values, group == "all" & i %in% plotted.i)
my.norm.values$colorGroup <- "Density from data"
my.norm.values$colorGroup[my.norm.values$poisson == "before"] <- "Simulated expression"
my.norm.values$panel <- "A"
my.norm.values.post <- subset(my.norm.values, poisson == "after")
my.sanity.prior <- subset(pdf.sanity.prior, i %in% plotted.i)
my.sanity.prior$colorGroup <- "Sanity"
my.sanity.prior$panel <- "A"
my.expected.values <- subset(all.expected.values.df, i %in% plotted.i)
my.expected.values$colorGroup <- "Original distribution"
my.expected.values$panel <- "A"
my.expected.values$mean <- my.expected.values$val
my.baredsc.post <- subset(baredsc.post, i %in% plotted.i)
my.baredsc.post$colorGroup <- "baredSC"
my.baredsc.post$colorGroup[my.baredsc.post$model == "1gauss"] <- "baredSC_1gauss"
my.baredsc.post$panel <- "B"
my.baredsc.post$value <- my.baredsc.post$mu
my.sanity.values <- subset(sanity.values, i %in% plotted.i)
my.sanity.values$colorGroup <- "Sanity"
my.sanity.values$panel <- "B"
my.pdf.baredsc.post <- subset(pdf.baredsc.post, i %in% plotted.i)
my.pdf.baredsc.post$colorGroup <- "baredSC"
my.pdf.baredsc.post$colorGroup[my.pdf.baredsc.post$model == "1gauss"] <- "baredSC_1gauss"
my.pdf.baredsc.post$panel <- "C"
my.pdf.sanity <- subset(pdf.sanity, i %in% plotted.i)
my.pdf.sanity$colorGroup <- "Sanity"
my.pdf.sanity$panel <- "C"
col.ribbon <- c("x", "low", "high", "colorGroup", "panel", "i")
my.ribbon <- do.call(rbind, list(my.df[, col.ribbon],
                                 my.sanity.prior[, col.ribbon]
))
col.lines <- c("x", "mean", "colorGroup", "panel", "i")
my.lines <- do.call(rbind, list(my.df[, col.lines],
                                my.sanity.prior[, col.lines],
                                my.expected.values[, col.lines],
                                my.pdf.baredsc.post[, col.lines],
                                my.pdf.sanity[, col.lines]
))
col.density <- c("value", "colorGroup", "panel", "i")
my.density <- do.call(rbind, list(my.norm.values.post[, col.density],
                                  my.baredsc.post[, col.density],
                                  my.sanity.values[, col.density]
                                  
))
my.panel.labeller <- c("Inferred\nexpression\ndistribution", "Density of\nexpected value\nfor each cell", "Posterior\nusing\nerror bars")
names(my.panel.labeller) <- LETTERS[1:3]
my.cond <-  c('Original distribution',
              'Density from data',
              'baredSC',
              'baredSC_1gauss',
              'Sanity')
ggplot() +
  geom_ribbon(data = my.ribbon,
              aes(x = x, ymin = low, ymax = high, fill = colorGroup), alpha = 0.3, color = NA) +
  geom_line(data = my.lines,
            aes(x = x, y = mean, color = colorGroup)) +
  geom_density(data = my.density,
               aes(x = value, color = colorGroup), fill = NA,
               trim = F, key_glyph = "path")  +
  geom_blank(data = data.frame(groupCol = unique(c(my.ribbon$colorGroup, my.lines$colorGroup, my.density$colorGroup))), 
             aes(color = groupCol, fill = groupCol)) +
  facet_grid(panel ~ i, scales='free', switch =  "y",
             labeller = labeller(i = label.i,
                                 panel = my.panel.labeller)) +
  theme_classic() +
  xlab(expression(paste("log(", lambda[g], ")"))) + 
  ylab("Density") +
  scale_fill_manual(name = "",
                    values = c('Original distribution' = "white",
                               'Density from data' = "white",
                               'baredSC' = "darkgreen",
                               'baredSC_1gauss' = 'green',
                               'Sanity' = "purple"),
                    breaks = my.cond) +
  scale_color_manual(name = "",
                     values = c('Original distribution' = "blue",
                                'Density from data' = "red",
                                'baredSC' = "darkgreen",
                                'baredSC_1gauss' = 'green',
                                'Sanity' = "purple"),
                     breaks = my.cond)
ggsave(paste0(output.prefix, "_aligned.pdf"),
       width = 7, height = 4, limitsize = F)
