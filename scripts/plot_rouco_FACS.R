options(stringsAsFactors = F)

# Dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
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
output.prefix <- commandArgs(TRUE)[3]

setwd(wd)

# Specify labels for limbtypes
label.group <- c(expression(paste("FL ", italic("Pitx1"^"+/+"))), expression(paste("HL ", italic("Pitx1"^"+/+"))), expression(paste("HL ", italic("Pitx1"^"Pen-/Pen-"))))
names(label.group) <- c("FLWT", "HLWT", "HLPendel")
color.group <- c("black", scales::hue_pal()(2))

# Read the data:
data <- do.call(rbind, lapply(names(label.group), function(limbtype){
  values <- read.csv(file.path(directory, paste0(limbtype, ".csv")))$FL1.AREA
  return(data.frame(limbtype = limbtype, value = values))
}))
g1 <- ggplot(data, aes(x = log10(value))) +
  geom_density(aes(color = limbtype), key_glyph = "path") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 4)) +
  theme_classic() + 
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  ggtitle("Fluorescence level") +
  ylab("Density") +
  xlab("log10(GFP)")
g2 <- ggplot(data, aes(x = log(1 + 1e-2 * value))) +
  geom_density(aes(color = limbtype), key_glyph = "path")  +
  coord_cartesian(ylim = c(0, 1.5)) + 
  theme_classic() +
  scale_color_manual(name = "Genotype",
                     values = color.group,
                     breaks = names(label.group),
                     labels = label.group) +
  ggtitle("Fluorescence level\nwith pseudo count") +
  ylab("Density") +
  xlab("log(1 + 0.01 GFP)")
g.all <- ggarrange(g1, g2, labels = c("A", "B"), legend = F)
ggsave(paste0(output.prefix, "_fig_Rouco_ab.pdf"), g.all,
       width = 6, height = 3)
pdf(paste0(output.prefix, "_gaussian_fit.pdf"))
for (limbtype in names(label.group)){
  h <- hist(log10(data$value[data$limbtype == limbtype]),
            breaks = 100, freq = F, main = limbtype,
            xlab = "log10(Fluorescence)")
  if (limbtype == "FLWT"){
    fit <- nls(density ~ (1 / (sigma1 * sqrt(2 * pi)) * exp(-(mids-mean1)**2/(2 * sigma1**2))),
               data=data.frame(density = h$density, mids = h$mids),
               start=list(mean1=1, sigma1=0.3),
               control = nls.control(maxiter = 200))  
    
  } else if (limbtype == "HLWT"){
    fit <- nls(density ~ (A1 / (sigma1 * sqrt(2 * pi)) * exp(-(mids-mean1)**2/(2 * sigma1**2)) +
                            A2 / (sigma2 * sqrt(2 * pi)) * exp(-(mids-mean2)**2/(2 * sigma2**2)) +
                            (1 - A1 - A2) / (sigma3 * sqrt(2 * pi)) * exp(-(mids-mean3)**2/(2 * sigma3**2))),
               data=data.frame(density = h$density, mids = h$mids),
               start=list(A1=0.3, mean1=1, sigma1=0.3,
                          A2=0.3, mean2=2, sigma2=0.3,
                          mean3=3, sigma3=0.3),
               control = nls.control(maxiter = 200))
    
  } else {
    fit <- nls(density ~ (A1 / (sigma1 * sqrt(2 * pi)) * exp(-(mids-mean1)**2/(2 * sigma1**2)) +
                            A2 / (sigma2 * sqrt(2 * pi)) * exp(-(mids-mean2)**2/(2 * sigma2**2)) +
                            (1 - A1 - A2)  / (sigma3 * sqrt(2 * pi)) * exp(-(mids-mean3)**2/(2 * sigma3**2))),
               data=data.frame(density = h$density, mids = h$mids),
               start=list(A1=0.3, mean1=1.3, sigma1=0.3,
                          A2=0.3, mean2=1.7, sigma2=0.3,
                          mean3=2.5, sigma3=0.3),
               control = nls.control(maxiter = 200)) 
    
  }
  lines(h$mids, predict(fit))
  print(coef(fit))
}
dev.off()
