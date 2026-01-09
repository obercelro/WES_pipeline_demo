options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("maftools", quietly = TRUE))
    BiocManager::install("maftools")
library(maftools)

maf_file <- "results/maf/Tumor_1000G.maf"
maf <- read.maf(maf = maf_file, isTCGA = FALSE)
print(maf)

dir.create("results/plots", showWarnings = FALSE)

png(filename = "results/plots/oncoplot.png", width = 800, height = 600)
oncoplot(maf = maf, top = 10, fontSize = 0.8)
dev.off()

png(filename = "results/plots/summary.png", width = 1000, height = 600)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

if ("TP53" %in% maf@gene.summary$Hugo_Symbol) {
    png(filename = "results/plots/lollipop_TP53.png", width = 800, height = 400)
    lollipopPlot(maf = maf, gene = "TP53", AACol = "Protein_Change", showMutationRate = TRUE)
    dev.off()
}

print("analysis complete. check results/plots/ for images.")
