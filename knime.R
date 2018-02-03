.libPaths("C:/Datos/R-3.4.3/library/")

# Load data
library(idendr0)
df <- read.csv("commercial.csv")
# Compute distances and hierarchical clustering
# Levenshtein Distance
d  <- adist(df$Sequence)
rownames(d) <- df$Header
hc <- hclust(as.dist(d))
idendro(hc, df, hscale = 3, vscale = 2, clusters = 15)



##R SERVER
library(Rserve)
Rserve(args='--vanilla')


plot(hc)