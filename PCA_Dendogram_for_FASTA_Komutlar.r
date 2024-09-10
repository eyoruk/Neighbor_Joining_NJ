                                                ## PCA_for_FASTA ##

> setwd('C:/Users/username/Desktop/PCA_NJ_DENDOGRAM_FASTA/')
> library(ape)
> library(phangorn)
> library(factoextra)
> dir()
[1] "ATP8.fasta"  "Rentrez.rmd"
> atp<-read.FASTA("ATP8.fasta")
> View(atp)
> dm  <- dist.ml(atp)
> pca<-prcomp(dm, scale = TRUE)
> fviz_eig(pca)

> xlsx::write.xlsx(pca$x,"pca.xlsx")
> dir()
[1] "ATP8.fasta"  "pca.xlsx"    "Rentrez.rmd"

> dm  <- dist.ml(atp)
> treeUPGMA  <- upgma(dm)
> plot(treeUPGMA, main="UPGMA")
> treeNJ  <- NJ(dm)
> library(phylocanvas)
> phylocanvas(treeNJ,treetype = "rectangular", alignlabels = T)
> class(atp)
[1] "DNAbin"
> atpx=as.phyDat(atp)
> parsimony(treeUPGMA, atpx)
> treeRatchet  <- pratchet(atpx, trace = 0, minit=100)
> parsimony(treeRatchet, atpx)
> treeRatchet  <- acctran(treeRatchet, atpx)
> if(inherits(treeRatchet, "multiPhylo")){
+     treeRatchet <- unique(treeRatchet)
+ }

> plotBS(midpoint(treeRatchet), type="phylogram")
                Error in plotBS(midpoint(treeRatchet), type = "phylogram") : 
                    You need to supply 'trees' or the tree needs support-values as node.label

> phylocanvas(treeRatchet,treetype = "rectangular", alignlabels = T)
> mt <- modelTest(atpx, model=c("K80", "HKY", "GTR"),
+                 control = pml.control(trace = 0))
> fit <- as.pml(mt, "BIC")

> bs <- bootstrap.pml(fit, bs=5, optNni=TRUE,
+                     control = pml.control(trace = 0))
                Error in edgeMatrix[k, 1:6] <- c(ab, cd, i, f) : 
                    number of items to replace is not a multiple of replacement length
                Error in edgeMatrix[k, 1:6] <- c(ab, cd, i, f) : 
                    number of items to replace is not a multiple of replacement length
