####################################################
#                                                  #
#   Script for Differential Expression Analyses    #
#   Written by Lucas R Moreira (Sep. 2022)         #
#                                                  #
####################################################

library(tidyverse)
library(limma)
library(edgeR)
library(biomaRt)
library(Glimma)
library(gplots)
library(dplyr)

#########################
#                       #
#   Prepare the data    #
#                       #
#########################

# Let's keep only the longest isoform for each gene
# files <- c("cd01-32.count.tsv","cd01-41.count.tsv","cd01-8mM.count.tsv","cd02-32.count.tsv","cd02-41.count.tsv","cd02-8mM.count.tsv","cd03-32.count.tsv","cd03-41.count.tsv","cd03-8mM.count.tsv","cd04-32.count.tsv","cd04-41.count.tsv","cd04-8mM.count.tsv","cd05-32.count.tsv","cd05-41.count.tsv","cd05-8mM.count.tsv")
# for(i in files){
#   print(i)
#   # read table
#   count <- read.delim(i,header=T)
#   # split ID names
#   gene_names <- count$Gene_ID
#   split <- str_split(gene_names,"[.]")
#   select_two_first <- function(x){x[c(1,2)]}
#   dd  <-  as.data.frame(do.call(rbind, lapply(split,select_two_first)))
#   count$Gene_ID <- dd$V1
#   count$Gene_Name <- dd$V2
#   count$Gene_Length <- count$End-count$Start
#   # keep longest isoform
#   count.ordered <- count[order(count$Gene_Name, -abs(count$Gene_Length) ), ] #sort by id and reverse of abs(value)
#   count.uniq <- count.ordered[ !duplicated(count.ordered$Gene_Name), ]              # take the first row within each id
#   # export table
#   name=str_split(i,"[.]")[[1]][1]
#   write.csv(count.uniq,paste0(name,".longestISO.count.csv"),row.names = F)
# }

files_longestISO <- c("cd01-32.longestISO.count.csv","cd01-41.longestISO.count.csv","cd01-8mM.longestISO.count.csv",
           "cd02-32.longestISO.count.csv","cd02-41.longestISO.count.csv","cd02-8mM.longestISO.count.csv",
           "cd03-32.longestISO.count.csv","cd03-41.longestISO.count.csv","cd03-8mM.longestISO.count.csv",
           "cd04-32.longestISO.count.csv","cd04-41.longestISO.count.csv","cd04-8mM.longestISO.count.csv",
           "cd05-32.longestISO.count.csv","cd05-41.longestISO.count.csv","cd05-8mM.longestISO.count.csv")
read.csv(files_longestISO[1],nrow=5,header=T)

x <- readDGE(files_longestISO, columns=c(2,7),sep = ",")
class(x)

samplenames <- unlist(strsplit(colnames(x), "[.]"))[seq(1,length(unlist(strsplit(colnames(x), "[.]"))),3)]
samplenames

# Group info
treatment <- as.factor(rep(c("32","41","37"),5))
species <- as.factor(rep(c("camel"),15))
individuals <- as.factor(c(rep("cd01",3),rep("cd02",3),rep("cd03",3),rep("cd04",3),rep("cd05",3)))
x$samples$group <- treatment
x$samples$species <- species                
x$samples$individuals <- individuals

x$samples

######################################################################
#                                                                    #
# Add annotation info (easy way when reference genome is on Ensembl) #
#                                                                    #
######################################################################

gene_info <- read.csv(files_longestISO[1],header=T)
annotation <- data.frame(TranscriptID=gene_info$Gene_ID,Gene_description=NA,Chr=gene_info$Reference,
                         Strand=gene_info$Strand,Start=gene_info$Start,End=gene_info$End,Length=gene_info$Gene_Length,GeneName=gene_info$Gene_Name)
head(annotation)

# Let's obtain the transcript ID
y <- readDGE(files_longestISO, columns=c(1,7),sep = ",")
transcriptid <- rownames(y)

#datasets <- listDatasets(ensembl)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#att <- listAttributes(ensembl)
phenotype_description <- select(ensembl, keys=transcriptid, columns=c("ensembl_transcript_id",
                                                                      "ensembl_gene_id","description",
                                                                      "phenotype_description"), 
                keytype="ensembl_transcript_id")
Gene_description <- numeric()
for(transcript in annotation$TranscriptID){
  pd <- phenotype_description[phenotype_description$ensembl_transcript_id==transcript,]
  xy <- data.frame(Gene_description=pd$description[1],phenotype=paste(pd$phenotype_description, collapse = "; "))
  Gene_description <- rbind(Gene_description,xy)
}
annotation$Gene_description <- Gene_description$Gene_description
annotation$Phenotype <- Gene_description$phenotype

go_terms <- select(ensembl, keys=transcriptid, columns=c("ensembl_transcript_id","ensembl_gene_id",
                                                         "external_gene_name","go_id","name_1006"), 
                   keytype="ensembl_transcript_id")
go.terms <- numeric()
for(transcript in annotation$TranscriptID){
  go.t <- go_terms[go_terms$ensembl_transcript_id==transcript,]
  godef <-  paste(go.t$name_1006, collapse = "; ")
  xx <- data.frame(GO_definition=godef)
  go.terms <- rbind(go.terms,xx)
}
annotation <- cbind(annotation,go.terms)

x$genes <- annotation

# Alternatively, the DGE object can be biuld from scratch
#z <- DGEList(counts = Counts, genes = annotation)

#########################
#                       #
#  Scale normalization  #
#                       #
#########################

# Normalization is required to ensure that the expression distributions of each sample are similar across the entire experiment.
# Normalization by the method of trimmed mean of M-values (TMM) (Robinson and Oshlack 2010)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#########################
#                       #
#     Pre-processing    #
#                       #
#########################

# counts per million (CPM) and log2-counts per million (log-CPM) transformations are used regularly although they do not account for gene length differences as RPKM and FPKM values do. 
# Assuming that there are no differences in isoform usage between conditions, differential expression analyses look at gene expression changes between conditions rather than comparing expression across multiple genes or drawing conclusions on absolute levels of expression. In other words, gene lengths remain constant for comparisons of interest and any observed differences are a result of changes in condition rather than changes in gene length.

cpm <- cpm(x)  # counts per million
lcpm <- cpm(x, log=TRUE) # log2-counts per million
#rpkm <- rpkm(x, gene.length = ) # reads per kilobase of transcript per million

#########################
#                       #
#      Filtering        #
#                       #
#########################

hist(rowSums(x$counts==0))

# By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size.keep.exprs <- filterByExpr(x, group=x$samples$group)
keep.exprs <- filterByExpr(x, group=x$samples$group)
x <- x[keep.exprs, keep.lib.sizes=FALSE]
dim(x)

# Keep genes with total counts more than 50.
# A <- rowSums(x$counts)
# isexpr <- A > 50
# 
# x <- x[isexpr,]

# Plot filtering outcome 

# L <- mean(x$samples$lib.size) * 1e-6
# M <- median(x$samples$lib.size) * 1e-6
# c(L, M)
# 
# summary(lcpm)
lcpm.cutoff <- log2(10/(median(x$samples$lib.size) * 1e-6) + 2/(mean(x$samples$lib.size) * 1e-6))
library(RColorBrewer)

nsamples <- ncol(x)
col <- colorRampPalette(brewer.pal(8, "Set2"))(nsamples)

par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.50), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.50), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

########################################
#                                      #
#  Unsupervised clustering of samples  #
#                                      #
########################################

lcpm <- cpm(x, log=TRUE)

pdf("PCA.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
plotMDS(lcpm,labels=x$samples$group, col=c(rep(c("blue"),3),rep(c("red"),3),rep(c("green"),3),rep(c("black"),3),rep(c("yellow"),3)))
title(main="A. Sample groups")
plotMDS(lcpm,labels=x$samples$group, col=rep(c("blue","red","green"),5),dim=c(3,4))
title(main="B. Treatment groups")
dev.off()

glMDSPlot(lcpm, labels=paste(x$samples$group, individuals, sep="_"), 
          groups=x$samples[,c(2,6)])

#########################
#                       #
#    Linear modeling    #
#                       #
#########################

# Create design matrix

design <- model.matrix(~0+treatment+individuals)
design

contr.matrix <- makeContrasts(
  LvN = treatment32 - treatment37, 
  NvH = treatment41 - treatment37, 
  LvH = treatment41 - treatment32, 
  levels = colnames(design))
contr.matrix


##############################
#                            #
#  Differential expression   #
#                            #
##############################

# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself.
# Typically, the voom-plot shows a decreasing trend between the means and variances resulting from a combination of technical variation in the sequencing experiment and biological variation amongst the replicate samples from different cell populations. 
# Experiments with high biological variation usually result in flatter trends, where variance values plateau at high expression values. Experiments with low biological variation tend to result in sharp decreasing trends.
# Moreover, the voom-plot provides a visual check on the level of filtering performed upstream. If filtering of lowly-expressed genes is insufficient, a drop in variance levels can be observed at the low end of the expression scale due to very small counts. If this is observed, one should return to the earlier filtering step and increase the expression threshold applied to the dataset.

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

dt1 <- decideTests(efit)
summary(dt1)

# For a stricter definition on significance, one may require log-fold-changes (log-FCs) to be above a minimum value. 
# The treat method (McCarthy and Smyth 2009) can be used to calculate p-values from empirical Bayes moderated t-statistics with a minimum log-FC requirement.

# tfit <- treat(vfit, lfc=1)
# dt2 <- decideTests(tfit)
# summary(dt2)

par(mfrow=c(1,1))
vennDiagram(dt1[,1:2], circle.col=c("turquoise", "salmon"))

low.vs.normal <- topTable(efit, coef=1, n=Inf, p.value=0.05)
normal.vs.high <- topTable(efit, coef=2, n=Inf, p.value=0.05)
head(low.vs.normal)
head(normal.vs.high)

p.value_all.genes <- efit$p.value

# all genes and their respective p-values
write.table(p.value_all.genes,"p.value_all.genes.txt",quote=F,sep = "\t",row.names = T)

# top 100 genes for N vs H
write.table(head(normal.vs.high,100),"top_100-NvsH.tsv",quote=F,sep = "\t",row.names = T)

# top 100 gene expression
top_100 <- cpm[rownames(cpm) %in% rownames(head(normal.vs.high,100)), ]
write.table(top_100,"top_100-NvsH_all-samples.tsv",quote=F,sep = "\t",row.names = T)

# all DE genes for N vs H
write.table(normal.vs.high,"camel-NvsH_all.tsv",quote=F,sep = "\t",row.names = F)

# all upregulated genes for N vs H
upregulated <- normal.vs.high[normal.vs.high$logFC>0,]
write.table(upregulated,"camel-NvsH_upregulated.tsv",quote=F,sep = "\t",row.names = F)

# all downregulated genes for N vs H
downregulated <- normal.vs.high[normal.vs.high$logFC<0,]
write.table(downregulated,"camel-NvsH_downregulated.tsv",quote=F,sep = "\t",row.names = F)

###############################
#                             #
#  Graphical representation   #
#                             #
###############################

pdf("Vulcano_plot-LvN.pdf",width = 6,height = 5)
plotMD(efit, column=1, status=dt1[,1], main=colnames(efit)[1])
dev.off()

pdf("Vulcano_plot-NvH.pdf",width = 6,height = 5)
plotMD(efit, column=2, status=dt1[,2], main=colnames(efit)[2])
dev.off()

glMDPlot(efit, coef=1, status=dt1, main=colnames(efit)[1],
         side.main="external_gene_name", counts=lcpm, groups=x$samples$group, launch=TRUE,
         html = "LvsN")

glMDPlot(efit, coef=2, status=dt1, main=colnames(efit)[2],
         side.main="external_gene_name", counts=lcpm, groups=x$samples$group, launch=TRUE,
         html = "NvsH")

glMDPlot(efit, coef=3, status=dt1, main=colnames(efit)[3],
         side.main="external_gene_name", counts=lcpm, groups=x$samples$group, launch=TRUE,
         html = "LvsH")

normal.vs.high.topgenes <- normal.vs.high$TranscriptID[1:30]
i <- which(v$genes$TranscriptID %in% normal.vs.high.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
pdf("Heatmap-NvH.pdf",width = 6,height = 5)
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$GeneName, labCol=x$samples$group, 
          col=mycol, trace="none", density.info="none",Rowv = T,dendrogram="column")
dev.off()

# ###############################
# #                             #
# #      GO-term enrichment     #
# #                             #
# ###############################
# 
# write.csv(rownames(normal.vs.high),"normal.vs.high.csv")
# write.csv(rownames(low.vs.normal),"low.vs.normal.csv")
# 
# ###############################
# #                             #
# #       Peak Annotation       #
# #                             #
# ###############################
# 
# atac_annot <- read.table("camel_atac_annot.txt",sep = "\t",header=T)
# 
# table_p <- as.data.frame(p.value_all.genes)
# table_p$Entrez.ID <- rownames(p.value_all.genes)
# 
# atac_annot_express <- merge(atac_annot,table_p, by="Entrez.ID", all=T)
# 
# colnames(atac_annot_express)
# 
# atac_annotation <- data.frame(Peak_ID=atac_annot_express$PeakID..cmd.annotatePeaks.pl.merged.peaks.saf..home.unix.lmoreira.vgb.reference_genomes.Camelus_dromedarius.CamDro2.dna_sm.toplevel.fa..gtf..home.unix.lmoreira.vgb.reference_genomes.Camelus_dromedarius.CamDro2.105.gtf..annStats.annStats.txt..gene.p.value_all.genes.txt.,
#                               Chr=atac_annot_express$Chr,Start=atac_annot_express$Start,End=atac_annot_express$End,
#                               Annotation=atac_annot_express$Annotation,Detailed_annotation=atac_annot_express$Detailed.Annotation,
#                               Distance_to_TSS=atac_annot_express$Distance.to.TSS,Nearest_promoter_ID=atac_annot_express$Nearest.PromoterID,
#                               Entrez_ID=atac_annot_express$Entrez.ID,Nearest_Unigene=atac_annot_express$Nearest.Unigene,
#                               Gene_name=atac_annot_express$Gene.Name,Gene_type=atac_annot_express$Gene.Type,LvN=atac_annot_express$LvN.y,NvH=atac_annot_express$NvH.y,LvH=atac_annot_express$LvH.y)
# 
# write.table(atac_annotation,"atac_annotation.tsv",quote=F,sep = "\t",row.names = F)
