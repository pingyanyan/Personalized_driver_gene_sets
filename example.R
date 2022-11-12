setwd("/pub6/temp/pyy/GA/identifying driver gene sets")
load("1.data preprocess/data/processed/ppi_exp.Rdata")
source("1.data preprocess/program/adjM.R")
#M <- adjM(ppi_exp, "max")##包含两部分网络图g和邻接矩阵M
#save(M, file = "1.data preprocess/data/processed/adjM.Rdata")
M <- adjM_W(ppi_exp,gene_case,"max")##利用基因表达构建共表达网络，包含两部分网络图g和邻接矩阵M
save(M, file = "1.data preprocess/data/processed/adjM_W.Rdata")


W <- NAM(M$M)##对邻接矩阵进行度标化
save(W, file = "1.data preprocess/data/processed/NAM.Rdata")

setwd("/pub6/temp/pyy/GA/identifying driver gene sets")
load("/pub5/xiaoyun/BigData/TCGA/GBM/result/cnv_M.Rdata")##cnv_M
colnames(cnv_M)<-substr(colnames(cnv_M),9,12)
pos<-which(apply(cnv_M,2,sum)==0)
cnv_M <- cnv_M[, - pos]
cnv_M<-as.matrix(cnv_M)

load("/pub5/xiaoyun/BigData/TCGA/GBM/result/mut_M.Rdata")##mut_M
colnames(mut_M)<-substr(colnames(mut_M),9,12)
load("/pub5/xiaoyun/BigData/TCGA/GBM/result/case_normal_exp.Rdata")##gene_case和gene_normal

library(igraph)
library(GA)
library(GSVA)
library(doParallel)
library(memoise)
library(clusterProfiler)
library(enrichplot)
setwd("/pub6/temp/pyy/GA/identifying driver gene sets/")
source("1.data preprocess/program/GA_genes and FC.R")

load("1.data preprocess/data/processed/hallmark_geneset.Rdata")##hallmark_geneset
file <- "1.data preprocess/data/original/h.all.v7.1.entrez.gmt"
gene_h<- read.gmt(file)
load("1.data preprocess/data/processed/NAM.Rdata") ##W
#load("1.data preprocess/data/processed/SM_0.3.Rdata")
load("1.data preprocess/data/processed/SM_W0.3.Rdata")
load("2.driver genesets/data/Cdrivergenes_by_SMW0.3.Rdata") ##Cdrivergenes

source("2.driver genesets/program/RWR.R")
source("2.driver genesets/program/f.R")
source("2.driver genesets/program/main.R")
source("2.driver genesets/program/GA_genes and FC.R")

#patient_id <-  "0882"

#dname <- "GBM1"
dflag <- "c"
#cnv <- cnv_M[, patient_id]
#mut <- mut_M[, patient_id]
meth=NULL
#genes <- GA_genes(dflag, cnv, mut, meth)
#genes<-intersect(genes,Cdrivergenes)
eflag<-"C"
#FC_M <- FC_exp(eflag, case = gene_case,normal=gene_normal,log=F)
#FC <- FC_M[, patient_id, drop = F]
genesets = hallmark_geneset
r = 0.3
u = 1e-9
popSize = 100
pcrossover = 0.9
pmutation = 0.01
elitism = 5
maxiter = 100
pthr=0.05

#mc <- makeCluster(4, type = "FORK", outfile = "/pub6/temp/pyy/GA/identifying driver gene sets/2.driver genesets/debug.txt")
#registerDoParallel(mc)
Dname <- "GBM2"
dflag<-"cm"

driver_set(Dname,flag = dflag, eflag = eflag, cnv_M = cnv_M, mut_M = mut_M, meth_M = NULL, gene_M = gene_case, gene_normal = gene_normal, SM = SM, gene_h = gene_h, genesets_all = genesets, Cdrivergenes = Cdrivergenes, pthr = pthr, pcrossover = pcrossover, pmutation = pmutation, elitism = elitism)


gaControl("binary" = list(selection = "gabin_tourSelection", crossover = "gabin_uCrossover"))



result <- GA_F(genes, FC, genesets, W, r, u, popsize, pcrossover, pmutation, elitism, maxiter, parallel = F)
file<-paste(patiend_id,".Rdata",sep = "")
save(result, file = file)
