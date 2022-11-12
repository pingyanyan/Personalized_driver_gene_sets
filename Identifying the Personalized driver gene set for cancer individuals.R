##GA_genes identifys the candidate driver genes for a cancer individual
GA_genes <- function(flag, cnv = NULL, mut = NULL, meth = NULL) {
    ##flag: indicate what type of data type for identifying the driver genes for a cancer individual,c for copy number,m for mutation,
    ########me for methylation,cm for copy number and mutation,cme for copy number and methylation,mme for mutation and methylation,
    ########cmme for copy number,muation and methylation
    ##cnv: is a 0-1 vector of copy number alteration for genes,1 represent copy number alteration and 0 not
    ##mut: is a 0-1 vector of mutation for genes,1 represent copy number alteration and 0 not
    ##meth:is a 0-1 vector of methylation for genes,1 represent hyper/hypo-methylation and 0 not
    if (flag == "c") {
        if (is.null(cnv))
            cat("Error:cnv is null")
        genes <- names(cnv[which(cnv == 1)])
    } else if (flag == "m") {
        if (is.null(mut))
            cat("Error:mut is null")
        genes <- names(mut[which(mut == 1)])

    } else if (flag == "me") {
        if (is.null(meth))
            cat("Error:meth is null")
        genes <- names(me[which(me == 1)])

    } else if (flag == "cm") {
        if (is.null(cnv) | is.null(mut))
            cat("Error:cnv or mut is null")
        genes <- unique(c(names(cnv[which(cnv == 1)]), names(mut[which(mut == 1)])))
    } else if (flag == "cme") {
        if (is.null(cnv) | is.null(meth)) {
            cat("Error:cnv or meth is null")
        }
        genes <- unique(c(names(cnv[which(cnv == 1)]), names(me[which(me == 1)])))
    } else if (flag == "mme") {
        if (is.null(meth) | is.null(mut))
            cat("Error:meth or mut is null")
        genes <- unique(c(names(mut[which(mut == 1)]), names(me[which(me == 1)])))

    } else if (flag == "cmmme") {

        if (is.null(cnv) | is.null(mut) | is.null(meth))
            cat("Error:cnv or mut or meth is null")
        genes <- unique(c(names(cnv[which(cnv == 1)]), names(mut[which(mut == 1)]), names(me[which(me == 1)])))

    }

    return(genes)
}
##FC_exp calculates the expression flod change of genes in cancer individuals
FC_exp <- function(eflag, case, normal = NULL,log=F) {

    ##eflag: indicate whether expression profiles of cancer samples or normal sample were used for calculating fold change of genes
    #######T for cancer samples和C for normal samples
    if (eflag == "C") {
        if (log) { 
          case <-2^case / rowMeans(2^normal)
        } else {
          case<-case/rowMeans(normal)

        }
    } else if (eflag == "T") {
        if (log) {
            case <- 2^case / rowMeans(2^case)
        } else {
            case <- case / rowMeans(case)

        }

    }
    return(case)


}
##############################################################################################################
##this part is random walk with restart
##adjM get the max component of ppi network and its adjacent matrix
adjM <- function(ppi, flag) {
    ##ppi:data.frame with two colunms for protein pairs
    ##flag="max" extract the max component of ppi network

    g <- graph_from_data_frame(ppi, directed = F)
    if (flag=="max") {
        gene_list <- groups(components(g))
        l <- which.max(unlist(lapply(gene_list, length)))
        genes <- gene_list[[l]]
        g<-induced_subgraph(g, genes)
    }
    M <- get.adjacency(g)
    M <- as(M, "matrix")
    return(list(g = g, M = M))
}
###adjM_W get the max component of ppi network and its weighted adjacent matrix
adjM_W <- function(ppi,gene_M,flag) {

    ##ppi:data.frame,包含蛋白质互作两列
    ##flag用于标识是否提取网络最大连接组分
      pos<-which(ppi[,1]==ppi[,2])
      if(length(pos)){
        ppi<-ppi[-pos,]
      }
      pos1<-which(ppi[,1]%in%rownames(gene_M))
      pos2<-which(ppi[,2]%in%rownames(gene_M))
      pos<-intersect(pos1,pos2)
      ppi<-ppi[pos,]
    g <- graph_from_data_frame(ppi, directed = F)
    if (flag=="max") {
        gene_list <- groups(components(g))
        l <- which.max(unlist(lapply(gene_list, length)))
        genes <- gene_list[[l]]
        g<-induced_subgraph(g, genes)
    }
    M <- as_adjacency_matrix(g)
    M <- as(M, "matrix")
    M[M!=0]<-1
    cor_M<-abs(cor(t(gene_M[rownames(M),])))
    M<-M*cor_M
   
    return(list(g = g, M = M))
}


##NAM tramsform the adjacent matrix of ppi network into matrix of transition probability through normalizing according to degrees of genes
NAM <- function(M) {
    ##M is the adjacent matrix of ppi network
    W <- sweep(M, 2, colSums(M), "/")##normalize according to degrees of genes
    return(W)
}

###steady_SM calculate the stable Similarity matrix of ppi network when random walk with restart stop
steady_SM <- function(transM, r) {
    ##transM:matrix of transition probability，the colnum sum is 1
    ##r: the restart rate
    I <- diag(1, dim(transM)[1])
    temp <- (I - (1 - r) * transM)
    temp_inv <- solve(temp, I)
    result <- r * temp_inv
    rownames(result) <- rownames(transM)
    colnames(result) <- colnames(transM)
    return(result)
}

##RWR1 is the random walk with restart to calculate the stable trasfer probility of genes in the ppi network driven by the subset of genes
RWR1 <- function(seeds,SM){
    ##seeds is the subset of candidate driver genes needed to be in ppi network
    ##SM:the stable Similarity matrix of ppi network when random walk with restart stop
    p0 <- matrix(0, nrow = dim(SM)[1], ncol = 1)
    rownames(p0) <- rownames(SM)
    p0[seeds,] <- 1 / length(seeds) #p0 initializaton
    p <- SM%*%p0
    return(p)
}
##########################################################################################################################
##this part is genetic algorithm, which is embeded with random walk with restart
##fff1 is the function for evaluating fitness index of subset of candidate genes 
fff1 <- function(driver_gene, DFscore, Cgenes, genesets,SM) {
    ##driver_gene: the subset of candidate genes to be evaluated in a cancer individual
    ##DFscore: the Escores of dysfunctional cancer hallmarks identify by GSEA
    ##cgenes:the common genes in ppi network and expression profiles
    ##genesets: the gene sets of cancer hallmarks
    ##SM:the stable Similarity matrix of ppi network when random walk with restart stop
    P <- RWR1(driver_gene,SM)
    P <- sort(P[Cgenes, ],decreasing = T)
    Hscore1 <- GSEA(P, TERM2GENE = genesets, verbose = FALSE, pvalueCutoff = 1)@result[names(DFscore),"NES"]

    
    COR <- cor(Hscore1, DFscore)
    result <- list(Hscore = cbind(Hscore1,DFscore), cor = COR)
    return(result)

}

##GA_F2 is the program for genetic algorithm
GA_F2 <- function(genes,DFscore, genesets, SM, popsize, pcrossover, pmutation, elitism, maxiter, parallel = NULL, Cgenes ) {
##genes is all candidate genes in the cancer individual
##DFscorethe Escores of dysfunctional cancer hallmarks identify by GSEA
##geneset: the gene sets of cancer hallmarks
##Cgenes:the common genes in ppi network and expression profiles
    genes <- intersect(genes, rownames(SM))
    if (length(genes) == 0) {
        cat("There are no genes in ppi network")
        return(NA)
    }

    GA <- ga(type = "binary", fitness = fff, nBits = length(genes), inital_gene = genes, Cgenes=Cgenes,DFscore = DFscore, SM = SM, genesets = genesets, popSize = popSize,
         names = genes, pcrossover = pcrossover, pmutation = pmutation, elitism = elitism, maxiter = maxiter, monitor = T, parallel = parallel)
    driver_gene <- names(which((GA@solution)[1,] == 1))
    result <- fff1(driver_gene, DFscore, Cgenes, genesets, SM)
    result$driver <- driver_gene
    result$GA<-GA
    return(result)
}
############################################################################################################################
##driver_set1 is the main function to identify the personlized driver gene sets in cancer individual
driver_set1 <- function(path,Dname, flag, eflag, cnv_M, mut_M, meth_M, gene_M, gene_normal, SM, genesets,Cdrivergenes, popsize, pcrossover, pmutation, elitism, maxiter, parallel = F,pthr) {
   ##path:the output path
   ##Dname:the disease name
   ##flag: indicate what type of data type for identifying the driver genes for a cancer individual,c for copy number,m for mutation,
   ########me for methylation,cm for copy number and mutation,cme for copy number and methylation,mme for mutation and methylation,
   ########cmme for copy number,muation and methylation
   ##eflag: indicate whether expression profiles of cancer samples or normal sample were used for calculating fold change of genes
   #######T for cancer samples和C for normal samples
   ##cnv_M: 0-1 matrix of copy number alteration of genes across cancer samples
   ##mut_M: 0-1 matrix of mutation of genes across cancer samples
   ##meth_M: 0-1 matrix of hyper/hypo-methylation of genes across cancer samples
   ##gene_M: the expression profile of cancer samples
   ##gene_normal:the expression profile of normal samples
   ##SM:the stable Similarity matrix of ppi network when random walk with restart stop
   ##Cdrivergenes
   ##popsize: population size in genetic algorithm
   ##pcrossover:crossover operator in genetic algorithm
   ##pmutation: mutation operator in genetic algorithm
   ##elitism: the number of individuals with high fitness to be retain 
   ##maxiter: the number of iterations in genetic algorithm
   ##parallel: logical value, whether canculate in parallel
   ##pthr: the significance threshold for identifying the dysfunctional cancer hallmark
    path <- paste(path, Dname, "/", flag, sep = "")
    if (!dir.exists(path)) {
        dir.create(path, recursive = T)
    }
    setwd(path)

    if (flag == "c") {
        patients <- colnames(cnv_M)
    } else if (flag == "m") {
        patients <- colnames(mut_M)

    } else if (flag == "me") {
        patients <- colnames(meth_M)

    } else if (flag == "cm") {
        patients <- intersect(colnames(cnv_M), colnames(mut_M))
    } else if (flag == "cme") {
        patients <- intersect(colnames(cnv_M), colnames(meth_M))
    } else if (flag == "mme") {
        patients <- intersect(colnames(mut_M), colnames(meth_M))

    } else if (flag == "cmmme") {
        patients <- intersect(colnames(cnv_M), intersect(colnames(mut_M), colnames(meth_M)))

    }

    if (eflag == "C") {
        FC_M <- FC_exp(eflag, case = gene_M, normal = gene_normal)
    } else if (eflag == "T") {
        FC_M <- FC_exp(eflag, case = gene_M)

    }

    for (patient_id in patients) {
        print(patient_id)

        if (flag == "c") {
            cnv <- cnv_M[, patient_id]
            genes <- GA_genes(flag, cnv = cnv)
        } else if (flag == "m") {
            mut <- mut_M[, patient_id]
            genes <- GA_genes(flag, mut = mut)

        } else if (flag == "me") {
            meth <- meth_M[, patient_id]
            genes <- GA_genes(flag, meth = meth)

        } else if (flag == "cm") {
            cnv <- cnv_M[, patient_id]
            mut <- mut_M[, patient_id]
            genes <- GA_genes(flag, cnv = cnv, mut = mut)
        } else if (flag == "cme") {
            cnv <- cnv_M[, patient_id]
            meth <- mut_M[, patient_id]
            genes <- GA_genes(flag, cnv = cnv, meth = meth)

        } else if (flag == "mme") {
            mut <- mut_M[, patient_id]
            meth <- meth_M[, patient_id]
            genes <- GA_genes(flag, mut = mut, meth = meth)

        } else if (flag == "cmmme") {
            cnv <- cnv_M[, patient_id]
            mut <- mut_M[, patient_id]
            meth <- meth_M[, patient_id]
            genes <- GA_genes(flag, cnv = cnv, mut = mut, meth = meth)

        }
        genes<-intersect(genes,Cdrivergenes)

        FC <- FC_M[, patient_id, drop = F]
        Cgene <- intersect(rownames(FC), rownames(SM))
        FC<-sort(FC[Cgene,], decreasing = T)
        DFscore <- GSEA(FC, TERM2GENE = genesets, verbose = FALSE, pvalueCutoff = 1)@result
        DF<-DFscore$NES
        names(DF)<-DFscore$ID
        DFscore<-DF[DFscore$pvalue<=pthr]

        gaControl("binary" = list(selection = "gabin_tourSelection", crossover = "gabin_uCrossover"))
        result <- GA_F2(genes, FC, genesets, SM, popsize, pcrossover, pmutation, elitism, maxiter, parallel = parallel,Cgenes=Cgenes,DFscore=DFscore)
        file <- paste(patient_id, ".Rdata", sep = "")
        print(file)
        save(result, file = file)
    }
}
