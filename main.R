GA_F <- function(genes, FC, genesets, W, r, u, popsize, pcrossover, pmutation, elitism, maxiter, parallel = F) {
    genes <- intersect(genes, rownames(W))
    if (length(genes) == 0) {
        cat("There are no genes in ppi network")
        return(NA)

    }
    GA <- ga(type = "binary", fitness = f, nBits = length(genes), inital_gene = genes, FC = FC, W = W, r = r, u = u, genesets = genesets, popSize = popSize,
         names = genes, pcrossover = pcrossover, pmutation = pmutation, elitism = elitism, maxiter = maxiter, monitor = F)
    driver_gene <- summary(GA)[11][[1]]
    result <- f1(driver_gene, FC, genesets, W, r, u)
    return(result)
}

GA_F1 <- function(genes, FC, genesets, SM,Cgene, popsize, pcrossover, pmutation, elitism, maxiter, parallel = NULL) {
    genes <- intersect(genes, rownames(SM))
    if (length(genes) == 0) {
        cat("There are no genes in ppi network")
        return(NA)

    }
    #common_gene<-intersect(rownames(FC),rownames(SM))
    #genes<-intersect(genes,rownames(SM))
    GA <- ga(type = "binary", fitness = ff, nBits = length(genes), inital_gene = genes, FC = FC, SM = SM, genesets = genesets,common_genes=Cgene,popSize = popsize,
         names = genes, pcrossover = pcrossover, pmutation = pmutation, elitism = elitism, maxiter = maxiter, monitor = F, parallel = parallel)
    driver_gene <- names(which((GA@solution)[1,] == 1))
    result <- ff1(driver_gene, FC, genesets, SM,Cgene)
    result$driver<-driver_gene
    return(result)
}




driver_set <- function(Dname,flag, eflag, cnv_M, mut_M, meth_M, gene_M, gene_normal, SM,gene_h,genesets_all,Cdrivergenes,pthr,pcrossover, pmutation, elitism) {
    path <- paste("2.driver genesets/data/", Dname, "/", flag, sep = "")
    if (!dir.exists(path)) {
        dir.create(path,recursive = T)
    }
    setwd(path)

    if (flag == "c") {
        patients <- intersect(colnames(cnv_M),colnames(gene_M))
    } else if (flag == "m") {
        patients <- intersect(colnames(mut_M), colnames(gene_M))

    } else if (flag == "me") {
        patients <- intersect(colnames(meth_M),colnames(gene_M))

    } else if (flag == "cm") {
        patients <- intersect(intersect(colnames(cnv_M), colnames(mut_M)),colnames(gene_M))
    } else if (flag == "cme") {
        patients <- intersect(intersect(colnames(cnv_M), colnames(meth_M)),colnames(gene_M))
    } else if (flag == "mme") {
        patients <- intersect(intersect(colnames(mut_M), colnames(meth_M)),colnames(gene_M))

    } else if (flag == "cmmme") {
        patients <- intersect(intersect(colnames(cnv_M), intersect(colnames(mut_M), colnames(meth_M))),colnames(gene_M))

    }

    if (eflag == "C") {
        FC_M <- FC_exp(eflag, case = gene_M, normal = gene_normal)
    } else if (eflag == "T") {
        FC_M <- FC_exp(eflag, case = gene_M)

    }
    #patients<-setdiff(patients,pp)
    for (patient_id in patients) {
        print(patient_id)
        genesets<-genesets_all
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
        file <- paste(patient_id, ".Rdata", sep = "")
        genes<-intersect(genes,Cdrivergenes)
        if(length(genes)==0){
           result=NULL
          save(result,file=file)
          next
          }
        FC <- FC_M[, patient_id, drop = F]
        
        Cgene <- intersect(rownames(FC), rownames(SM))
        FC1<-sort(FC[Cgene,], decreasing = T)
        GSEA_FC <- GSEA(FC1, TERM2GENE =gene_h, verbose = FALSE, pvalueCutoff = 1)
        save(GSEA_FC,file=paste(patient_id,"_GSEA_FC.Rdata",sep=""))
        DFscore<-GSEA_FC@result
       # names(DF)<-DFscore$ID
        Hname <- rownames(DFscore[DFscore$pvalue <= pthr & DFscore$NES,])
        if(length(Hname)>=3){
          
            genesets<- lapply(lapply(DFscore[Hname, "core_enrichment"], strsplit, split = "/"), unlist)##显著富集在癌症hallmark的核心基因集
            names(genesets) <- Hname
        }
        l<-length(genes)
        if (l== 1) {
            result <- ff1(genes, FC, genesets, SM, Cgene)
            result$driver <- genes
            save(result, file = file)
            next

        } else if (l<= 6) {
            popsize<-l
            maxiter <- 1
            parallel<-F
        } else if (l<=10) {
            popsize <- 20
            maxiter <- ceiling(2 ^ l / 20)
            parallel <- 4
        } else if (l<=15) {
            popsize <- 50
            maxiter <- 150
            parallel <- 4
        } else {
            popsize <- 100
            maxiter <- 150
            parallel <- 4
        }
        gaControl("binary" = list(selection = "gabin_tourSelection", crossover = "gabin_uCrossover"))
        result <- GA_F1(genes, FC, genesets, SM, Cgene,popsize, pcrossover, pmutation, elitism, maxiter, parallel = parallel)
        
       
        save(result, file = file)
    }
}

GA_F2 <- function(genes,DFscore, genesets, SM, popsize, pcrossover, pmutation, elitism, maxiter, parallel = NULL, Cgenes ) {
##genes是癌症个体的遗传改变基因集合
##DFscore是癌症个体hallmark激活得分
##geneset：是hallmark基因集合，是一个两列的d
##Cgenes是表达谱与互作网络共同基因
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
driver_set1 <- function(path,Dname, flag, eflag, cnv_M, mut_M, meth_M, gene_M, gene_normal, SM, genesets,Cdrivergenes, popsize, pcrossover, pmutation, elitism, maxiter, parallel = F,pthr) {
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