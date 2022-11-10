GA_genes <- function(flag, cnv = NULL, mut = NULL, meth = NULL) {
    ##flag决定用什么类型数据获得候选驱动基因集：取值有c,m,me,cm,cme,mme,cmme
    ##cnv是拷贝数0-1向量
    ##mut是突变0-1向量
    ##meth是甲基化0-1向量
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

FC_exp <- function(eflag, case, normal = NULL,log=F) {

    ##eflag指示用癌症群体表达计算FC还是正常样本计算FC，取值T和C
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
##随机游走函数
RWR1 <- function(seeds,SM){
    ##seeds是用于随机游走的种子，要求完全出现在网络中
    ##r是重启概率
    ##u是稳态判断条件
    ##W是网络转移概率矩阵
    #初始概率向量p0
    p0 <- matrix(0, nrow = dim(SM)[1], ncol = 1)
    rownames(p0) <- rownames(SM)
    p0[seeds,] <- 1 / length(seeds) #p0的初始化
    p <- SM%*%p0
    return(p)
}


fff1 <- function(driver_gene, DFscore, Cgenes, genesets,SM) {
    ##该函数用于计算驱动基因集的适应值评价指标
    ##driver_gene表示一个个体中识别的驱动基因集合
    ##FC是指一个个体转录组基因表达变化fold change向量
    ##genesets表示癌症hallmark基因集
    ##W是用于随机游走的蛋白质网络的转移矩阵
    ##r表示随机游走重启概率
    ##u表示随机游走稳态阈值

    P <- RWR1(driver_gene,SM)
    P <- sort(P[Cgenes, ],decreasing = T)
    Hscore1 <- GSEA(P, TERM2GENE = genesets, verbose = FALSE, pvalueCutoff = 1)@result[names(DFscore),"NES"]

    
    COR <- cor(Hscore1, DFscore)
    result <- list(Hscore = cbind(Hscore1,DFscore), cor = COR)
    return(result)

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
