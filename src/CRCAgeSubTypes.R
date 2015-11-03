
CollectAggregateData <- function(base.dir)
  {
    ## Clinical features : start with CLI_ followed by clinical feature name ex: CLI_age.
    ## Clustering results : start with CLUS_ followed by platform_method ex: CLUS_mRNAseq_cHierarchical.
    ## Somatic mutation genes : start with SMG_ followed by version number( mutsig2.0,cv,2cv)_gene name ex: SMG_mutsig.2CV_FAM47C.
    ## Somatic mutation genes expression : start wit SMG_ followed by gene name_mRNA ex: SMG_KRT3_mRNA.
    ## Marker genes in each mRNAseq clustering subtype : star with mRNA_ followed by CNMF_gene name_difference_cluster number ex: mRNA_CNMF_FAM66E_.0.6_2 (In each cluster, the top 5 up regulated and top 5 down regulated genes were selected).
    ## Copy number alterations :
    ##     Copy number focal events : start with Amp_/Del_ followed by cytoband name ex: Amp_1q32.1/Del_1p36.32.
    ##     Copy number arm level events : start with CN_ followed by arm_Amp/Del ex: CN_10p_Amp/CN_10p_Del.
    ##     Copy number alterative gene with expression: start with Amp/Del_ followed by gene name, cytoband names which contains this alterative gene and _mRNA ex: Amp_SOX2_3q26.32_mRNA/Del_PARK2_6q24.3_mRNA.

    file.path <- paste0(base.dir, "COADREAD-TP.samplefeatures.txt")
    aggregate.features <- read.delim(file.path, header=TRUE, row.names=1,
                                     stringsAsFactors=FALSE, sep='\t', na.strings="NA")

    
    return(aggregate.features)
  }

CollectCNV <- function(base.dir){
  file.path <- paste0(base.dir, "amp_genes.conf_99.txt")
  cnv <- read.delim(file.path, header=T, stringsAsFactors=FALSE, sep='\t', na.strings='NA')
  return(cnv)
}

  
BoxPlotClusters <- function(dd.plot)
  {
    require(reshape)
    require(ggplot2)
    
    colnames(dd.plot) <- c(gsub("CLUS_", '', colnames(dd.plot)))
    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "type", "cluster")

    ## remove NA
    dd.plot <- dd.plot[!is.na(dd.plot$cluster),]

    ## generate factors
    dd.plot$cluster <- factor(dd.plot$cluster)
    p <- ggplot(dd.plot, aes(x=cluster, y=age)) + geom_boxplot(aes(colour=cluster))
    p <- p + geom_jitter(aes(alpha=0.4), colour='gray50',  show_guide=FALSE)
    p <- p + facet_wrap(~ type, ncol=3)
    print(p)
    
  }

ECDFClusters <- function(dd.plot, colrs)
  {
    require(reshape)
    require(ggplot2)
    
    colnames(dd.plot) <- c(gsub("CLUS_", '', colnames(dd.plot)))
    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "type", "cluster")

    ## remove NA
    dd.plot <- dd.plot[!is.na(dd.plot$cluster),]

    ## generate factors
    dd.plot$cluster <- factor(dd.plot$cluster)
    p <- ggplot(dd.plot, aes(colour=cluster, x=age)) + stat_ecdf() + scale_color_manual(values=colrs)
    p <- p + facet_wrap(~ type, ncol=3)
    print(p)
    
  }

BoxPlotCNVAge <- function(dd.plot, values)
  {
    require(reshape)
    require(ggplot2)
    
    ## convert CNV values to binary
    for(cc in colnames(values)){
      if(grepl("Amp", cc)){
        dd.plot[,cc] <- dd.plot[,cc] >= values['thd', cc]
      }else{
        dd.plot[,cc] <- dd.plot[,cc] <= values['thd', cc]
      }
    }
      
    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "CNV", "altered")

    ## trun binary value to factor
    dd.plot$altered <- factor(dd.plot$altered)
    
    p <- ggplot(dd.plot, aes(x=altered, y=age)) + geom_boxplot(color=c('blue', 'red'))
    p <- p + geom_jitter(aes(alpha=0.4), colour='gray50', show_guide=FALSE)
    p <- p + facet_wrap(~ CNV, ncol=3)
    print(p)
    
  }

ECDFCNVAge <- function(dd.plot, values, colrs)
  {
    require(reshape)
    require(ggplot2)
    for(cc in colnames(values)){
      if(grepl("Amp", cc)){
        dd.plot[,cc] <- dd.plot[,cc] >= values['thd',cc]
      }else{
        dd.plot[,cc] <- dd.plot[,cc] <= values['thd', cc]
      }
    }

    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "CNV", "altered")

    ## turn binary values to factor
    dd.plot$altered <- factor(dd.plot$altered)

    p <- ggplot(dd.plot, aes(colour=altered, x=age)) + stat_ecdf()
    p <- p + scale_colour_manual(values=colrs, name="Status", labels=c("WT", "CNV"))
    p <- p + facet_wrap(~ CNV, ncol=3)
    print(p)
    
  }

CorTest <- function(dd, gene)
  {
    ## Spearman rank correlation
    res <- cor(dd["CLI_age"], dd[,gene], method='spearman', use='complete.obs')
    return(res)
  }

CNVWilcoxTest <- function(dd, cnv)
  {
    ## Del or Amp? 
    del <- grepl("Del", cnv) 
    thd = quantile(dd.cn[,cnv], seq(0,1,0.1))["80%"]
    if(del){
      thd = quantile(dd.cn[,cnv], seq(0,1,0.1))["20%"]
    }
    ## plot(density(dd[,cnv], bw=0.5, na.rm=TRUE), main=cnv, col='red', lwd=2)
    ## abline(v=thd, lty=3, col='gray')
    if(!del){
      res <- wilcox.test(dd[dd[,cnv] <= thd, "CLI_age"], dd[dd[,cnv] > thd, "CLI_age"], alternative='two.sided', paired=FALSE)
    }
    else{res <- wilcox.test(dd[dd[,cnv] >= thd, "CLI_age"], dd[dd[,cnv] < thd, "CLI_age"], alternative='two.sided', paired=FALSE)}

    return(list(pvl=res$p.value, thd=thd))
  }


###########
## Main
###########
data.dir <- "/zenodotus/dat01/betellab_store/dob2014/Shah/BROAData/"
mutation.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.MutSigNozzleReport2CV.Level_4.2014101700.0.0/")
patient.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD.Merge_Clinical.Level_1.2015020400.0.0/")
aggregate.analysis.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.Aggregate_AnalysisFeatures.Level_4.2014101700.1.0/")
cnv.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0/")

## colrs <- c("#2121D9", "#9999FF", "#D92121", "#21D921", "#FFFF4D", "#FF9326")
colrs <- c("#000000", "#0000D4", "#009500", "#FFFF00", "#FF9100", "#E40000", "#546723")

plot2file=TRUE
if(plot2file)
{
  pdf(paste0("../results/SubType_", Sys.Date(), ".pdf"))
  
}

## Collect Aggregate data
dd <- CollectAggregateData(aggregate.analysis.dir)

#########################
## Clustering analysis
## TCGA pipeline uses non-negative matrix factorization and
## hierarchical clustering
##########################
cluster.names <- colnames(dd)[grep("CLUS", colnames(dd))]
dd.reduce <- dd[,c("CLI_age", cluster.names)]
dd.reduce <- dd.reduce[!is.na(dd.reduce[,"CLI_age"]),]


##remove all data points with no cluster assignments
age.index <- grep("CLI_age", colnames(dd.reduce))
index <- apply(dd.reduce, 1, function(x) all(is.na(x[-age.index])))
dd.reduce <- dd.reduce[!index,]

## BoxPlotClusters(dd.reduce)
ECDFClusters(dd.reduce, colrs)

#########################
## Copy number alterations
##########################
cn.names <- colnames(dd)[grep("(Amp|Del).*[0-9]+$", colnames(dd), perl=T)]
dd.cn <-  dd[,c("CLI_age", cn.names)]
dd.cn <- dd.cn[!is.na(dd.cn[,"CLI_age"]),]

##remove all samples with no CN values
age.index <- grep("CLI_age", colnames(dd.cn))
index <- apply(dd.cn, 1, function(x) all(is.na(x[-age.index])))
dd.cn <- dd.cn[!index,]

cnv.cor <- sapply(cn.names, function(x) CNVWilcoxTest(dd.cn[,c("CLI_age", x)], x))
dd.cn.sel <- cnv.cor[,cnv.cor['pvl',] <=0.1]
dd.plot <- dd.cn[,c("CLI_age", colnames(dd.cn.sel))]
## BoxPlotCNVAge(dd.plot, dd.cn.sel)
ECDFCNVAge(dd.plot,dd.cn.sel, colrs)

if(plot2file){
  dev.off()
}
