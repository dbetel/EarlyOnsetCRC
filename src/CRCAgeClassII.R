
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


WilcoxTest <- function(dd, gene)
  {
    res <- wilcox.test(dd[dd[,gene] == 0, "CLI_age"], dd[dd[,gene] == 1, "CLI_age"], alternative='two.sided', paired=FALSE)
    return(res$p.value)
  }

given.n <- function(x){
  return(c(y=median(x) * 1.05, label=length(x)))
}

BoxPlotGenesAge <- function(dd.plot)
  {
    require(reshape)
    require(ggplot2)
    
    colnames(dd.plot) <- c(gsub("SMG_mutsig.2CV_", '', colnames(dd.plot)))
    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "gene", "mutated")
    
    p <- ggplot(dd.plot, aes(x=mutated, y=age)) + geom_boxplot(color=c('blue', 'red'))
    p <- p + geom_jitter(aes(alpha=0.4), colour='gray50', show_guide=FALSE)
    p <- p + stat_summary(fun.data= given.n, geom='text', fun.y=median, size=3.55)
    p <- p + facet_wrap(~ gene, ncol=3)
    print(p)
    
  }

ECDFPlotGenesAge <- function(dd.plot, colrs)
  {
    require(reshape)
    require(ggplot2)
    
    colnames(dd.plot) <- c(gsub("SMG_mutsig.2CV_", '', colnames(dd.plot)))
    dd.plot <- melt(dd.plot, id=c("CLI_age"))
    colnames(dd.plot) <- c("age", "gene", "mutated")
    
    p <- ggplot(dd.plot, aes(age, colour=mutated)) + stat_ecdf()
    p <- p + scale_colour_manual(values=colrs, name="Status", labels=c("WT", "Mutated"))
    p <- p + facet_wrap(~ gene, ncol=3)
    print(p)
    
  }

###########
## Main
###########
data.dir <- "/zenodotus/dat01/betellab_store/dob2014/Shah/BROAData/"
mutation.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.MutSigNozzleReport2CV.Level_4.2014101700.0.0/")
patient.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD.Merge_Clinical.Level_1.2015020400.0.0/")
aggregate.analysis.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.Aggregate_AnalysisFeatures.Level_4.2014101700.1.0/")
colrs <- c("#000000", "#0000D4", "#009500", "#FFFF00", "#FF9100", "#E40000", "#546723")


## Collect Aggregate data
dd <- CollectAggregateData(aggregate.analysis.dir)

somatic.mutation.genes <- colnames(dd)[grep("SMG_mutsig.2CV", colnames(dd))]
dd.reduce <- dd[,c("CLI_age", somatic.mutation.genes)]
dd.reduce <- dd.reduce[!is.na(dd.reduce[,"CLI_age"]),]

##remove all data points with no mutations
age.index <- grep("CLI_age", colnames(dd.reduce))
index <- apply(dd.reduce, 1, function(x) all(is.na(x[-age.index])))
dd.reduce <- dd.reduce[!index,]

## turn mutation numbers to binary values
dd.reduce[,-age.index] <- dd.reduce[,-age.index] > 0 

sig.genes <- sapply(somatic.mutation.genes, function(x) WilcoxTest(dd.reduce[,c("CLI_age", x)], x))

dd.plot <- dd.reduce[,c("CLI_age", names(sig.genes[sig.genes <=0.1]))]

plot2file=TRUE
if(plot2file){
  pdf(paste0("../results/MutAge_", Sys.Date(), ".pdf"))
}

## BoxPlotGenesAge(dd.plot)
ECDFPlotGenesAge(dd.plot, colrs)

if(plot2file){
  dev.off()
}

