


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

CollectMethylData <- function(base.dir)
  {
    
    file.path <- paste0(base.dir, "COADREAD-TP.membership.txt")
    methy <- read.delim(file.path, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t", na.string="NA")
    return(methy)
    
  }
###########
## Main
###########
data.dir <- "/zenodotus/dat01/betellab_store/dob2014/Shah/BROAData/"
methylation.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.Methylation_Clustering_CNMF.Level_4.2014101700.0.0/")

aggregate.analysis.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.Aggregate_AnalysisFeatures.Level_4.2014101700.1.0/")


## Collect Aggregate data
dd <- CollectAggregateData(aggregate.analysis.dir)

## Collect sample mutation rates
meth <-  CollectMethylData(methylation.dir)

index <- sapply(gsub('-[0-9]+$','',rownames(meth)), function(x) grep(x, rownames(dd)))

index <- index[sapply(index, length)==1]
meth$age <- dd[index, "CLI_age"]
colrs <- scale_colour_brewer()
## convert cluster membership to factors
membership_names <- grep("membership", colnames(meth))
meth[,membership_names] <- lapply(meth[,membership_names], factor)


given.n <- function(x){
  return(c(y=median(x) * 1.05, label=length(x)))
}

plot2file <- TRUE
if(plot2file){
  pdf(paste0("../results/CRCAgeMethylation_", Sys.Date(), ".pdf"))
}

## Plot
PlotMethylclusters <- function(meth.dd, clust){
  require(ggplot2)
  require(scales)
  p <- ggplot(meth, aes_string(x=clust, y="age")) + geom_boxplot(aes_string(colour=clust),show_guide=FALSE) + scale_colour_brewer(palette="Set1")
  P <- p + scale_y_continuous(breaks=pretty_breaks(n=10))
  p <- p + geom_jitter(aes(alpha=0.4),
                     colour='gray60',
                       show_guide=FALSE)
  p <- p + stat_summary(fun.data= given.n, geom='text', fun.y=median)
  labs(title="Methylation clusters by age", x="clusters", y="Age")
  print(p)
}

sapply(colnames(meth)[membership_names], function(x) PlotMethylclusters(meth, x))

if(plot2file){
  dev.off()
}
