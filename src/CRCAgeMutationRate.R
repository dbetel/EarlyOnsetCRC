
CollectSigMutatedGenes <- function(base.dir)
  {
    ## Column Description
    ## N = number of sequenced bases in this gene across the individual set
    ## n = number of (nonsilent) mutations in this gene across the individual set
    ## npat = number of patients (individuals) with at least one nonsilent mutation
    ## nsite = number of unique sites having a non-silent mutation
    ## nsil = number of silent mutations in this gene across the individual set
    ## n1 = number of nonsilent mutations of type: *CpG->T
    ## n2 = number of nonsilent mutations of type: *Cp(A/C/T)->mut
    ## n3 = number of nonsilent mutations of type: A->mut
    ## n4 = number of nonsilent mutations of type: *CpG->(G/A)
    ## n5 = number of nonsilent mutations of type: indel+null
    ## n6 = number of nonsilent mutations of type: double_null
    ## p_classic = p-value for the observed amount of nonsilent mutations being elevated in this gene
    ## p_ns_s = p-value for the observed nonsilent/silent ratio being elevated in this gene
    ## p_cons = p-value for enrichment of mutations at evolutionarily most-conserved sites in gene
    ## p_joint = p-value for clustering + conservation
    ## p = p-value (overall)
    ## q = q-value, False Discovery Rate (Benjamini-Hochberg procedure)

    file.path <- paste0(base.dir,"sig_genes.txt" )
    mutated.genes <- read.table(file.path, header=TRUE, stringsAsFactors=FALSE, sep='\t')
    mutated.genes <- mutated.genes[mutated.genes[,'q'] < 0.1, ]
    return(mutated.genes)
  }

CollectClinicalData <- function(base.dir)
  {

    file.path <- paste0(base.dir, "COADREAD.merged_only_clinical_clin_format.txt")
    patients <- read.table(file.path, header=TRUE, stringsAsFactors=FALSE, sep='\t')
    return(patients)
  }

CollectPatientMutData <- function(base.dir)
  {
    ## see https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
    ## for description of Mutation Annotation Format
    file.path <- paste0(base.dir, "COADREAD-TP.final_analysis_set.maf")
    patient.mutation <- read.delim(file.path, header=TRUE, stringsAsFactors=FALSE, sep='\t')

    ## filter by Variant_Classification == "Silent"
    return(patient.mutation)
  }



CollectPatientCountRates <- function(base.dir)
  {
    ## name = patient id
    ## nmut = number of mutation
    ## n_coding = number of coding mutation
    ## n_flank  
    ## n_coding_nonsilent = number of non-silent mutations
    ## fracflank
    ## callscheme_name = ["coding only", "exome+100bp flanks"]
    ## cov_idx 
    ## N_tot = Grand total from all samples ?
    ## N_c
    ## nsil_tot = non-silent total
    ## nnon_tot = nonsymnon total ?
    ## n_tot = total
    ## rate_sil= nsil_tot/N_tot
    ## rate_non = nnon_tot/N_tot
    ## rate_tot = n_tot/N_tot
    ## log_rate_tot = log10(n_tot/N_tot)
    ## n_c
    ## n_ind - indels ?
    ## N_ind
    ## rate_c
    ## rate_ind

    file.path <- paste0(base.dir, "patient_counts_and_rates.txt")
    patient.mutation <- read.delim(file.path, header=TRUE, stringsAsFactors=FALSE, sep='\t')
    return(patient.mutation)
  }


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

CollectCNVData <- function(base.dir)
  {
    ## The all lesions file summarizes the results from the GISTIC run.
    ## It contains data about the significant regions of amplification
    ## and deletion as well as which samples are amplified or deleted in each of these regions.
    ## The identified regions are listed down the first column, and the samples are listed across the first row, starting in column 10.

    ##  Columns 1-9 present the data about the significant regions as follows:
    ## Unique Name: A name assigned to identify the region.
    ## Descriptor: The genomic descriptor of that region.
    ## Wide Peak Limits: The 'wide peak' boundaries most likely to contain the targeted genes. These are listed in genomic coordinates and marker (or probe) indices.
    ## Peak Limits: The boundaries of the region of maximal amplification or deletion.
    ## Region Limits: The boundaries of the entire significant region of amplification or deletion.
    ## Q values: The Q value of the peak region.
    ## Residual Q values: The Q value of the peak region after removing ('peeling off') amplifications or deletions that overlap other, more significant peak regions in the same chromosome.
    ## Broad or Focal: Identifies whether the region reaches significance due primarily to broad events (called 'broad'), focal events (called 'focal'), or independently significant broad and focal events (called 'both').
    ## Amplitude Threshold: Key giving the meaning of values in the subsequent columns associated with each sample.
    
    file.path <- paste0(base.dir, "all_lesions.conf_99.txt")
    all.lesions <- read.delim(file.path, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t", na.string="NA", check.names=FALSE)
    return(all.lesions)
    
  }
###########
## Main
###########
data.dir <- "/zenodotus/dat01/betellab_store/dob2014/Shah/BROAData/"
mutation.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.MutSigNozzleReport2CV.Level_4.2014101700.0.0/")
patient.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD.Merge_Clinical.Level_1.2015020400.0.0/")
aggregate.analysis.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.Aggregate_AnalysisFeatures.Level_4.2014101700.1.0/")
gistic.dir <- paste0(data.dir, "gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0/")

## Collect Aggregate data
dd <- CollectAggregateData(aggregate.analysis.dir)

## Collect sample mutation rates
rates <-  CollectPatientCountRates(mutation.dir)

## find sample overlap
index <- sapply(rownames(dd), function(x) grep(x, rates$name))
temp <- sapply(index, length)
index <- index[temp==1]
rates[,'age'] <- 0
rates[,'id'] <- names(index)
rates[unlist(index),'age'] <- dd[names(index), 'CLI_age']

## silent mutation rates
PlotCorrelations <- function(x,y, title, xlab, ylab)
  {
    ct <- cor.test(x,y, method='spearman')
    ct.adj <- p.adjust(ct$p.value)
    mdl <- lm(y ~ log(x))
    plot(log(x), y, pch=19, col='gray', main=title, xlab=xlab, ylab=ylab)
    abline(reg=mdl, col='red', lty=2)
    txt <- paste("Spearman corr =",ct$estimate, "\n")

    legend("topright",
           legend=do.call("expression", list(bquote(Spearman==.(round(ct$estimate, 5))),
               bquote(p-val==.(round(ct$p.value, 4)))
               )))
    ## mtext(txt, side=3, padj=2)
  }

plot2file <- FALSE
if(plot2file){
  pdf(paste0("../results/CRCAgeMutation_", Sys.Date(), ".pdf"))
}

PlotCorrelations(rates$rate_sil,rates$age, "Age vs. Silent mutation rates", xlab="Silent mutation rates (log scale)", ylab="age")
PlotCorrelations(rates$rate_non,rates$age, "Age vs. Nonsynonymous mutation rates", xlab="Nonsynonymous mutation rates (log scale)", ylab="age")
PlotCorrelations(rates$rate_tot,rates$age, "Age vs. Total mutation rates", xlab="Total mutation rates (log scale)", ylab="age")


###############
## CNV
###############
cnv <- CollectCNVData(gistic.dir)
cnv.column.index <- sapply(rownames(dd), function(x) {
    ## match rownames to a column in cnv.
    index <- grep(x, colnames(cnv))
})

cnv.column.index <- cnv.column.index[sapply(cnv.column.index, length) == 1]

cnv.row.index <- grep("CN values", rownames(cnv), invert=TRUE)
cnv.counts <- sapply(names(cnv.column.index), function(x) {## count all CNV events
                                                           sum(cnv[cnv.row.index,cnv.column.index[[x]]] > 0 )})

cnv.data <- cbind(dd[names(cnv.counts), "CLI_age"], cnv.counts)
colnames(cnv.data) <- c("Age", "CNV_events")
cnv.data <- cnv.data[!is.na(cnv.data[,"Age"]),]
plot(cnv.data[,2], cnv.data[,"Age"], pch=19, col='gray',
     main="CNV events vs. Aga", xlab="number of CNV events", ylab="Age")


cnv.data <- as.data.frame(cnv.data, stringsAsFactors=FALSE)
cnv.data$group <- "Old"
cnv.data$group[cnv.data$Age <=50] <- "Young"

given.n <- function(x){
  return(c(y=median(x) * 1.05, label=length(x)))
}


require(ggplot2)
require(scales)
p <- ggplot(cnv.data, aes(x=group, y=CNV_events)) + geom_boxplot(colour=c('blue', 'red'))
p <- p + geom_jitter(aes(alpha=0.4), colour='gray60', show_guide=FALSE)
P <- p + scale_y_continuous(breaks=pretty_breaks(n=10))
p <- p + stat_summary(fun.data= given.n, geom='text', fun.y=median)
p <- p + labs(title="CNV events by age", x=NULL, y="Number of CNV events")
print(p)


if(plot2file){
  dev.off()
}

