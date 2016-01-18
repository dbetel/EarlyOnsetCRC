
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
require("dplyr")
require("tidyr")

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

plot_dd <- select(rates, name, nmut, rate_tot, rate_sil, rate_non, age) %>%
    mutate(age_range=cut(age, breaks =c(20,50,60,70,80,100), labels = c('20-49', '50-59', '60-69', '70-79', '80<='),
             right=FALSE)) %>%
    gather(key=mutation_type, value=mutation_rate, rate_sil,rate_non,rate_tot, na.rm=TRUE)

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

plot2file <- TRUE
if(plot2file){
  pdf(paste0("../results/CRCAgeMutationII_", Sys.Date(), ".pdf"))
}

PlotAgeMutDist <- function(ddx){
  require(ggplot2)
  require(scales)

  p <- ggplot(ddx, aes(x=age_range, y=1e6*mutation_rate, colour=mutation_type,fill=mutation_type ))
  p <- p + geom_point(position=position_jitterdodge(dodge.width=0.9), show_guide=FALSE)
  p <- p + geom_boxplot(fill='white', position=position_dodge(width=0.9), outlier.shape = NA, alpha=0.5)
  p <- p + scale_y_log10(breaks=c(1,10,100))
  p <- p + scale_colour_brewer(palette= 'Set1' ,type='qual', labels=c('Silent', 'nonsynonymous', 'total'))
  p <- p + theme_bw() + labs(title='Mutation rate by age group', x= 'Age group', y=expression('Mutation rate per million bases'))

  print(p)

}

PlotAgeMutDist(plot_dd)

## add hypermutated 'age_group'
## According to TCGA Nature, 487, 19 July 2012 hypermutated phenotype
## is defined as >12 per million bases
plot_dd2<- plot_dd
levels(plot_dd2$age_range) <- c(levels(plot_dd2$age_range),'hypermutated')
plot_dd2$age_range[plot_dd2$mutation_rate*1e6 >12] = 'hypermutated'

PlotAgeMutDist(plot_dd2)

if(plot2file){
  dev.off()
}

