
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

CollectPatientList <- function(gene, patients, n.pat)
  {
    ## Warning: Do not use. 
    require("assertthat")
    pat.subset <- patients[patients$Hugo_Symbol==gene,]
    index1 <- pat.subset$Validation_Status == "Unknown" &
                         pat.subset$Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',
                                                                'In_Frame_Ins', 'Missense_Mutation','Nonsense_Mutation',
                                                                'Translation_Start_Site', 'Nonstop_Mutation',
                                                                'De_novo_Start_OutOfFrame')
    index2 <- pat.subset$Validation_Status == "Valid"

    browser()
    patient.ids <- unique(pat.subset$patient)

    ##Validate that number of patients is the same as n.pat
    are_equal(length(patient.ids), n.pat)

    return(patient.ids)
    
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

RunSVM <- function(dd.train, permute=FALSE)
  {
    require(e1071)
    ##remove all data points with no mutations
    age.index <- grep("CLI_age", colnames(dd.train))
    index <- apply(dd.train, 1, function(x) all(is.na(x[-age.index])))
    dd.train <- dd.train[!index,]

    ## turn mutation numbers to binary values
    dd.train[,-age.index] <- dd.train[,-age.index] > 0 
    if(permute){
      dd.train[,"CLI_age"] <- sample(snp.dat[,"CLI_age"],size=length(dd.train[,"CLI_age"]), replace=FALSE)
      ttl <- 'Classification of randomized age group'
    } else{
      ttl <- 'Classification by patient age'
    }

    label <- rep('old', nrow(dd.train))
    label[dd.train$CLI_age<=50] <- 'young'
    label <- factor(label)
    dd.train <- cbind(label, dd.train)

    model <- svm(label ~., data=dd.train[,-grep("CLI_age", colnames(dd.train))], type='C-classification',
                 kernel='linear', scale=TRUE, cost=0.1, cross=10)

    ## same as model$fitted
    pred <- predict(model, data=dd.train[,-grep("CLI_age", colnames(dd.train))])
    ## Check agreement with prediction
    tab <- table(pred=pred, true=label)
    return(tab)
    
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

#########################
## DO NOT USE: Does not output
## correct results
#########################
##1. Collect significantly mutated genes
## genes <- CollectSigMutatedGenes(mutation.dir)

##2. Identify the patients with nonsilent mutations in the collected genes
## pat <- CollectPatientMutData(mutation.dir)
## TODO : Need to identify the patient ids with mutated genes using 'pat'
## 'pat' contains all mutations so need to filter by Variant_classification for
## functional mutations.
## As validation of correct filtering the 'npat' field in genes has to match the unique number
## of patient ids that carry mutations in this gene.
## Build list(gene1=c(patient ids), gene2=c(patient-ids))
## gene.patient.list <- lapply(genes$gene, function(x) CollectPatientList(x, pat, genes[genes$gene==x,'npat']) )
## names(gene.patient.list) <- genes$gene
##3. Collected Patients clinical data.
##4. Build a matrix of n- patients m-mutated genes and last column is patient age

################################

## Collect Aggregate data
dd <- CollectAggregateData(aggregate.analysis.dir)

somatic.mutation.genes <- colnames(dd)[grep("SMG_mutsig.2CV", colnames(dd))]
dd.reduce <- dd[,c("CLI_age", somatic.mutation.genes)]
dd.reduce <- dd.reduce[!is.na(dd.reduce[,"CLI_age"]),]
## RunSVM(dd.reduce)

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

