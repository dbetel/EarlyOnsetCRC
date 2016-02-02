library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
source("multiplot.R")

## Load and join CRC and Risk factor data
crc <- read.table("../data/CRC_county_clinical_2015-10-20.csv", sep='\t', header=TRUE, stringsAsFactors=FALSE, comment.char='#',
                  colClass=c('character','character','character','numeric','numeric','numeric','numeric','numeric','numeric'))

cdc <- read.csv("../data/Risk_factors_rates.csv", sep=',', header=TRUE, stringsAsFactors=FALSE,
                colClass=c('character','character','character','numeric','numeric','numeric','numeric'))


colnames(crc) <- gsub('\\.+', '_', colnames(crc))
colnames(cdc) <- gsub('\\.+', '_', colnames(cdc))
dd <- dplyr::inner_join(crc, cdc, by="FIPS")


plot2file <- FALSE
if(plot2file){
    pdf(paste0('../results/YoungCRCEpidemiology_', Sys.Date(), '.pdf'))
}


## Group CRC counts by anatomical location and sum the age groups
dd1 <- group_by(dd,  Anatomical_location) %>%
  select(Counts_Young, Counts_Old, Young_Population, Old_Population) %>%
  summarize( Old = sum(Counts_Old), Young=sum(Counts_Young),
            Young_Population=sum(Young_Population), Old_Population=sum(Old_Population))

## Chi-sequare test for differences in anatomical location
## between young and old.
## transform dd1 - Note that chiseq is same if dd1 is defined as follows
## clnames <- select(dd1,Anatomical_location)
## dd1 <- t(select(dd1,Old, Young))
## colnames(dd1) <- clnames$Anatomical_location
ct1 <- chisq.test(as.data.frame(dd1[, c('Old', 'Young')]))
summary(ct1)

## add fraction column 
dd1 <- mutate(dd1, fraction=(Young/Young_Population)/ (Old/Old_Population))

## Remove  C18.9-Colon, NO and C19.9-Rectosigmoid junction locations
dd1 <- dd1[-grep('^C18.[89]',dd1$Anatomical_location),]

p1 <- ggplot(as.data.frame(dd1), aes(x=Anatomical_location, y=fraction, fill=Anatomical_location)) + geom_bar(stat='identity')
p1 <- p1 + scale_fill_brewer(palette="Spectral")
p1 <- p1 + theme(axis.text.x = element_text(angle=45, hjust=1),
               panel.background = element_blank(),
               panel.grid.major= element_line(color='gray', linetype=2, size=0.5), legend.position='none')
p1 <- p1 + labs(title = "Ratios of E-CRC rates over T-CRC rates by anatomical location" , x = 'Anatomical Location',
              y='E-CRC/T-CRC rates (age-adjusted)')
print(p1)

## Group CRC by county and calculate mean of Smokers, Obese, Drinking and take sum of age groups
## Young/old summary should be sum by county
dd2 <- group_by(dd, State, County.y) %>%
  select(FIPS, State, County.x,Counts_Old, Counts_Young,Young_Population,Old_Population, X_Smokers, X_Obese, X_Excessive_Drinking,X_Diabetes_Age_Adjusted) %>%
  summarize(smokers = median(X_Smokers), obese=median(X_Obese), drinking=median(X_Excessive_Drinking),
            diabetes = median(X_Diabetes_Age_Adjusted),
            young = sum(Counts_Young), old = sum(Counts_Old),
            Young_Population=median(Young_Population), Old_Population=median(Old_Population))

dd2 <-  mutate(dd2, fraction = (young/Young_Population)/(old/Old_Population))

PlotClinical <- function(ddx, clinical){
  ddx <- filter(ddx, !is.infinite(fraction))
  p <- ggplot(ddx, aes_string(x=clinical, y='fraction')) + geom_point(shape=16, na.rm=TRUE) + stat_smooth(method=loess, na.rm=TRUE)
  p <- p + labs(y='Age adj. young/old rates')
  return(p)
}

PlotYoungOldClinical <- function(ddx, clinical ,isYoung=TRUE){
  ddx <- select(ddx, State, County.y, smokers, obese, diabetes, drinking, young, old,
                Young_Population, Old_Population)

  if(isYoung){
    ddx <- mutate(ddx, rate = young/Young_Population)
  } else{
    ddx <- mutate(ddx, rate = old/Old_Population)
  }
  
  fit <- lm(as.formula(paste('rate ~' ,clinical)), data=ddx)
  
  p <- ggplot(ddx,aes_string(x=clinical, y='rate')) + geom_point(shape=16, na.rm=TRUE) + stat_smooth(method=lm, na.rm=TRUE)
  xxlab <- ifelse(isYoung, paste(clinical, 'young'), paste(clinical, 'old'))
  xvals <- as.data.frame(select_(ddx, .dots=list(clinical)))[,2]
  
  lbl <- paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 4),
                    "\nSlope =",signif(fit$coef[[2]], 4),
                     "\np =",signif(summary(fit)$coef[2,4], 4))
  p <- p + annotate('text', label=lbl,
                    x= min(xvals,na.rm=TRUE) + diff(range(xvals, na.rm=TRUE))/8 ,
                    y=max(ddx[,"rate"]) - diff(range(ddx[,"rate"]))/9 ,
                    size=2.0) 
  p <- p + labs(x=xxlab, y='Age adj. rate') + theme_bw()
  
  return(p)
    
}

plots <- lapply(c('smokers', 'obese','diabetes', 'drinking'), function(x) PlotClinical(dd2, x))

multiplot(plotlist=plots, cols=2)

youngplots <- lapply(c('smokers', 'obese',  'diabetes', 'drinking'), function(x) PlotYoungOldClinical(dd2, x, isYoung=TRUE))
oldplots <- lapply(c('smokers', 'obese', 'diabetes','drinking'), function(x) PlotYoungOldClinical(dd2, x, isYoung=FALSE))

multiplot(plotlist=c(youngplots, oldplots), cols=2)

if(plot2file){
  grid.newpage()
  library(gridExtra)
  grid.table(as.data.frame(dd1))  
  dev.off()
}
