library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(grid)

source("multiplot.R")


LoadSheet <- function(filename, sheetNumber)
{
  require('xlsx')
  sheet <- read.xlsx(filename, sheetIndex=sheetNumber, startRow=2, header=TRUE, stringsAsFactors=FALSE)

  ## remove rows with no data
  sheet <- sheet[-which(is.na(sheet$Std.Pop)),]
  ## header names
  colnames(sheet)[1:3] <- c('Location', 'Age.Group', 'Year')
  colnames(sheet) <- gsub('\\.+', '_', colnames(sheet))
  return(sheet)
}

SummarizeAgeYearCRCRates <- function(df, eg)
{
  ## Summarize table by grouping age-group, and year, and sum over num of cases.
  ## dt = data_frame, eg = ethnic group
  df <- select(df, Location, Age_Group, Year, Count, Pop, Std_Pop) %>%
    group_by(Age_Group, Year) %>%
      summarize(count = sum(Count), population=unique(Pop), all_pop = unique(Std_Pop)) %>%
        ## mutate(Ethnic_Group= eg, rate=count/population)
        mutate(Ethnic_Group= eg)

  
  ## Add young and old class and summarize
  df2 <- mutate(df, Age = ifelse((Age_Group %in% c('20-24 years' ,'25-29 years' ,'30-34 years' ,
                        '35-39 years','40-44 years', '45-49 years')), 'E-CRC' , 'T-CRC')) %>%
    select(Age_Group, Year, count, population, all_pop, Age, Ethnic_Group) %>%
    group_by(Age, Year) %>%
    summarize(count = sum(count), population=sum(population), all_pop=sum(all_pop),
              Ethnic_Group=unique(Ethnic_Group), Age_Group = toString(Age_Group)) %>%
                mutate(rate=count/population)
  
  return(df2)
      
}

PlotData <- function(dd_x, regression_x, sf, x_pos, legend.title=NULL){
    ## dd_x - data to plot
    ## regression_x regression data_frame with Age, Etnic_Group, fit (lm model obj), Slope , Intercept
    ## sf = scaling factor
    ## x_pos = x-axis position (year)
    ## y_pos <-  regression_x %>% group_by(Age, Ethnic_Group) %>%
    ##     do(predict.lm(fit, newdata=data.frame(Year=x_pos)))
    y_pos <- regression_x %>% do(as.data.frame(predict.lm(.$fit, newdata=data.frame(Year=x_pos))))
    dd_x <- mutate(dd_x, rate=rate*sf)
    p <- ggplot(dd_x,aes( x=Year, y=rate, colour= Ethnic_Group, shape=Age))
    p <- p + geom_point(size=4)
    p <- p + scale_colour_hue(l=50)
    p <- p + geom_smooth(method=lm) ## note that regression is done on the colour + share grouping in ggplot
    p <- p + scale_color_manual(name=legend.title, values=c('#1b9e77','#d95f02', '#7570b3'))
    p <- p + guides(shape=FALSE) ## turn off shape legend
    p <- p + annotate('text', x=2003, y=c(round(y_pos))[[1]], size=6,
                      label=paste(regression_x$Age, regression_x$Ethnic_Group, format(round(regression_x$Slope, 3))))
    p <- p + theme_bw()

    p <- p + theme(aspect.ratio=1, axis.title.x=element_text(size=18, vjust=-0.6), axis.title.y=element_text(size=18),
                   axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
                   legend.position="bottom", legend.title= element_text(size=18, face='bold'),
                   legend.text=element_text(size=18, face='bold'))

    p <- p + labs(title=NULL,
                  y=paste("Age adjusted rates of",legend.title,"cases (per 100K)"))
    return(p)
    ## print(p)

}
############
## Main
###########
## Load and join CRC and Risk factor data
data.file <- '../data/Rates_2000-2011.xlsx'
black <- LoadSheet(data.file, 1)
hispanic <- LoadSheet(data.file, 2)
white <- LoadSheet(data.file, 3)

white <- SummarizeAgeYearCRCRates(white, "White")
black <- SummarizeAgeYearCRCRates(black, "Black")
hispanic <- SummarizeAgeYearCRCRates(hispanic, "Hispanic")

dd <- bind_rows(white, black, hispanic)

plot2file <- TRUE
if(plot2file){
    pdf(paste0('../results/CRC_Rates_2000To2011_', Sys.Date(), '.pdf'), width=20, height=15)
}

scale_factor <- 1e5
##Fit regression lines
regressions <- dd %>% group_by(Age, Ethnic_Group) %>% do(fit=lm(rate*scale_factor ~ Year, data=.)) %>%
  mutate(Slope= summary(fit)$coeff[2], Intercept = summary(fit)$coeff[1])


p <- PlotData(dd, regressions, scale_factor, 2005, 'CRC')
## print(p)
eCRC <- PlotData(filter(dd, Age=='E-CRC'), filter(regressions, Age=='E-CRC'), scale_factor, 2005, 'E-CRC')
tCRC <- PlotData(filter(dd, Age=='T-CRC'), filter(regressions, Age=='T-CRC'), scale_factor, 2005, 'T-CRC')
Title_multiplot(plotlist=list(eCRC, tCRC), cols=2,
          ttl="Trends of Incidence of E-CRC and T-CRC by Ethnic Group from 2000-2011")
if(plot2file){
  dev.off()
}
