if(as.double(R.Version()$minor) != 1.3){
  stop("Wrong version of R. Requires R3.1.3", call.=FALSE)
}


library(dplyr)
library(ggplot2)
library(ggcounty)
library(maptools)
library(gpclib)
gpclibPermit()
source("multiplot.R")
source("LoadPopCensus.R")

PlotMap <- function(dd, fill_field, sfg, ttl)
  {
    ## plot map using data dd field fill_field and
    ## scale_fill_gradient sfg
    ## 'dd' must have FIPS column
    us <- ggcounty.us()
    gg <- us$g 

    gg <- gg + geom_map(data=dd, map=us$map,
                        aes_string(map_id='FIPS', fill=fill_field),
                        color="gray", size=0.125)
    gg <- gg + theme(legend.text = element_text(size=9))
    gg <- gg + sfg + labs(title=ttl)
    return(gg)
    
  }


PlotState <- function(dd, state, fill_field, sfg, ethnic)
  {
    st <- ggcounty(state, fips=TRUE)
    gg <- st$gg
    print(state)
    
    gg <- gg + geom_map(data=dd, map=st$map,
                        aes_string(map_id='FIPS', fill=fill_field), color='gray')
    gg <- gg + sfg
    gg <- gg + theme(legend.position='bottom', legend.box='horizontal',
                     legend.text = element_text(size = 9, angle=45),
                     legend.key.size = unit(1.5, 'lines'))
    gg <- gg + labs(title=paste(ethnic, state, fill_field))
    return(gg)
  }


CRCFractionCorr <- function(crc_age, crc_eth, title, x_name, y_name)
  {
   
    ## correlate CRC ratios of young and old rates
    ## between two CRC populations

    ## Filter NaN, Inf values of young_rate/old_rate
    crc_grp1 <- filter(crc_eth, (!is.na(fraction)) & (!is.infinite(fraction)) & (fraction > 0 ))
    crc_grp2 <- filter(crc_age, (!is.na(fraction)) & (!is.infinite(fraction)) & (fraction > 0 ))


    ## join(x=crc_eth y=crc_age)
    combn <- inner_join(crc_grp1, crc_grp2, by=c('County','FIPS')) %>%
      select(FIPS, County,
             young_crc_eth= young.x, old_crc_eth= old.x,
             young_crc_ref=young.y, old_crc_ref= old.y,
             young_pop_eth=young_pop.x, old_pop_eth=old_pop.x,
             young_pop_ref=young_pop.y, old_pop_ref=old_pop.y,
             ## young_rate, old_rate,
             fraction_ethnic= fraction.x, fraction_ref= fraction.y) %>%
       mutate(eth_young_rate = young_crc_eth/young_pop_eth, ref_young_rate = young_crc_ref/young_pop_ref)

    ## compute RMS
    dis <- combn$fraction_ethnic - combn$fraction_ref
    rms <- sqrt(mean(dis^2))

    ## Wilcoxon singed-rank test
    test <- wilcox.test(combn$eth_young_rate, combn$ref_young_rate, paired=TRUE, alt='greater')
    
    print(dim(combn))
    p <- ggplot(combn, aes(x=fraction_ref, y=fraction_ethnic)) + geom_point()
    ## p <- p + geom_smooth(method=lm, fullrange=TRUE)
    p <- p + geom_abline(intercept = 0,slope=1, colour = "red")
    ## p <- p + annotate("text", label=paste0("rms=", round(rms,4)), x=min(combn$fraction_ref)+0.1,y=max(combn$fraction_ethnic)-0.1, size=3.5)
    p <- p + annotate("text", label=paste0("p-val=", format(test$p.value, scientific=T, digits=3)),
                      x=min(combn$fraction_ref)+0.1,y=max(combn$fraction_ethnic)-0.1, size=3.5)
    p <- p + theme(aspect.ratio=1, axis.title.x=element_text(size=9, vjust=-0.6), axis.title.y=element_text(size=9))
    p <- p + labs(title=title,
                  x=paste(x_name, "E-CRC"),
                  y=paste(y_name, "E-CRC"))
    
    return(p)
  }

Quantiles2Value <- function(...)
  {
    values <- c(...)
    
    quant <-  quantile(values, prob=seq(0,1,0.01))
    reduced_quant <-  rle(quant)[['values']]
    names(reduced_quant) <- as.numeric(sub('%', '',names(reduced_quant)))
    return(reduced_quant)
    
  }

RateECDF <- function(...)
  {
    ## Generate ecd function from arbitray number of vector values
    ## returns ecdf objects (i.e. function) that takes value as input and returns
    ## the quantile
    values <- c(...)

    ecdfFn <- ecdf(values)
       
    return(ecdfFn)
  }

#######################
## use the following package to
## design color palette
## library(colorspace)
## pal <- choose_palette()
#######################
plot2file <- TRUE
if(plot2file){
    pdf(paste0('../results/CRCAgeMaps_', Sys.Date(), '.pdf'), width=20, height=15)
}


######################
## Load CRC data for
## all and ethnic groups
######################
## yong_rate and old_rate are rates per age population not total county population
## therefore, no need to pop_census data
crc <- LoadCRCData("../data/CRC_county_clinical_2015-10-20.csv")
crc_black <- LoadCRCData("../data/CRC_county_Black_2015-10-20.csv")
crc_hispanic <- LoadCRCData("../data/CRC_county_Hispanic_2015-10-20.csv") 
crc_white <-  LoadCRCData("../data/CRC_county_White_2015-10-20.csv")

#######################
## Census data - obsolete
######################
pop_census <- LoadCensus("../data/us_census_2011.csv")
pop_census <- select(pop_census, FIPS, Geography, young_population, old_population, TOTAL)

## crc_age_adj <- left_join(crc, pop_census, by="FIPS") %>%
##   select(FIPS, County, young, old, young_population, old_population) %>%
##   mutate(young_rate = 1e5*(young/young_population), old_rate = 1e5*(old/old_population))

#############
## Age adjusted rates
#############
scale_factor <- 1e5
crc_age_adj <- mutate(crc, young_rate = scale_factor*(young/young_pop), old_rate = scale_factor*(old/old_pop))
crc_black <- mutate(crc_black, young_rate = scale_factor*(young/young_pop), old_rate = scale_factor*(old/old_pop))
crc_hispanic <- mutate(crc_hispanic, young_rate = scale_factor*(young/young_pop), old_rate = scale_factor*(old/old_pop))
crc_white <- mutate(crc_white, young_rate = scale_factor*(young/young_pop), old_rate = scale_factor*(old/old_pop))

require(scales)
ecdf_crc_young_rate_Fn <- RateECDF(crc_age_adj$young_rate, crc_black$young_rate, crc_hispanic$young_rate)
ecdf_crc_old_rate_Fn <- RateECDF(crc_age_adj$old_rate, crc_black$old_rate, crc_hispanic$old_rate)
q2v_young <- Quantiles2Value(crc_age_adj$young_rate, crc_black$young_rate, crc_hispanic$young_rate)
q2v_old <- Quantiles2Value(crc_age_adj$old_rate, crc_black$old_rate, crc_hispanic$old_rate)

#####################
## fill gradient color palettes
#####################
young_rate_max <-  max(max(crc_age_adj$young_rate), max(crc_black$young_rate), max(crc_hispanic$young_rate))
old_rate_max <-  max(max(crc_age_adj$old_rate), max(crc_black$old_rate), max(crc_hispanic$old_rate))

sfg_old <- scale_fill_gradient(low = "#FFFFFF", high = "#A20043", ##low = "#FEECEF", high = "#A20043",
                               trans=trans_new('old_trans', function(x) {ecdf_crc_old_rate_Fn(x) },
                                 function(x) {as.integer(q2v_old[findInterval(100*x, as.numeric(names(q2v_old)), all.inside=T)]) },
                                 breaks=extended_breaks(n=4) ),
                               limits=c(0,old_rate_max),
                            name=expression('Age adj. rates'))

sfg_young <- scale_fill_gradient(low = "#FFFFFF", high = "#A75F00", ## low = "#F7EFEA", high = "#A75F00",
                               trans=trans_new('young_trans', function(x) ecdf_crc_young_rate_Fn(x),
                                 function(x) {as.integer(q2v_young[findInterval(100*x, as.numeric(names(q2v_young)), all.inside=T)]) },
                                 breaks=extended_breaks(n=3) ),
                                 limits=c(0, young_rate_max),
                                 name=expression('Age adj. rates'))

######################
## plot entire US young
## and old rates and ethnic population
######################
populations_plots <- list(young=PlotMap(crc_age_adj, 'young_rate', sfg_young,ttl= 'Young rates - All'),
                          hispanic=PlotMap(left_join(crc_hispanic, pop_census, by='FIPS') %>%
                            mutate(Hispanic=(young_pop + old_pop)/TOTAL), 'Hispanic',
                            sfg=scale_fill_gradient(low = "#fff5eb", high = "#7f2704",name='fraction'),ttl= 'Young Hispanic pop'),
                          black=PlotMap(left_join(crc_black, pop_census, by='FIPS') %>%
                            mutate(Black=(young_pop + old_pop)/TOTAL), 'Black',
                            sfg=scale_fill_gradient(low = "#f7fbff", high = "#08306b",name='fraction'), ttl='Young Black pop'),
                          old = PlotMap(crc_age_adj, 'old_rate', sfg_old, ttl='Old rates - All'),
                          white=PlotMap(left_join(crc_white, pop_census, by='FIPS') %>%
                            mutate(White=(young_pop + old_pop)/TOTAL), 'White',
                            sfg=scale_fill_gradient(low = "#f7fcf5", high = "#00441b",name='fraction'), ttl='Old White Pop')
                          )

## multiplot(plotlist=populations_plots, cols=2)


######################
## plot entire US young
## and old rates and ethnic rates
######################
national_crc_rate_plots <- list(
                          hispanic_young=PlotMap(crc_hispanic, 'young_rate',sfg_young, ttl='Hispanic E-CRC rates'),
                          white_young=PlotMap(crc_white, 'young_rate', sfg_young, ttl='White E-CRC rates'),
                          black_young=PlotMap(crc_black, 'young_rate', sfg_young, ttl='Black E-CRC rates'),
                                
                          hispanic_old=PlotMap(crc_hispanic, 'old_rate',sfg_old, ttl='Hispanic CRC rates'),
                          white_old=PlotMap(crc_white, 'old_rate',sfg_old, ttl='White CRC rates'),
                          black_old=PlotMap(crc_black, 'old_rate', sfg_old, ttl='Black CRC rates')
                          )

multiplot(plotlist=national_crc_rate_plots, cols=2)

######################
## plot by state young
## and old rates
######################
## states <- c('California', 'Connecticut', 'Michigan',
##             'Georgia', 'Hawaii', 'Iowa', 'Kentucky',
##             'Louisiana', 'New Jersey', 'New Mexico', 'Washington', 'Utah')

states <- c('California', 'Connecticut', 'New Jersey')

young_state_plots <- lapply(states, function(stt) {PlotState(crc_white, stt, 'young_rate', sfg_young, 'White')})
old_state_plots <- lapply(states, function(stt) {PlotState(crc_white, stt, 'old_rate', sfg_old, 'White')})
black_young_state_plots <- lapply(states, function(stt) {PlotState(crc_black, stt, 'young_rate', sfg_young, 'Black')})
black_old_state_plots <- lapply(states, function(stt) {PlotState(crc_black, stt, 'old_rate', sfg_old, 'Black')})
hispanic_young_state_plots <- lapply(states, function(stt) {PlotState(crc_hispanic, stt, 'young_rate', sfg_young, 'Hispanic')})
hispanic_old_state_plots <- lapply(states, function(stt) {PlotState(crc_hispanic, stt, 'old_rate', sfg_old, 'Hispanic')})

sapply(1:length(states), function(i) multiplot(old_state_plots[[i]], black_old_state_plots[[i]], hispanic_old_state_plots[[i]],
                                               young_state_plots[[i]],black_young_state_plots[[i]],hispanic_young_state_plots[[i]],
                                               cols=2))

###################
## 1. Plot old vs. young rates
## 2. Plot old and young densities
###################
## plot old_rate vs. young_rate
gg_points <- ggplot(crc_age_adj, aes(x=old_rate, y=young_rate)) + geom_point(color='red')
## print(gg_points)

crc_age_adj_reduced <- crc_age_adj %>% select(FIPS, County, young_rate, old_rate) %>%
  tidyr::gather(group, rate, old_rate:young_rate, na.rm=TRUE)

gg_density <- ggplot(crc_age_adj_reduced, aes(x=rate))
gg_density <- gg_density + geom_histogram(aes(y=..density..,fill=group), binwidth=25, alpha=0.3)
gg_density <- gg_density + geom_density(aes(colour = group)) 
gg_density <- gg_density + scale_color_manual(values=c('red', 'orange')) + scale_fill_manual(values=c('gray10', 'gray60'))
gg_density <- gg_density + theme_bw()
## print(gg_density)


########################
## Correlations between
## young rate/ old rate
## by ethnic groups
########################
## Correlate age-corrected rates and ethnic groups
ethnic_plots <- list(
                     ## CRCFractionCorr(crc_age_adj, crc_hispanic, title="Hispanic", 'All','Hispanic'),
                     ## CRCFractionCorr(crc_age_adj, crc_black, title="Black", 'All', 'Black'),
                     ## CRCFractionCorr(crc_age_adj, crc_white, title="White", 'All', 'White'),
                     ## CRCFractionCorr(crc_age_adj, crc, title="All", 'All', 'All'),
                     CRCFractionCorr(crc_white, crc_hispanic, title="White vs. Hispanic", 'White', 'Hispanic'),
                     CRCFractionCorr(crc_white, crc_black, title="White vs. Black", 'White', 'Black')
                     ## CRCFractionCorr(crc_hispanic, crc_black, title="Hispanic vs. Black", 'Hispanic', 'Black')
                     )
    
multiplot(plotlist=ethnic_plots, cols=1)

##############
## close plot device
##############
if(plot2file){
    dev.off()
}
