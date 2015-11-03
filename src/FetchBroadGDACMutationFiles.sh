#!/bin/bash

BASEURL="http://gdac.broadinstitute.org/runs/analyses__2014_10_17/reports/cancer/"
OUTDIR="/zenodotus/dat01/betellab_store/dob2014/Shah/BROAData"

######################### 
## use `curl` command to download
## entire TCGA data
##########################
# for case in $(curl -s $BASEURL | grep href | sed 's/.*href="//' | sed 's/".*//' | grep '^[a-zA-Z].*' | sed 's/[/]//'); 
# do
#     echo $case
#     URL="$BASEURL/$case/MutSigNozzleReport2CV/sig_genes.txt"
    
#     ## check HTTP status
#     status=`curl -s -o /dev/null -I -w "%{http_code}" $URL`
#     echo $status
#     if test $status -ne 404
# 	then
# 	## fetch file
# 	curl $URL -o $OUTDIR/${case}_sig_genes.txt
#     fi

# done    

#########################
## Use firehose_get for specific data set
#########################
FIREHOSE_GET="/home/dob2014/lib/firehose_get/firehose_get"

cd $OUTDIR

## collect all data from COADREAD study (Colorectal Adenocarcinoma)
$FIREHOSE_GET -batch stddata latest COADREAD
# $FIREHOSE_GET -batch stddata latest COAD

## collect all analyses from COADREAD study
$FIREHOSE_GET -batch analyses latest COADREAD
# $FIREHOSE_GET -batch analyses latest COAD

## collect all 'mutation' analyses from all studies
## March 26. This download did not complete - wget error 8
# $FIREHOSE_GET -batch -o mut analyses latest
