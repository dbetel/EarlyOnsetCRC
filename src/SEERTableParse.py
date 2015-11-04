import sys
if sys.version_info < (2, 7):
    print "must use python 2.7 or greater"
    sys.exit(0)

import argparse
import re
import collections
import sys

## difflib for differences between strings
import difflib 
sys.path.append('/home/dob2014/lib/xlrd-0.9.2')
import xlrd

def OutputData(fips,cnty, lctn, cc, ppltn, ages, ofl):
    '''
    fips: FIPS geographical code
    cnty: county name
    lctn: anatomical location
    cc: list of counts by ages
    ppltn: list of population by ages
    ages: list of age ranges
    ofl: output file
    '''
    young_count = sum(cc[0:6])
    old_count = sum(cc[6:])
    young_rate = (100000 * young_count) / sum(ppltn[0:6])
    old_rate = (100000 * old_count) / sum(ppltn[6:])

    print(fips, cnty, lctn, young_count, old_count, young_rate, old_rate)
    ofl.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n' %(fips, cnty, lctn, young_count, old_count,
                                                       young_rate, old_rate, sum(ppltn[0:6]), sum(ppltn[6:])))
        


def Array2Dictionary(aa):
    dd = collections.OrderedDict()
    for i in range(len(aa)):
        if not aa[i]:
            continue
            
        if aa[i] in dd.keys():
            dd[ aa[i] ].append(i)
        else:
            dd[aa[i]] = [i]

    return dd
              

def GetCellsData(rows, column_number, datasheet):
    '''
    get values from rows in column_number in datasheet
    '''
    column = datasheet.col(column_number)
    values = map(column.__getitem__, rows)
    values = [item.value for item in values]
    return values

        
#########
## Main
#########
'''
usage:  python SEERTableParse.py ../data/frequency3.xlsx 1 ../data/CRC_County_clinical_2015-10-07.csv

input : inputfile sheetnumber outputfile
'''
## Command-line arguments
parser = argparse.ArgumentParser(description='Parse SEER data tables')
parser.add_argument('i', nargs=1, help='Excel file')
parser.add_argument('s', nargs=1, help='sheet number')
parser.add_argument('o', nargs=1, help='<sample>.csv')

try:
    args=parser.parse_args()
except IOError:
    parser.error()


book = xlrd.open_workbook(args.i[0])
out_file = open(args.o[0], 'w')

## write header
out_file.write("## " + " ".join(sys.argv) + "\n")
out_file.write('FIPS\tCounty\tAnatomical location\tCounts Young\tCounts Old\tRates Young\tOld Rates\tYoung Population\tOld Population\n')

## read in Excel sheet
sh = book.sheet_by_index(int(args.s[0])-1)

## 1. read in counties and make dictonary of column indicies
counties = sh.row_values(1)
counties_dict = Array2Dictionary(counties)

## 2. read age column 0 and dictionary of row indicies
ages = sh.col_values(0)
ages_dict = Array2Dictionary(ages)

## 3. read anatomical location and make dictionary of row indicies
locations = sh.col_values(1)
locations_dict = Array2Dictionary(locations)
# remove anatomical location not starting with 'C[0-9]+'
reObj = re.compile('^C[0-9]+')
for key in locations_dict.keys():
    if not reObj.match(key):
        del locations_dict[key]

## 4. For each county (skipping Unknown)
reUnknwn = re.compile("Unknown")
for county in counties_dict.keys():
    if reUnknwn.search(county):
        continue
    else:
        if re.search("\([0-9]+", county):
            fips = str(county[county.index("(") + 1:county.rindex(")")])
        else:
            fips='NA'
        ## for each unique location
        for location in locations_dict.keys():
            counts = GetCellsData(locations_dict[location], counties_dict[county][1], sh)
            population = GetCellsData(locations_dict[location], counties_dict[county][2], sh)
            age_values = [str(ages[i]) for i in locations_dict[location]]

            ## if population is 0 for either old or young group skip this county.
            if sum(population[0:6]) == 0 or sum(population[6:]) == 0:
                continue
            else:       
                OutputData(fips, str(county), str(location), counts, population, age_values, out_file)
        
out_file.close()
