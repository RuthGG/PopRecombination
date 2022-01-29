# gap.txt

* from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
* table explannation: all tables - gap http://genome.ucsc.edu/cgi-bin/hgTables
* 0-based

# cytoband.txt

* from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
* all tables - cytoband http://genome.ucsc.edu/cgi-bin/hgTables
* 0-based

# inversionsAnnotation 

* Mario sent it in a mail, also available in CSV format. CSV version has Deleted HsInv0030 and HsInv1122 (deletions)  

# repeatMasker.txt

* from UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1227000939_PXdA8FuNSAiCxfofGBuZawDKW5tO&clade=mammal&org=Human&db=hg19&hgta_group=rep&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr17%3A40%2C565%2C779-42%2C502%2C878&hgta_outputType=primaryTable&hgta_outFileName=genomicSuperDups.txt
* table explannation: group Repeats, track RepeatMasker, table rmsk ("describe table schema")

# genomicSuperDups.txt

* from UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1227000939_PXdA8FuNSAiCxfofGBuZawDKW5tO&clade=mammal&org=Human&db=hg19&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=0&hgta_regionType=genome&position=chr17%3A40%2C565%2C779-42%2C502%2C878&hgta_outputType=primaryTable&hgta_outFileName=repeatMasker.txt
* table explannation: group Repeats, track SegmentalDups, table gneomicSuperDups ("describe table schema")

# Recombination rate data in hg19 from Spence and Song (SpenceSong_hg19_recMaps/)

Source / tutorial: https://medium.com/@paudelanjanchandra/download-google-drive-files-using-wget-3c2c025a8b99

```
mkdir tmp/
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv" -O hg19_maps.tar.gz && rm -rf /tmp/cookies.txt
tar -zxvf hg19_maps.tar.gz 
```
## Processed Spence maps (SpenceSong_hg19_recMaps_processed/)

### CEU_SRR.bed

SRRs calculated for CEU population from Spence and Song, just to make a test about how extreme is SRR>10

# chrSizes.txt

* Manually generated from: https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
* Steps: generate simple headers separated with tabs, delete all single spaces, delete all commas, replace ^ by chr, delete chr in first column name

# Bherer map (Bherer_Refined_genetic_map_b37/)

* manually downloaded from https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/blob/master/Refined_genetic_map_b37.tar.gz
* manually extracted

## Processed Bherer maps (Bherer_Refined_genetic_map_b37_processed/)

### BhererGeneticDistance.txt

```
# while in the main dir
for i in $(ls data/raw/Refined_genetic_map_b37/sexavg_chr*); do tail -n1 $i; done > data/use/BhererGeneticDistance.txt
```

### BhererAllChroms.txt

```
# while in the main dir
for i in $(ls data/raw/Refined_genetic_map_b37/sexavg_chr*); do tail -n+2 $i >> data/use/BhererAllChroms.txt; done 
```
