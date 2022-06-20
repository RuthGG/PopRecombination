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
* 0-based

## genomicSuperDups_coords.bed

```
cd data/
awk  -v OFS="\t" '{print $2, $3, $4}' genomicSuperDups.txt > genomicSuperDups_coords.bed
awk -v OFS="\t" 'NR>1{print $8, $9, $10}' genomicSuperDups.txt >> genomicSuperDups_coords.bed

``` 

# Recombination rate data in hg19 from Spence and Song (SpenceSong_hg19_recMaps/)

BED files are 0-based, half.open intervals (python style)

Source / tutorial: https://medium.com/@paudelanjanchandra/download-google-drive-files-using-wget-3c2c025a8b99

```
mkdir tmp/
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv" -O hg19_maps.tar.gz && rm -rf /tmp/cookies.txt
tar -zxvf hg19_maps.tar.gz 
```
## Processed Spence maps (SpenceSong_hg19_recMaps_processed/)

Al chromosome bed files calculated for each population from Spence and Song

```
cd data/
mkdir SpenceSong_hg19_recMaps_processed/
for POP in $(ls SpenceSong_hg19_recMaps/); do 
head -n 1 SpenceSong_hg19_recMaps/${POP}/${POP}_recombination_map_hg19_chr_1.bed > SpenceSong_hg19_recMaps_processed/${POP}_recombination_map_hg19_allChr.bed
for FILE in $(ls SpenceSong_hg19_recMaps/${POP}/${POP}_recombination_map_hg19_chr_*.bed); do
tail -n+2 $FILE >> SpenceSong_hg19_recMaps_processed/${POP}_recombination_map_hg19_allChr.bed
done
done

```

# chrSizes.txt

* Manually generated from: https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
* Steps: generate simple headers separated with tabs, delete all single spaces, delete all commas, replace ^ by chr, delete chr in first column name

# Bherer map (Bherer_Refined_genetic_map_b37/)

These are 0-based indexed, half-open intervals (python style)

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
cd data/
mkdir Bherer_Refined_genetic_map_b37_processed/
head -n 1 Bherer_Refined_genetic_map_b37/sexavg_chr1.txt > Bherer_Refined_genetic_map_b37_processed/sexavg_allChr.txt
for FILE in $(ls Bherer_Refined_genetic_map_b37/sexavg_chr*); do 
tail -n+2 $FILE >> Bherer_Refined_genetic_map_b37_processed/sexavg_allChr.txt
done 
head -n 1 Bherer_Refined_genetic_map_b37/female_chr1.txt > Bherer_Refined_genetic_map_b37_processed/female_allChr.txt
for FILE in $(ls Bherer_Refined_genetic_map_b37/female_chr*); do 
tail -n+2 $FILE >> Bherer_Refined_genetic_map_b37_processed/female_allChr.txt
done 
```

To make them bed format equivalent to Spence maps:

```
head -n 1 SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed > Bherer_Refined_genetic_map_b37_processed/sexavg_allChr.bed
for FILE in $(ls Bherer_Refined_genetic_map_b37/sexavg_chr*); do 
awk -v OFS="\t" 'NR == FNR {arr[FNR] = $2;next} { print $1, $2 ,arr[FNR+1], $3/1000000}' $FILE $FILE | tail -n+2 |head -n-1 >> Bherer_Refined_genetic_map_b37_processed/sexavg_allChr.bed 
done
head -n 1 SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed > Bherer_Refined_genetic_map_b37_processed/female_allChr.bed
for FILE in $(ls Bherer_Refined_genetic_map_b37/female_chr*); do 
awk -v OFS="\t" 'NR == FNR {arr[FNR] = $2;next} { print $1, $2 ,arr[FNR+1], $3/1000000}' $FILE $FILE | tail -n+2 |head -n-1 >> Bherer_Refined_genetic_map_b37_processed/female_allChr.bed 
done


```