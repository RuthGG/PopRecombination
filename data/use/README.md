# inversionsAnnotation 

* Mario sent it in a mail, transformed into CSV in use/. Deleted HsInv0030 and HsInv1122 (deletions)

# BhererGeneticDistance.txt

```
# while in the main dir
for i in $(ls data/raw/Refined_genetic_map_b37/sexavg_chr*); do tail -n1 $i; done > data/use/BhererGeneticDistance.txt
```

# BhererAllChroms.txt

```
# while in the main dir
for i in $(ls data/raw/Refined_genetic_map_b37/sexavg_chr*); do tail -n+2 $i >> data/use/BhererAllChroms.txt; done 
```
