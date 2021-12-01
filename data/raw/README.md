gap.txt
from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
table explannation: all tables - gap http://genome.ucsc.edu/cgi-bin/hgTables
0-based

cytoband.txt
from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
all tables - cytoband http://genome.ucsc.edu/cgi-bin/hgTables
0-based

inversionsAnnotation - Mario sent it in a mail, transformed into CSV in use/

repeatMasker.txt
from UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1227000939_PXdA8FuNSAiCxfofGBuZawDKW5tO&clade=mammal&org=Human&db=hg19&hgta_group=rep&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr17%3A40%2C565%2C779-42%2C502%2C878&hgta_outputType=primaryTable&hgta_outFileName=genomicSuperDups.txt
table explannation: group Repeats, track RepeatMasker, table rmsk ("describe table schema")

genomicSuperDups.txt
from UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1227000939_PXdA8FuNSAiCxfofGBuZawDKW5tO&clade=mammal&org=Human&db=hg19&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=0&hgta_regionType=genome&position=chr17%3A40%2C565%2C779-42%2C502%2C878&hgta_outputType=primaryTable&hgta_outFileName=repeatMasker.txt
table explannation: group Repeats, track SegmentalDups, table gneomicSuperDups ("describe table schema")

	# # Dowload recombination rate data - [RECRATE.BED + GENSIZES -- but multichromosome] 
	# # *****************************************************************************

	# 	# Check source/tutorial before reusing this comand!!!!!
	# 	# Source / tutorial: https://medium.com/@paudelanjanchandra/download-google-drive-files-using-wget-3c2c025a8b99

	# 	wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17KWNaJQJuldfbL9zljFpqj5oPfUiJ0Nv" -O ${PARENT}/data/hg19_maps.tar.gz && rm -rf /tmp/cookies.txt
	# 	tar -zxvf ${PARENT}/data/hg19_maps.tar.gz -C ${PARENT}/data/
