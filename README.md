# exonCNV
A tool for calling Exon CNVs by NGS data.

## Install
It can be installed from github:

```
devtools::install_github('sanadamakomi/exonCNV')
library(exonCNV)
```

## Pipeline to call CNVs
'DExonCNV_v1.0.R' in exec/ is a batch tool to run exonCNV.

Before using you need install package:

```
installed.packages("GetoptLong")
```

Please use 'Rscript' to run batch tool. See more details by command:

```
Rscript DExonCNV_v1.0.R -help

# see mode

Rscript DExonCNV_v1.0.R -modeHelp depth
```

## Input file format

### annotation file

It can be download from [UCSC database](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz). Or you can use hg19_refGene.txt from annovar tool database.

### Index File

It's a TAB-delimited text file which has two columns without header.

|Gene_Symbol_in_annotation_file|Transcript_Name_in_annotation_file|
|-|-|
|DMD|NM_004009|
|SMN1|NM_000344|

### Update

2019.3.27 Debug callGenderByCov(): Use "SRY" gene to call gender.

2021.3.24 Debug depthOfRegion(): Calculate coverage by chromosome.

2021.6.3  Debug outputVcf(): Optimize output file.

2022.3.30  Debug filterExonCNV(): Filter cn = 0.

2022.7.27  Edit batchCall(): Add parameter geneName to call cnv in genes.
