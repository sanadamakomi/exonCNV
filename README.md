# exonCNV
A tool for calling Exon CNVs(Copy number variations)

## Install
It can be installed from github:

```
devtools::install_github('sanadamakomi/exonCNV')
library(exonCNV)
```

## Pipeline to call CNVs:
'DExonCNV_v1.0.R' in exec/ is a batch tool to run exonCNV.

Before using, need install package:

```
installed.packages("GetoptLong")
```

Please use 'Rscript' to run. See more details by command:

```
Rscript DExonCNV_v1.0.R -help

# see mode

Rscript DExonCNV_v1.0.R -modeHelp depth
```
