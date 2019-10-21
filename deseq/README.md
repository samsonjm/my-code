# DESeq2.R: An R script to run DESeq2 on mapped read data.
#### Jonathan Samson

DESeq2.R runs analysis on the HISAT2 pipeline output.  It is used to create
heatmaps, expression plots, and other visual representations of the data.

## Installation
DESeq2.R requires R 3.4.4 and git installed on the system.

### Data files.
The "data" folder contains files that DESeq2.R must read in order to obtain 
the count data.  Additionally, the library paths and working directory
must be set.  All of these paths are at the top of the file, lines 14-20.

## Running
The first time DESeq2.R is used, the appropriate packages must be installed.
These are located in the section that goes from line 23 to line 38.  After
they are installed, this section should be commented out to save time.
These packages are xml2, igraph, DESeq2, graph, RCy3, EBSeq, biomaRt, topGO,
blockmodeling, mnormt, multtest, Hmisc, and psych.

The functions and data files should then be read into the workspace.  This 
allows R to call the functions.  If a function is read into memory and then
changed, it must be read into memory again to update the changes.

Finally, below the functions are the calls.  These are separated based on the
data set that they use.  The data sets are defined as follows:

Total: 401 + 411 + Col0 Shoot + PLT + Han

Base: 401 + 411 + Col0 Shoot

411: 411

401: 401

401shoot: 401 + Col0 Shoot

411shoot: 411 + Col0 Shoot

Exp: 401 + 411

The total data set was not typically used, as the PLT and Han samples were not
a part of this experiment, but was included for completeness.  Its section has
been commented out to reflect this.

After the separate calls for the data sets, there is a section for comparing
this pipeline with the kallisto -> sleuth pipeline.  The sleuth files will
not be provided here, but can be generated using that pipeline.  This section
performs pearson correlation between the two pipelines to see whether they 
provide similar results.
