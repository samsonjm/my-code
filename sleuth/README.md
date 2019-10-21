# sleuth.R: An R script to run sleuth on pseudomapped read data.
#### Jonathan Samson

sleuth.R runs analysis on the kallisto output.  It is used to create
visual representations of the data.

## Installation
sleuth.R requires R 3.4.4 and git installed on the system.

### Data files.
The "data" folder contains the design.txt file that sleuth.R must read 
in order to separate out the kallisto files based on condition. 

## Running
The first time sleuth.R is used, the paths must be set, and 
the appropriate packages must be installed.  
The paths are in lines 5-7, and the package installation from  line 10 to 
line 25.  Important is to change the path to reflect the location of the
design.txt file, located in the data folder.  After the packages
 are installed, their section should be commented out to save time.
These packages are rhdf5, Rcpp, devtools, scales, tidyselect, htmltools,
plyr, pkgconfig, reshape2, ggfortify, stringi, tidyverse, purrr,
gplots, and sleuth.

The functions and data files should then be read into the workspace.  This 
allows R to call the functions.  If a function is read into memory and then
changed, it must be read into memory again to update the changes.

Below the functions are the calls.  These are separated based on the
data set that they use.  The data sets are defined as follows:

Total: 401 + 411 + Col0 Shoot + PLT + Han

Base: 401 + 411 + Col0 Shoot

411: 411

401: 401

401sh: 401 + Col0 Shoot

411sh: 411 + Col0 Shoot

Syn: 401 + 411

The total data set was not typically used, as the PLT and Han samples were not
a part of this experiment, but was included for completeness.

The bottom of the calls section has a couple calls that could be used for
comparing specific conditions (For example, the 411 4h Mock sample with the
411 4h Dex sample, as written currently).  This section was not yet used
to generate extensive data.

At the end, there is a section that uses the sleuth data to plot a PCA
for conditions, plot bootstraps, and launch a live sleuth interface. 
These were all set up manually and not automated, which is why they are
in a separate function.  The scripts provided here may not have extensive use
for future research, unless new data sets are generated.
