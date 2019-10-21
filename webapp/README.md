# diffexprwebapp: A dash (by plotly) web application for viewing expresion profiles
#### Jonathan Samson

diffexprwebapp creates a server that allows users to explore the expression
profiles of particular genes.  It also allows for correlation analysis to find
genes that are similar to that particular gene.


## Installation
diffexprwebapp requires python3.7, miniconda (or equivalent), and git 
installed on the system.

### Installing Miniconda
The first step to installation is to install miniconda (if not already present).
The installer and installation instructions can be found 
[here.](https://docs.conda.io/en/latest/miniconda.html)

After installation, channels must be added so conda can find the appropriate
packages:
```
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Load the conda environment.
Use the included dash.yml files to create a conda environment:
```
conda env create -f dash.yml
```

### Enter the appropriate conda environment.
Once you've created the conda environment, enter it:
```
conda activate dash
```
To deactivate an active environment, use:
```
conda deactivate
```

### Set the host and port to run the application on
The last line of the file runs the server.  Select the host and port
to run the server on, which will determine who can access it and how.
As packaged, it will run on the server, making it available to anyone
with a wired connection at WUR, through port 54323.

## Running
The script can be run with the following command:
```
python diffexpressiontool.py
```

If run from the bioinformatics servers, and the host and port remain
unchanged, it can be accessed from the WUR network through a LAN connection
(must be plugged in, wireless will prevent it from working) by going to
bioinformatics.nl:54324

## Options Available
When running the application, the following options are availble for selection:

### File
The files for different conditions.  This allows the user to decide which set
of data they want to use.

### Loci
The gene loci of the genes they wish to view.  Any number of loci can be selected.

### Correlation
Whether they wish to see which loci correlate with the selected gene loci.
Note that correlation will only work if there is one loci selected.  The
correlation is based on the slider below.

### Correlation based on Absolute Value
Whether this correlation should be only positive (unchecked) or positive and
negative (checked).  The negative correlation will also show genes which are
regulated in the opposite pattern to the seleted loci.  This box only comes
into effect if the correlation box is also checked.
