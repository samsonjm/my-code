# runkallisto: A script to run kallisto pseudomapping on RNA-Seq data.
#### Jonathan Samson

runkallisto allows for the input of multiple RNA-seq files and runs a
short pipeline to allow pseudomapping of all files.  It also takes a
sample name conversion file in to give human-readable names to the
mapped files.

## Installation
runkallisto requires python3.7, miniconda (or equivalent), kallisto, and git 
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

### Optional: Create a new conda environment.
If you don't already have a conda environment you wish to work in, a new one
should be created:
```
conda create -n envname python
```

### Enter the appropriate conda environment.
Once you've selected a conda environment, enter it:
```
conda activate envname
```
To deactivate an active environment, use:
```
conda deactivate
```

### Install kallisto in the environment.
Kallisto 0.44 must now be installed in the conda environment.
```
conda install kallisto=0.44
```
Say yes to installing the additional packages.

### Clone the git repository into a new folder.
Finally, the git repository needs to be forked, then cloned into a new folder.
Cloning can be done via SSH:
```
git clone git@git.wur.nl:samso008/runkallisto
```
or HTTPS:
```
git clone https://git.wur.nl/samso008/runkallisto
```

## Running
Prior to running the script, a kallisto index must be built.  This can be done
using kallisto index:
```
Usage: kallisto index [arguments] FASTA-files

Required argument:
-i, --index=STRING          Filename for the kallisto index to be constructed 

Optional argument:
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
    --make-unique           Replace repeated target names with unique names
```

Additionally, a tab-separated file must be created that has the raw RNA-seq
file names in the first column, and output file names in the second.  It is
recommended that the output file names be informative, giving the 
conditions of the sample for further analysis.

The program runs on a single command with the following input:

```
Usage: kallisto.py [-h] [-o /PATH/TO/OUTPUT/] [-i /PATH/TO/INPUT/]
                   [-s /PATH/TO/SAMPLE/FILE] [-t #]

optional arguments:
  -h, --help            show this help message and exit
  -o , --output_folder 
                        The directory to output the psudoalignment to.
  -i , --input_folder   The directory containing the raw read files. These
                        files should be in .fasta or .fasta.gz format.
  -s , --sample_names   A tab-separated file containing, the first column
                        being the file name of the raw data, and the second
                        being the condition of the sample. The condition
                        should be in the following format:
                        Construct_S(ample)Number_Time_Treatment.
  -t , --num_threads    The number of threads to use when running kallisto
                        quant.
```

The process follows these steps:

### Step 1: Gather the RNA-Seq file names.
runkallisto goes into the input folder and parses the file names there.
The file names are then output into a file.

### Step 2: Create a dictionary of file names and sample names.
The user-created sample name file is read so the output files are easier to 
sort and analyze.  

### Step 3: Kallisto is run.
The script gathers all input files with the same condition (based on 
the input file names) and runs kallisto on them to pseudomap the reads to 
the user-generated index.

#### Optional arguments:
Optionally, you may provide the number of threads for kallisto to use. (option
-t, default is 20).

### User friendliness:
If the pipeline is going to be run multiple times using the same input folder,
output folder, or samplename file, the default arguments in the argument_parser
function (the very first function of the script) can be changed.  This prevents
the user from having to input the arguments each time the script is run.
