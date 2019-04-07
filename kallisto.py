#! /usr/bin/env python
"""
Author: Jonathan Samson

A script to take RNA-seq data and run it through kallisto and sleuth to give
data to be interpreted.

Keyword Arguments:
argv[1] -- the kallisto directory
argv[2] -- The directory that contains the RNA-seq raw data
argv[3] -- A .csv file with run_names, sample_names

A kallisto index must be built before running the script.
"""

### Import functions
from sys import argv
import re
import subprocess
import os.path

### Definte Functions
def gen_file_list(input_directory, output_file):
    """ Creates a file with the names of all the RNA-seq files.

    Keyword Arguments:
    input_directory -- The directory that the RNA-seq files are stored (in
        .fastq or .fastq.gz format)
    """
    try:
        if not os.path.exists(output_file):
            cmd = 'ls {}*.fastq* > {}'.format(input_directory, output_file)
            subprocess.check_call(cmd, shell=True)
    except IOError:
        print("Error in gen_file_list; the file already exists. " +\
              "\nContinuing with the previously-generated file")

def run_kallisto_quant(wrk_dir, raw_file, conversion, raws):
    """ Runs kallisto quant on each raw read file.

    Keyword Arguments:
        wrk_dir -- The kallisto directory to output to, containing the index
        raw_file -- The file containing a list of raw read file names
        conversion -- Dictionary containing read_name: sample_name
        raws -- Directory with raw read files
    """
    for line in raw_file:
        reg = re.search(r"(\w*)_(\w{6}(-\w*)*)(-|\_)\w{8}-\w{8}\_(\w*)\_\w*\.fastq.*",\
                               line)
        run_name = reg.group(2)
        out_folder = conversion[run_name]
        try:
            if not os.path.exists(wrk_dir + out_folder):
                """ Standard deviation (-s 30) is from community forum posts.
                Length is from "Protocol for use with NEBNext Ultra
                Directional RNA Library Prep Kit for Illumina".  Company claims to use
                this protocol, and doesn't site any deviations from typical (though the
                protocol has addendums for different fragment lengths).
                """
                cmd = ('kallisto quant -t 20 -i ' + wrk_dir + '/411index --single '+\
                    "-l 200 -s 30 -b 100 -o " + wrk_dir + out_folder +\
                    " " + raws + "*" + run_name + "*")
                subprocess.check_call(cmd, shell=True)
        except IOError:
            print("Error in run_kallisto_quant; the output file for " +\
                  "{} already exists. Continuing with".format(out_folder) +\
                  "previously generated file.")

def main():
    """This function is run when the file is called.
    """
    # Define variables
    raw_files = "./seq_files" # A file to hold raw read file names
    sample_names = {} # A dictionary of filenames: samplenames
    kallisto_dir = argv[1] # The directory with the kallisto index
    raw_dir = argv[2] # The directory with the raw read files
    conv_file = argv[3] # The file with run_name,sample_name

    # Create a file holding the names of the raw file reads
    gen_file_list(raw_dir, raw_files)

    # Fills the dictionary of run_names, sample_names
    with open(conv_file) as converter:
        for line in converter:
            reg = re.search("(.*)\\t(.*)\\n", line)
            key = reg.group(1)
            value = reg.group(2)
            sample_names[key] = value

    # Runs kallisto quant on each of the raw reads files
    with open(raw_files) as raw_path:
        run_kallisto_quant(kallisto_dir, raw_path, sample_names, raw_dir)

# Main executable
if __name__ == "__main__":
    main()
