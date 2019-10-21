#! /usr/bin/env python3
"""
Author: Jonathan Samson

Uses kallisto quant to pseudoalign raw RNA-Seq counts to
a reference.

A kallisto index must be built before running the script.
This is often done using a cdna file.
"""

from pathlib import Path
import argparse
import re
import subprocess


def argument_parser():
    """ Parses the arguments of the function.

    Returns:
            parser.parse_args() -- The parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Runs kallisto quant " +
                                    "on the files in the input folder.",
                                    epilog="Kallisto index must be run " +
                                    "prior to running this program, and " +
                                    "the index must be stored in the " +
                                    "output folder.  The Araport11" +
                                    "_411_cdna.fasta file is typically " +
                                    "used for the index.")
    output = parser.add_argument("-o", "--output_folder",
                               metavar="/PATH/TO/OUTPUT/",
                               help="The directory to output the " +
                               "psudoalignment to.",
                               default="/mnt/scratch/samso008/" +
                               "Project103470/kallisto/")
    input_dir = parser.add_argument("-i", "--input_folder",
                               metavar="/PATH/TO/INPUT/",
                               help="The directory containing the raw read " +
                               "files.  These files should be in .fasta or " +
                               ".fasta.gz format.",
                               default="/mnt/scratch/samso008/" +
                               "Project103470/RawData/")
    samples = parser.add_argument("-s", "--sample_names",
                               metavar="/PATH/TO/SAMPLE/FILE",
                               help="A tab-separated file containing, the " +
                               "first column being the file name of the " +
                               "raw data, and the second being the " +
                               "condition of the sample.  The condition " +
                               "should be in the following format: " +
                               "Construct_S(ample)Number_Time_Treatment.",
                               default="./samplenames.csv")
    threads = parser.add_argument("-t", "--num_threads", metavar="#", type=int,
                               help="The number of threads to use when running " +
                               "kallisto quant.",
                               default=20)
    usage = parser.format_usage()
    parser.usage = usage
    for argument in [output, input_dir, samples, threads]:
        argument.metavar=""
    return parser.parse_args()

def generate_directory_file(input_directory, output_file):
    """ Creates a file containing the names of all the RNA-seq files.

    Keyword Arguments:
    input_directory -- (Path) The /PATH/TO/DIRECTORY/ that the
                       RNA-seq files are stored (in .fastq or
                       .fastq.gz format).
    output_file -- (Path) The /PATH/TO/FILE to output the
                   RNA-seq filenames to.
    """
    try:
        if not output_file.exists():
            cmd = f"ls {input_directory}*.fastq* > {output_file}"
            subprocess.check_call(cmd, shell=True)
    except IOError:
        print("Error in generate_directory_file; the file already exists. " +
              "\nContinuing with the previously-generated file")

def run_kallisto_quant(wrk_dir, raw_file, conversion, raws, threads):
    """ Runs kallisto quant on each raw read file.

    Keyword Arguments:
        wrk_dir -- (Path) The /PATH/TO/DIRECTORY/ containing the
                   index.
        raw_file -- (sting iterator) The /PATH/TO/FILE containing a
                    list of raw read file names.
        conversion -- {(string)read_name:(string)sample_name}
                      Dictionary containing the file name of the
                      raw reads as the keys, and human
                      readable file names as the values.
        raws -- (Path) The /PATH/TO/DIRECTORY with raw read files.
        threads -- (int) The number of threads to use.
    """
    for line in raw_file:
        reg = re.search(
            r"(\w*)_(\w{6}(-\w*)*)(-|\_)\w{8}-\w{8}\_(\w*)\_\w*\.fastq.*",
            line)
        run_name = reg.group(2)
        out_folder = conversion[run_name]
        try:
            folder = wrk_dir / out_folder
            if not folder.exists():
                # Standard deviation (-s 30) is from community forum
                # posts.  Length is from "Protocol for use with NEBNext
                # Ultra Directional RNA Library Prep Kit for Illumina".
                # Company claims to use this protocol, and doesn't site
                # any deviations from typical (though the protocol had
                # no addendums for different fragment lengths).
                cmd = (f"kallisto quant -t {threads} -i {wrk_dir}" +
                       "/411index --single " +
                       "-l 200 -s 30 -b 100 " +
                       f"-o {wrk_dir}/{out_folder}" +
                       f" {raws}/*{run_name}*")
                subprocess.check_call(cmd, shell=True)
        except IOError:
            print("Error in run_kallisto_quant; the output file for " +
                  f"{out_folder} already exists. Continuing with " +
                  "previously generated file.")

def main():
    """ Runs when the module is called. """
    # Define variables
    arguments = argument_parser()
    kallisto_dir = Path(arguments.output_folder)
    raw_dir = Path(arguments.input_folder)
    conv_file = Path(arguments.sample_names)
    quant_threads = arguments.num_threads
    raw_files = Path("./seq_files")  # A file to hold raw read file names
    sample_names = {}  # A dictionary of filenames: samplenames

    # Create a file holding the names of the raw file reads
    generate_directory_file(raw_dir, raw_files)

    # Fills the dictionary of run_names, sample_names
    with open(conv_file) as converter:
        for line in converter:
            reg = re.search("(.*)\\t(.*)\\n", line)
            key = reg.group(1)
            value = reg.group(2)
            sample_names[key] = value

    # Runs kallisto quant on each of the raw reads files
    with open(raw_files) as raw_path:
        run_kallisto_quant(kallisto_dir, raw_path, sample_names, raw_dir,
                           quant_threads)


if __name__ == "__main__":
    main()
