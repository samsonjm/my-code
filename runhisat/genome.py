#! /usr/bin/env pyhon3
"""
Author: Jonathan Samson

Runs the HISAT2 -> StringTie -> StringTie merge -> Cuffdiff ->
PrepDE.py -> DESEQ pipeline.
"""

# Import statements
import subprocess
import re
import argparse
from pathlib import Path

def argparser():
    """ Parses the arguments given when the script is called.

    Returns:
        arguments -- (parsed ArgumentParser) The parsed arguments.
    """
    parser = argparse.ArgumentParser(prog="genome.py",
                                     description="Runs a HISAT2 pipeline "+
                                     "on the raw read files to prepare " +
                                     "the counts for use with DESeq.",
                                     epilog="The pipeline runs HISAT2, " +
                                     "StringTie, Cuffdiff, and PrepDE.  " +
                                     "The output file can then be utilized " +
                                     "by the DESeq2.R script to produce " +
                                     "graphs and other statistical data.")
    parser.add_argument("-w", "--wrk_dir", metavar="/PATH/TO/DIRECTORY/",
                        default="/mnt/scratch/samso008/Project103470/",
                        help="The /PATH/TO/DIRECTORY/ to work from.")
    parser.add_argument("-r", "--raw_folder", metavar="RawFolder/",
                        default="RawData/",
                        help="The folder within the working directory " +
                        "containing the raw reads files.")
    parser.add_argument("-s", "--sample_name_converter",
                        metavar="/PATH/TO/FILENAME",
                        default="./samplenames.csv",
                        help="The /PATH/TO/SAMPLE_NAME_CONVERTER " +
                        "containing the conversion " +
                        "from the raw read file name to the human readable " +
                        "file name.")
    parser.add_argument("-R", "--raw_file_names", metavar="PATH/TO/FILENAME",
                        default="./seq_files",
                        help="The /path/to/seq_filename that holds the raw " +
                        "read file names.")
    parser.add_argument("-o", "--out_folder", metavar="/PATH/TO/DIRECTORY",
                        default="/mnt/scratch/samso008/Project103470/hisat2/",
                        help="The PATH/TO/OUTPUT_FOLDER/ for the counts.")
    parser.add_argument("-a", "--annotation_file", metavar="/PATH/TO/FILENAME",
                        default="/mnt/scratch/samso008/Project103470/" +
                        "reference/Araport11_GFF3_genes_" +
                        "transposons_synthetic.201606.gff",
                        help="The /PATH/TO/ANNOTATION.GFF annotation file. " +
                        " The file can be in either .gff or .gtf format.")
    parser.add_argument("-m", "--merged_file", metavar="/PATH/TO/FILENAME",
                        default="/mnt/scratch/samso008/Project103470/" +
                        "hisat2/merged.gtf",
                        help="The /PATH/TO/MERGED.GTF file that gets ouput " +
                        "by the stringtie merge call.")
    arguments = parser.parse_args()
    return arguments

def separate_tuple(tup, wd):
    """ Returns a string of comma separated filenames with extension.

    Keyword Arguments:
    tup -- (tuple) The tuple that is to be separated.
    wd -- (Path) The working directory.

    Returns:
    separated -- (string) The comma separated file names with .gtf
                 extension.
    """
    separated = ""
    for element in tup:
        if separated == "":
            separated += str(element) + ".gtf"
        else:
            separated += " " + str(element) + ".gtf"
    return separated

def gen_file_list(input_directory, output_file):
    """ Creates a file with the names of all the RNA-seq files.

    Keyword Arguments:
    input_directory -- (Path) The directory that the RNA-seq files are stored (in
        .fastq or .fastq.gz format).
    output_file -- (Path) The name of the output file to be created.

    The output file is line-delimited.
    """
    if not output_file.exists():
        cmd = f"ls {input_directory}/*.fastq* > {output_file}"
        subprocess.check_call(cmd, shell=True)

def call_hisat2index(genome, indexname, wrk_dir):
    """ Calls hisat2-build to build an index from the listed genome.

    Keywork Arguments:
    genome -- (Path) .fasta file containing the genome to be indexed
              Araport_genome_synthetic.fa is usually used.
    indexname -- (string) the desired name of the index.
    wrk_dir -- (Path) The directory containing the working folders.
    """
    if not wrk_dir.joinpath("hisat2", indexname).exists():
        cmd = ("hisat2-build " + str(wrk_dir) + "reference/" + str(genome) +
               " " + str(wrk_dir) + "hisat2/" + indexname)
        subprocess.check_call(cmd, shell=True)

def call_hisat2(file_raws, conversion, out, folderraw):
    """ Calls hisat2, returns list of all .sam files.

    Keyword Arguments:
    file_raws -- (open file) File containing a list of raw read file
                 names.
    conversion -- {read_name:sample_name} What to convert read names
                  into.
    out -- (Path) Folder containing the HISAT2 indexes and sample
           outputs.
    folderraw -- (Path) Folder containing the raw reads.

    Returns:
    sam_list -- [Path] a list of all .sam files generated.
    """
    sam_list = []
    list_of_lines = []  # The list of lines in the raw file
    for line in file_raws:
        list_of_lines.append(line)
    file_raws.seek(0)
    for line in file_raws:
        current_raw = "" # a comma-separated list of the replicate reads
        reg = \
            re.search(r"((\w*)_(\w{6}(-\w*)*)(-|\_)\w{8}-\w{8}\_(\w*)\_\w*\.fastq.*)\n",
                  line)
        run_name = reg.group(3)
        out_folder = conversion[run_name]
        sam_file = out.joinpath(out_folder).with_suffix(".sam")
        if not sam_file.exists():
            current_raw += "/" + reg.group(1)
            for additional in list_of_lines:
                newreg = \
                    re.search(r"((\w*)_(\w{6}(-\w*)*)(-|\_)\w{8}-\w{8}\_(\w*)\_\w*\.fastq.*)\n*",
                                   additional)
                if newreg.group(3) == run_name and additional != line:
                    current_raw += "," + str(folderraw) + "/" +  newreg.group(1)
            # Arguments for hisat2:
            #   -x: path/to/indexfile
            #   -q: input files are FASTQ
            #   --pen-noncansplice: default 12
            #   --dta: alignments tailored for transcript assemblers
            #   -p: number of threads to launch
            #   -S: file for .sam output
            #   --no-temp-splicesite: prevents use of splice sites found by
            #                           earlier reads
            #   -U: comma separated list of raw read files
            cmd = ("hisat2 -x /mnt/scratch/samso008/Project103470/hisat2/411index"+
                " -q --dta -p 20 " +
                "--no-temp-splicesite -U " +
                str(folderraw) + current_raw + " -S " + str(sam_file))
            subprocess.check_call(cmd, shell=True)
        sam_list.append(sam_file.with_suffix(""))
    return (sam_list)

def call_samtoolssort(sam_files, out):
    """ Calls samtools which will sort and convert HISAT2 output to bam.

    Keyword Arguments:
    sam_files -- [Path]  a list of the .sam files (without extension)
                 generated by HISAT2.
    out -- (Path) the folder that has the .sam files, .bam files will
           output here.
    """
    for sam in sam_files:
        if not out.joinpath(sam.with_suffix(".bam")).exists():
            cmd = ("samtools sort -@ 20 -o " +
                   str(out.joinpath(sam.with_suffix(".bam "))) +
                   str(out.joinpath(sam.with_suffix(".sam"))))
            subprocess.check_call(cmd, shell=True)

def call_stringtie(bam_files, workingdir):
    """ Calls StringTie to assemble and quantify the transcripts.

    Keyword Arguments:
    bam_files -- [Path] a list of .bam files (without extension) from
                 samtools.
    workingdir -- (Path) The working directory.
    """
    for name in bam_files:
        try:
            if not Path(name).with_suffix(".gtf").exists():
                cmd = ("stringtie -p 20 -e -G " +
                       str(workingdir) + "/reference/" +
                    "Araport11_GFF3_genes_" +
                    "transposons_synthetic.201606.gff -B -o " + str(name) +
                    ".gtf -l " + str(name) + " " + str(name) + ".bam")
                subprocess.check_call(cmd, shell=True)
        except IOError:
            print("The .gtf file has already been created (stringtie)")

def stringtie_merge(guide, merged, input_list):
    """calls StringTie --merge to merge .gtf files.

    Keyword Arguments:
    guide -- (Path) the reference annotation in .GTF/.GFF3 format.
    merged -- (Path) the output file for the merged transcripts (.GTF).
    input_list -- [string] a list of .GTF/GFF files to merge.
    """
    try:
        if not merged.exists():
            cmd = "stringtie --merge -G " + str(guide) + " -o " +\
                  str(merged) + " " + input_list
            subprocess.check_call(cmd, shell=True)
    except IOError:
        print("The merged .gtf file (from stringtie --merge) already exists")

def group_conditions(file_names):
    """ Sorts the files into conditions, then returns a dictionary.

    Keyword Arguments:
    file_names -- [strings] list of the names of the .sam files,
                  without extension.

    Returns:
        grouped -- {condition(string):filename, filename(string)} the
                   dictionary of filenames grouped by condition.
    """
    grouped = {}
    condition = ""
    for name in file_names:
        # Reg ex to find condition
        reg = re.search(r".*\/(.*)-.*S.*\_(.*\_[a-zA-Z]+)", str(name))
        condition = reg.group(1) + reg.group(2)
        # Add name to dictionary with condition as key
        if condition in grouped:
            grouped[condition] = grouped[condition] + ", " + str(name) +\
                    ".sam"
        else:
            grouped[condition] = str(name) + ".sam"
    return grouped

def call_cuffdiff(groups, transcripts, out):
    """ Calls CuffDiff to determine differential expression.

        groups -- {condition(string):filename, filename(string)} the
                   dictionary of filenames grouped by condition.
        transcripts -- (Path) a transcripts.gtf file produced by
                       StringTie -merge.
        out -- (Path) the output folder.
    """
    samples = "" # String of samples <rep1-1.sam, rep1-2.sam> <rep2.sam>
    for condition in groups:
        samples = samples + "<" + groups[condition] + ">"
    cmd = "cuffdiff -o " + str(out) + " " + str(transcripts) + " " + samples
    subprocess.check_call(cmd, shell=True)

def call_prepde(out):
    """ Calls PrepDE.py to prepare the data from CuffDiff for DESEQ.

    Keyword Arguments:
        out -- (Path) The output folder for the previous steps.
    """
    cmd = "cd " + str(out) + " | python2.7 prepDE.py -i " + str(out) + "/lst.txt"
    subprocess.check_call(cmd, shell=True)

def call_samtoolsindex(bam_files):
    """ Calls samtools, sorts HISAT2 output, out is *.bam.bai.

    Keyword Arguments:
    bam_files -- [strings] a list of the .sam files (without extension)
                 generated by HISAT2.
    out -- (Path) the folder that has the .sam files, .bam.bai files
           will output here.
    """
    for bam in bam_files:
        if not bam.with_suffix(".bam.bai").exists():
            cmd = ("samtools index " + str(bam) + ".bam " +\
                     str(bam) + ".bam.bai")
            subprocess.check_call(cmd, shell=True)

def main():
    """ Runs when the script is called.
    """
    # Define Variables
    args = argparser()
    wk_dir = Path(args.wrk_dir)
    rawfolder = Path(args.raw_folder)
    samplenameconverter = Path(args.sample_name_converter)
    raw_files = Path(args.raw_file_names)
    out_folder = Path(args.out_folder)
    sample_names = {} # A dictionary of filenames:samples
    gtf_guide = Path(args.annotation_file)
    merged_gtf = Path(args.merged_file)
    cond_dict = {} # A ditcionary of condition:(sample1, sample2, ..., sampleN)
    files = ()

    # Generates a list of raw reads in a file
    gen_file_list(rawfolder, raw_files)

    # Fills the dictionary with run_names, sample_names
    with open(samplenameconverter) as converter:
        for line in converter:
            reg = re.search("(.*)\\t(.*)\\n", line)
            key = reg.group(1)
            value = reg.group(2)
            sample_names[key] = value

    # Call HISAT2, files is list of .sam files without extension
    with open(raw_files) as raw_path:
        files = call_hisat2(raw_path, sample_names, out_folder, rawfolder)

    # Call samtools to sort, convert, and index HISAT2 output to .bam
    # and .bam.bai files that StringTie can use
    call_samtoolssort(files, out_folder)
    call_samtoolsindex(files)

    # Call StringTie on the samtools output
    call_stringtie(files, wk_dir)

    # Creates a dictionary, putting the sample conditions together
    cond_dict = group_conditions(files)

    # Calls StringTie --merge to merge the .gtf files for CuffDiff
    separated_files = separate_tuple(files, out_folder)
    stringtie_merge(gtf_guide, merged_gtf, separated_files)

    # Calls CuffDiff to determine differential expression
    call_cuffdiff(cond_dict, merged_gtf, out_folder)

    # Calls PrepDE.py to convert CuffDiff output to what DESEQ can use
    call_prepde(out_folder)

if __name__ == "__main__":
    main()
