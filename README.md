## Fusion sequence

fuseq is a command line tool that provides a visual representation of where the breakpoints are located in the fusion gene sequences analyzed by [Genomon](https://genomon-project.github.io/GenomonPagesR/).
That breakpoints can be obtained by [Blat](https://genome.ucsc.edu/cgi-bin/hgBlat).
So, fuseq first collects the readnames and sequences from a Genomon result and passes them as input data for Blat.
Blat outputs multiple candidates, but fuseq only extracts candidates that match the Genomon result and finally writes the breakpoint information to a file.

### Requires

- Python >= 3.6
- Blat

### Installation

```bash
$ git clone https://github.com/kinkalow/fuseq.git
$ cd fuseq
$ python3 setup.py install --user
$ fuseq --version  # check if the installation was successful
```

### Basic Usage

fuseq command requires two arguments.
The first argument specifies the output root directory of Genomon (genomon_root_dir).
The second argument specifies a fuseq output root directory for storing breakpoint information files and fuseq working directories.

```bash
$ fuseq <genomon_root_dir> <fuseq_root_dir>
```

fuseq can optionally change input data.
The default input files are fusion gene detection files with filtering (merge_fusionfusion_filt.txt).
The files are in genomon_root_dir/post_analysis/\<single or multiple samples\>.
If multiple sample directories exist, all files in them are used as input data.
If you want to input only certain file, the following option allows you to change the target.

```bash
# The target is a user-specified file
# The first line of the file must be a header line
$ fuseq <genomon_root_dir> <fuseq_root_dir> --fusion-file </your/path/to/fusion/gene/detection/file>
```

By default, all lines in the files are given as input.
If you want to focus only specific lines, use the following option.

```bash
# The target is user-specified lines
# In the following example, the target lines are 2, 4, 5, 6, and 8
# Line numbers are separated by commas
# A Hyphen is used for sequence numbers: '4-6' is the same as '4,5,6'
# Line numbers start with 1, and the first line does not need to be specified
# because it represents the header line
$ fuseq <genomon_root_dir> <fuseq_root_dir> --lines '2,4-6,8'
```

The default is to input filtered files.
If you want to change unfiltered files, use the following option.

```bash
# The targets are unfiltered files named merge_fusionfusion.txt
# They are in the same directory as filted files named merge_fusionfusion_filt.txt
fuseq <genomon_root_dir> <fuseq_root_dir> --no-use-filt
```

fuseq handles three steps: Data collection for blat input, blat execution, and extracting data matching Genomon results from the blat output.
The final output depends on blat options or filtering settings.

```bash
# Change the options passed to blat command
# This option affects step2 results
$ fuseq <genomon_root_dir> <fuseq_root_dir> --blat-opt '-minScore=20 -minMatch=1'

# Expand the search range when filtering on breakpoint
# The following setting extends the search from [tStart:tEnd] to [tStart-5:tEnd+10]
# where tStart and tEnd are alignment start and end positions obtained at blat computation, respectively
# This option affects step3 results
$ fuseq <genomon_root_dir> <fuseq_root_dir> --start 5 --end 10
```

To allow recomputation with different options, fuseq has a restart function which can be started from each step.

```bash
# Restart from blat computation (step2)
$ fuseq <genomon_root_dir> <fuseq_root_dir> --restart-blat

# Restart from filtering blat results (step3)
$ fuseq <genomon_root_dir> <fuseq_root_dir> --restart-filter
```

However, the restart requires a working directory that fuseq outputs at runtime.
By default, the working directory will be deleted last.
Adding the following option will leave the working directory in \<fuseq_root_dir/single or multiple samples\>.

```bash
# Don't delete working directory
$ fuseq <genomon_root_dir> <fuseq_root_dir> --no-delete-work
```

Blat computation requires assembly database.
By default, fuseq uses GRCh37 installed on [Shirokane](https://gc.hgc.jp/en/) supercomputer.
The database can optionally be changed.

```bash
# Change reference data
$ fuseq <genomon_root_dir> <fuseq_root_dir> --reference </your/path/to/reference/genome>
```

The process of step1 (collecting input files for blat) is time-consuming.
The following option allows you to use multiprocessing to process step1 in parallel.

```bash
# Parallel processing of step1
# The default number for multiprocessing is 4
$ fuseq <genomon_root_dir> <fuseq_root_dir> --coll-procs <Number of multiprocessing>
```

See help for short names of options and more options.

```bash
$ fuseq --help
```

Note that fuseq is developed based on Genomon 2.6.3, so fuseq requires the output direcotry structure and file format of that version.

### Output format

The fuseq output root directory (fuseq_root_dir) has sample directories whose names are the same as Genomon output samples (<genomon_root_dir/post_analysis/samples>).
Each sample directory consists of a breakpoint information file (fusion_sequence.txt), input data (fusion.txt and star), and a working directory (work).

Here is a description of the breakpoint information file.
Its output consists of 5 lines for each readname and has the following format:

```
readname fLineNr=<line number of fusion.txt>                     |
chromesomeA breakpointA strandA chromesomeB breakpointB strandB  | information in Genomon output files
original fusion sequence                                         |
[sequence for chromesomeA and breakpointA] [string number range] [breakpoint range] [chromesomeA] [strandA2]  | information from blat computation
[sequence for chromesomeB and breakpointB] [string number range] [breakpoint range] [chromesomeB] [strandB2]  | strandA2 and strandB2 are not necessarily equal to strandA and strandB, respectively

next readname
...
```

The following is a concrete example of showing breakpoint information about the fusion of a chr3 sequence and another chr3 one.

```
SomeReadName fLineNr=2
3 197602647 + 3 197592293 -
AGCTGTACCTAAATTAACAATGGCGAAATGCAGGCGAAATGTGGAAAATTTCCTAGAAGCTTGCAGAAAAATTGGTGTACCTCAGAGTGATGACAGACCTAATGCTCTATTAAGTTCACCTGCAAC
       CCTAAATTAACAATGGCGAAATGCAGGCGAAATGTGGAAAATTTCCTAGAAGCTTGCAGAAAAATTGGTGTACCTCAG                                          [8,85]   [197602569,197602646] chr3 +
                                                                                  CAGAGTGATGACAGACCTAATGCTCTATTAAGTTCACCTGCAAC [83,126] [197592291,197592334] chr3 +
       <-- seq for chrA(=3) and bpA(=197602647) ------------------------------------>
                                                                                  <-- seq for chrB(=3) and bpB(=197592293) -->
```

The first three lines show the data from Genomon result, and the fourth and fifth lines show the Blat computation result.
The sequence on line 4 corresponds to chromesomeA, breakpointA and strandA on line 2, and the sequence on line 5 corresponds to chromesomeB, breakpointB and strandB.
The lines 4 and 5 print only if the chromesomes obtained by Genomon and Blat are consistent and the breakpoint obtained by Genomon is close to the breakpoint range obtained by Blat.
For data where the results of Genomon and Blat do not match, their data are collected and output at the end of a file (fusion_sequence.txt).
