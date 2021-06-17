## Fusion sequence

fuseq is a command line tool that provides breakpoint information in fusion gene sequences analyzed by [Genomon](https://genomon-project.github.io/GenomonPagesR/).
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
The second argument specifies a fuseq output root directory (fuseq_root_dir) primarily for storing a breakpoint information file.
For more information about the fuseq output structure, see the output section.

```bash
$ fuseq <genomon_root_dir> <fuseq_root_dir>
```

An input data can be modified by an option.
The default input file is fusion gene detection file with filtering (merge_fusionfusion_filt.txt).
It is in a sample directory in \<genomon_root_dir/post_analysis\>.
If multiple sample directories exist, all fusion files in the samples are treated as input data.
However, the next option changes the input to a specific file instead of all files.

```bash
# The target is a user-specified file
# The first line of the file must be a header line
$ fuseq <genomon_root_dir> <fuseq_root_dir> --fusion-file </your/path/to/fusion/gene/detection/file>
```

By default, all read lines in a fusion file are given as input.
Only certain reads can be inputted by using the following option.

```bash
# The target is user-specified lines
# In the following example, the target lines are 2, 4, 5, 6, and 8
# Line numbers are separated by commas
# A Hyphen is used for sequence numbers: '4-6' is the same as '4,5,6'
# Line numbers start with 1, and the first line does not need to be specified
# because it represents a header line
$ fuseq <genomon_root_dir> <fuseq_root_dir> --lines '2,4-6,8'
```

Only data matching a specified readname or sequence can be extracted.

```bash
# Filter with readname
$ fuseq <genomon_root_dir> <fuseq_root_dir> --readname 'A_READNAME'

# Filter with sequence
$ fuseq <genomon_root_dir> <fuseq_root_dir> --sequence 'A_SEQUENCE'
```

The default input file is a filtered file (merge_fusionfusion_filt.txt).
The following option can change to an unfiltered file (merge_fusionfusion.txt).

```bash
# The target is an unfiltered file
# As with a filtered file, an unfiltered file is also in a sample directory
fuseq <genomon_root_dir> <fuseq_root_dir> --no-use-filt
```

fuseq handles three steps: Data collection for blat input, blat execution, and extracting data matching Genomon results from the blat output.
Since step1 is time consuming, there is an option to compute step1 in parallel using multiprocessing.

```bash
# Parallel processing of step1
# The default number for multiprocessing is 4
$ fuseq <genomon_root_dir> <fuseq_root_dir> --collection-processes <Number of multiprocessing>
```

The final result depends on the blat command options in step2 and the filtering settings for breakpoints in step3.
These settings can be adjusted using the following options.

```bash
# Change the options passed to blat command
# The leading hyphen required for blat options can be omitted
# The first letter of an argument value must not start with a hyphen
# This option affects step2
$ fuseq <genomon_root_dir> <fuseq_root_dir> --blat-options 'minScore=20'            # OK
$ fuseq <genomon_root_dir> <fuseq_root_dir> --blat-options '-minScore=20'           # NG
$ fuseq <genomon_root_dir> <fuseq_root_dir> --blat-options 'minScore=20 -maxGap=3'  # OK
$ fuseq <genomon_root_dir> <fuseq_root_dir> --blat-options 'minScore=20 maxGap=3'   # OK

# Expand the filtering range for breakpoints
# This option affects step3
# The following setting extends the filtering range from [tStart:tEnd] to [tStart-5:tEnd+10]
# The tStart and tEnd are alignment start and end positions obtained at blat computation, respectively
# Only breakpoints within that range are extracted
$ fuseq <genomon_root_dir> <fuseq_root_dir> --start 5 --end 10
```

fuseq has a restart option that starts from each step.
The restart is used when recomputing step2 and step3 with different balt option or filtering settings.

```bash
# Restart from blat computation (step2)
$ fuseq <genomon_root_dir> <fuseq_root_dir> --restart-blat

# Restart from filtering (step3)
$ fuseq <genomon_root_dir> <fuseq_root_dir> --restart-filter
```

However, the restart requires a working directory that fuseq outputs at runtime.
By default, the working directory will be deleted last.
Adding the following option will leave the working directory.

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

On Shirokane, collection (step1) and blat (step2) can be run as array job.
The maximum number of tasks for an array job can be set with the following options.
Note that --shirokane option is rquired when computing on Shirokane.

```bash
# Compute on Shirokane
# Available for Shirokane
$ fuseq <genomon_root_dir> <fuseq_root_dir> --shirokane

# Change the maximum number of array job tasks for collection (step1)
# Available for Shirokane
$ fuseq <genomon_root_dir> <fuseq_root_dir> --shirokane --collection-tasks 100

# Change the maximum number of array job tasks for blat (step2)
# Available for Shirokane
$ fuseq <genomon_root_dir> <fuseq_root_dir> --shirokane --blat-tasks 100
```

See help for short names of options and more options.

```bash
$ fuseq --help
```

Note that fuseq is developed based on Genomon 2.6.3, so fuseq requires the output direcotry structure and file format of that version.

### Configuration file

A configuration file ([fuseq.cfg](https://github.com/kinkalow/fuseq/blob/main/fuseq.cfg)) is used to change the default options of fuseq command.
It can be downloaded from the root of this repository.
The file is loaded by placing it in one of the following locations.
  1. The current directory where fuseq is executed
  1. The path specified by the FUSEW_CFG environment variable
  1. ~/.config/fuseq.cfg

fuseq checks for configuration files in order from the top and reads only the first one found.

The contents of the configuration file are described below.
The lines starting with # represent comments.
The key names must be the same as the option names.
However, the names omit the leading double hyphen and replace the rest of the hyphen with an underbar.
The values are the default value of the options.

### Output

The output root directory for fuseq has a single or multiple sample directories.
The directory names are the same as Genomon output samples (\<genomon_root_dir/post_analysis/samples\>).
Each sample directory consists of a breakpoint information file (fusion_sequence.txt), input files (fusion.txt and star), and a working directory (work_restart).
The input files are symbolic links to a fusion file and a star directory in a Genomon output directory unless --lines option is used.
The working directory includes the intermediate files required for restart.

Here is a description of the breakpoint information file.
Its output consists of 6 lines for each readname and has the following format:

```
[readname] [fusionLineNr=<line number of fusion.txt>]                                          |
[chromesomeA] [breakpointA] [strandA] [gene overlappingA] [exon-intron junction overlappingA]  |  information in Genomon output
[chromesomeB] [breakpointB] [strandB] [gene overlappingB] [exon-intron junction overlappingB]  |
[fusion sequence]                                                                              |
[sequence for chromesomeA and breakpointA] [sequence range] [breakpoint range] [chromesomeA] [strandA2]  |  information from blat computation
[sequence for chromesomeB and breakpointB] [sequence range] [breakpoint range] [chromesomeB] [strandB2]  |

next readname
...
```

Note that strandA2 and strandB2 are not necessarily equal to strandA and strandB, respectively.
The following is a concrete example of showing breakpoint information about the fusion of a chr3 sequence and another chr3 one.

```
SomeReadName fusionLineNr=2
3 197602647 + LRCH3 LRCH3.start
3 197592293 - ENST00000536618.1;ENST00000414675.2;ENST00000438796.2;ENST00000441090.2;ENST00000425562.2 ENST00000536618.1.end;ENST00000441090.2.end;ENST00000414675.2.end;ENST00000425562.2.end;ENST00000438796.2.end
AGCTGTACCTAAATTAACAATGGCGAAATGCAGGCGAAATGTGGAAAATTTCCTAGAAGCTTGCAGAAAAATTGGTGTACCTCAGAGTGATGACAGACCTAATGCTCTATTAAGTTCACCTGCAAC
       CCTAAATTAACAATGGCGAAATGCAGGCGAAATGTGGAAAATTTCCTAGAAGCTTGCAGAAAAATTGGTGTACCTCAG                                          [8,85]   [197602569,197602646] chr3 +
                                                                                  CAGAGTGATGACAGACCTAATGCTCTATTAAGTTCACCTGCAAC [83,126] [197592291,197592334] chr3 +
       <-- seq for chrA(=3) and bpA(=197602647) ------------------------------------>
                                                                                  <-- seq for chrB(=3) and bpB(=197592293) -->
```

The first three lines show the data from Genomon result, and the fourth and fifth lines show the Blat computation result.
The sequence on line 5 corresponds to chromesomeA and breakpointA, and the sequence on line 6 corresponds to chromesomeB and breakpointB.
The lines 5 and 6 print only if the chromesomes obtained by Genomon and Blat are consistent and the breakpoint obtained by Genomon is close to the breakpoint range obtained by Blat.
Therefore, for a fusion sequence where blat and Genomon results do not match, the output consists of three lines.
These three lines are output at the end of the file (fusion_sequence.txt).
