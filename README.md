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
The second argument specifies a fuseq output directory for storing files related to the breakpoint information of fusion sequences.

```bash
$ fuseq <genomon_root_dir> <fuseq_out>
```

The default input files for fuseq are fusion gene detection files with filtering (merge_fusionfusion_filt.txt).
The files are in genomon_root_dir/post_analysis/\<single or multiple samples\>.
Note that the all files in the sample directories are used as input data.
If you want to input only certain file instead of all files, the following option allows you to change the target.

```bash
# The target file will be user-specified file
$ fuseq <genomon_out> <fuseq_out> --fusion-file </your/path/to/fusion/gene/detection/file>

# The target files will be merge_fusionfusion.txt instead of merge_fusionfusion_filt.txt
# merge_fusionfusion.txt is in the same directory as merge_fusionfusion_filt.txt
fuseq <genomon_root_dir> <fuseq_out> --no-use-filt
```

fuseq handles three steps: Data collection for blat input, blat execution, and extracting data matching Genomon results from the blat output.
The final output depends on blat options or filtering settings, so fuseq has restart options that starts from each step.

```bash
# Change the options passed to blat command
$ fuseq <genomon_out> <fuseq_out> --blat-opt '-minScore=20 -minMatch=1'

# Expand the search range when filtering on breakpoint
# The following setting extends the search from [tStart:tEnd] to [tStart-5:tEnd+10]
# where tStart and tEnd are the start and end positions of alignment obtained at blat computation, respectively
$ fuseq <genomon_out> <fuseq_out> --start 5 --end 10

# Restart from blat computation (2step)
$ fuseq <genomon_out> <fuseq_out> --restart-blat

# Restart from filtering blat results (3step)
$ fuseq <genomon_out> <fuseq_out> --restart-filter
```

However, the restart requires a working directory that fuseq outputs at runtime.
By default, the working directory will be deleted last.
Adding the following option will leave the working directory in <fuseq_out>.

```bash
# Don't delete working directory
$ fuseq <genomon_out> <fuseq_out> --no-delete-work

# Start with blat computation without deleting working directory
$ fuseq <genomon_out> <fuseq_out> --no-delete-work --restart-blat
```

Blat computation requires assembly database.
By default, fuseq uses GRCh37 installed on [Shirokane](https://gc.hgc.jp/en/) supercomputer.
The database can optionally be changed.

```bash
# Change reference data
$ fuseq <genomon_out> <fuseq_out> --reference </your/path/to/reference/genome>
```

See help for short names of options and more options.

```bash
$ fuseq --help
```

Note that fuseq is developed based on Genomon 2.6.3, so fuseq requires the output direcotry structure and file format of that version.

### Output format

fuseq outputs the breakpoint information of the fusion sequences for each read name.
Its output format is at least 3 lines for each read as follows:

```
readname
chromesomeA breakpointA strandA chromesomeB breakpointB strandB
one or more lines of fusion sequences

next readname
...
```

Note that only when Genomon fusion sequence is also found in Blat, fusion sequence lines are displayed in multiple lines to show the break position.
For example, the following shows a sequence on chromosome 9 (=chromosomeA) is fused to a sequence on chromosome 2 (=chromosomeB).
The upper sequence is the information for chromosomeA and the lower sequence is the information for chromosomeB.

```
SOME_READNAME
9 131833825 + 2 216261859 -
GTTCTGCTGTGGCTCTGCCCTTTCCAGGTTGAGAGGCCTG
                                    CCTGGGGTGGTGCTCCTCTCCCAGGAGACTGTGAGCACTCCAGTGTCAGGGTTTGCCTCCAGATGCAAGTTTGTTGGTGGAGACAATGGT
<--- sequence for chromesomeA(=9) ----->
                                    <--- sequence for chromesomeB(=2) ------------------------------------------------------->
```

A single-line sequence indicates that Genomon fusion sequence is not found in Blat candidates.
The list found in Blat is displayed before the list not found.
