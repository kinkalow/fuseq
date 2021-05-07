## Fusion sequence

fuseq is a command line tool that provides breakpoint information for fusion gene sequences analyzed by [Genomon](https://genomon-project.github.io/GenomonPagesR/).
That breakpoints can be obtained by [Blat](https://genome.ucsc.edu/cgi-bin/hgBlat).
So, fuseq first collects the readnames and their sequences from the Genomon analysis results and passes them as input data for Blat.
Blat outputs multiple candidates, but fuseq only extracts candidates that match Genomon results and finally writes breakpoint information to a file.

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

The fuseq command specifies the output root directory of Genomon.
By default, the target is all data written in the fusion gene detection file (merge_fusionfusion_filt.txt) in that directory.
If you want to target specific data, write the target data in a file and optionally specifying the file.
Note that fuseq is developed based on Genomon 2.6.3, so fuseq uses the file format of that version.
Besides, by default, the assembly database required at Blat runtime uses GRCh37 installed on [Shirokane](https://gc.hgc.jp/en/).
The database can optionally be changed if necessary.

```bash
# The first argument specifies the output root directory of Genomon.
# The second argument specifies a root directory to store the result files of fusion sequences.
$ fuseq <genomon_out> <fuseq_out>

# Use --reference option to change the reference genome file.
$ fuseq <genomon_out> <fuseq_out> --reference </your/path/to/reference/genome>

# Use --fusion-file option to use a different file instead of a fusion gene detection file in Genomon output directory.
$ fuseq <genomon_out> <fuseq_out> --fusion-file </your/path/to/fusion/gene/detection/file>

# Use --blat-opt option to change Blat options.
$ fuseq <genomon_out> <fuseq_out> --blat-opt '-minScore=20 -minMatch=1'

# Use --restart-blat option to restart from Blat computation.
# On restart, the working directory created by fuseq must exist in <fuseq_out>.
# By default, the working directory will be deleted, but you can leave it by using the --no-delete-work option.
$ fuseq <genomon_out> <fuseq_out> --restart-blat --no-delete-work

# See help for more options.
$ fuseq --help
```

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
