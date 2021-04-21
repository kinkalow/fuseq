## Fusion sequence

fuseq is a command line tool that allows you to get breakpoint information for fusion gene sequences obtained by [Genomon](https://genomon-project.github.io/GenomonPagesR/).
That breakpoint information can be obtained by [Blat](https://genome.ucsc.edu/cgi-bin/hgBlat).
So, fuseq first examines the readnames and their sequences from the Genomon analysis results and passes them as input data for Blat.
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

Since fuseq has been developed for Genomon(v2.6.3) installed on [Shirokane](https://gc.hgc.jp/en/) supercomputer, fuseq uses the output data of that Genomon version by default.
Also, the assembly database required at Blat runtime uses GRCh37 installed on Shirokane.
However, these output data files and database can be changed by options.

```bash
# The first argument specifies the output root directory of Genomon.
# The second argument specifies a root directory to store the result files of fusion sequences.
$ fuseq <genomon_out> <fuseq_out>

# Use reference option to change the reference genome file.
$ fuseq <genomon_out> <fuseq_out> --reference </your/path/to/reference/genome>

# See help for more options.
$ fuseq --help
```

### Output format

fuseq outputs the breakpoint information of the fusion sequences for each read name.
Its output format is at least 3 lines for each read as follows:

```
readname
chromesome1 breakpoint1 strand1 chromesome2 breakpoint2 strand2
one or more lines of fusion sequences

next readname
...
```

Note that only when Genomon fusion sequence is also found in Blat, fusion sequence lines are displayed in multiple lines to show the break position.
A single-line sequence indicates that Genomon fusion sequence is not found in Blat candidates.
The list found in Blat is displayed before the list not found.
