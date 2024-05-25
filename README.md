# Equinox-Alignment
Equinox Alignment is an alignment tool that matches sequence reads with reference genome data, and functions similarly to [BWA MEM](https://github.com/lh3/bwa).

## Installation
Install in terminal:
```
git clone https://github.com/ianfmcn/Equinox-Alignment
cd equinox; make
```
Install using PIP:
```
pip install Equinox-Alignment
```

## Getting Started
Equinox is made to have similar usage to [BWA MEM](https://bio-bwa.sourceforge.net/bwa.shtml), which takes as inputs a reference genome file (.fa) and fastq read files (.fastq), and outputs a sam file (.sam).
```
equinox reference.fa read1.fq read2.fq > output.sam
```
Example fastq files can be found here: [Fastq Files](https://drive.google.com/drive/folders/1PVqUAGe60cw056kn5xN-kyurCVh9kETV?usp=sharing)

Basic test files can be found in the example_files folder.

## Complete Usage Instructions
```
equinox [-m matchScore] [-s mmPenalty] [-d indelPenalty] [-b bandWidth]\
        [-g gapOpenPen] [-e gapExtPen] [-k minSeedLen] [-c maxOcc] [-R RGline]\
        [-t cutOutput] [-C commentFAST] [-v verboseLevel] db.prefix reads.fq [mates.fq]
```

## Usage Instructions for Progress Report
```
Make parsing.py from reference_code executable via chmod +x parsing.py
./parsing.py can be run with the arguments [-m, -s, -d, -g, -e, -b]

Note: file parsing to obtain reads from .fa/.fq files has not been fully implemented yet; running ./parsing.py will return alignment of a placeholder pair of sequences.

The alignment code for BandedAlignment can be tested for single end alignment by running that py script directly using the example files. (It currently contains a placeholder main function to parse arguments.)
python3 ./reference_code/BandedAlignment.py ./example_files/test_reference.fa ./example_files/test_sequence.fq 5 -m 1 -s -1 -d -1 > test_banded.txt
This test_banded.txt is also in the folder example_files.
```

## Credits
