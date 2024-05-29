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
Equinox is made to have similar usage to [BWA MEM](https://bio-bwa.sourceforge.net/bwa.shtml), which takes as inputs a reference genome file (.fa) and fastq read files (.fq), and outputs a sam file (.sam).
```
equinox reference.fa read1.fq read2.fq -o output.sam
```
Basic test files can be found in the example_files folder.

## Complete Usage Instructions
```
equinox -m matchScore -s mmPenalty -d indelPenalty [-b bandWidth]\
         reference.fa reads.fq [mates.fq] -o output.sam
```

## Usage Instructions for Progress Report
```
Make parsing.py from reference_code executable via chmod +x parsing.py
./parsing.py can be run with the arguments [-m, -s, -d, -b]

Note: file parsing to obtain reads from .fa/.fq files has not been fully implemented yet; running ./parsing.py will return alignment of a placeholder pair of sequences.

The alignment code for BandedAlignment can be tested for single end alignment by running that py script directly using the example files. (It currently contains a placeholder main function to parse arguments.)
python3 ./reference_code/BandedAlignment.py ./example_files/test_reference.fa ./example_files/test_sequence.fq 5 -m 1 -s -1 -d -1 > test_banded.txt
This test_banded.txt is also in the folder example_files.

The alignment code for Local Alignment can be tested for single end alignment by running the py script described below. (test, test1, and test2.txt are the test files for this alignment, with the first string being the read, and the second the genome.)
py ./reference_code/local.py ./example_files/test.txt -m 1 -s -1 -d -1 -a

Note 05/26/24:
We have worked on this project further, and many of the files and test files have been altered. The code can no longer be run as described above. Please run the commands below to see our project in action.

IN PYTHON:
python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -o ./example_files/test_local.txt
python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -b 5 -o ./example_files/test_banded.txt

EXECUTABLE:
chmod +x parsing.py
./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -o ./example_files/test_local.txt
./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -b 5 -o ./example_files/test_banded.txt
```

## Credits
