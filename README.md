# Equinox-Alignment
Equinox Alignment is an alignment tool that matches sequence reads with reference genome data, and functions similarly to [BWA MEM](https://github.com/lh3/bwa).

## Installation
In terminal:
```
git clone https://github.com/ianfmcn/Equinox-Alignment
cd Equinox-Alignment
```
Install Required Packages:\
In terminal:
```
pip install Biopython
pip install Pandas
pip install argparse
```
## Complete Usage Instructions
```
python3 parsing.py reference.fa read.fq -m matchScore -s mmPenalty -d indelPenalty [-b bandWidth] -o output.sam
```

## Usage Examples
IN PYTHON:\
For Local Alignment:
```
python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -o ./example_files/test_local.sam
```
For Banded Alignment:
```
python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -b 3 -o ./example_files/test_banded.sam
```

EXECUTABLE:\
First, run:
```
chmod +x parsing.py
```
Then, run either of the following
```
./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -o ./example_files/test_local.sam
./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -b 5 -o ./example_files/test_banded.sam
```
