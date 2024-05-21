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
## Complete Usage Instructions
```
equinox [-k minSeedLen] [-w bandWidth] [-c maxOcc] [-A matchScore] [-B mmPenalty]\
        [-O gapOpenPen] [-E gapExtPen] [-R RGline] [-T cutOutput] [-C commentFAST]\
        [-v verboseLevel] db.prefix reads.fq [mates.fq]
```

## Credits
