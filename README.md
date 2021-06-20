# merge fasta
Python script for merging few alignments in FASTA format. 

It is a script for phylogenetic alignments merging. 
Here it reads fasta filesfrom IN paths, merge it by name, and write to the OUT file. 
Additionally you can create NEXUS file with the length parameters of your sequences.

## Options:
    -h, --help         show this help message and exit
  
    -in [IN [IN ...]]  paths to all alignment files in fasta format, which you
                     want merge
                     
    -out OUT           path to output fasta file with merged sequences
  
    -nexus NEXUS       path to output nexus file with length parameters of your
                     sequences

## Example of usage: 
```bash
python merge_fasta.py -in ./test/input_1.fasta ./test/input_2.fasta -out ./test/output.fasta -nexus ./test/output.nexus
```

## Required:
* python >3.6
* pandas


## Autor:
Daria Kupaeva: d.kupaeva@gmail.com
