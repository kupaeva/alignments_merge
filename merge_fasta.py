import argparse
import pandas as pd
import sys


# this function parse arguments from command line
def parser():
    parser = argparse.ArgumentParser(description=
                                     'It is a script for phylogenetic alignments'
                                     ' merging. Here it reads fasta files from IN'
                                     ' paths, merge it by name, and write to the'
                                     ' OUT file. Additionally you can create nexus'
                                     ' file with the length parameters of your'
                                     ' sequences.\n\n',
                                     epilog='Example of usage: '
                                     '<python merge_fasta.py -in ./test/input_1.fasta'
                                     ' ./test/input_2.fasta -out ./test/output.fasta'
                                     ' -nexus ./test/output.nexus>')
    parser.add_argument('-in', nargs='*',
                        help='paths to all alignment files in fasta format, which you want merge')
    parser.add_argument('-out', help='path to output fasta file with merged sequences')
    parser.add_argument('-nexus', help='path to output nexus file with length parameters of your sequences')
    args = parser.parse_args()
    return args


# this function read and processed fasta-files from IN
def file_read(path):
    dna_list = {}
    with open(path, 'r') as file:
        for line in file:
            line = line.rstrip('\n')
            if '>' in line:
                if line not in dna_list:
                    dna_list[line] = str()
                    temp_taxon = line
                else:
                    print(path, 'contain duplicated name', line)
                    raise NameError()
            else:
                dna_list[temp_taxon] = dna_list[temp_taxon] + line
    dna_table = pd.DataFrame.from_dict(dna_list, orient='index')
    dna_table.reset_index(inplace=True)
    dna_table.columns = ['taxon', path.split('\\')[-1]]
    return dna_table


args = parser()
path_input = vars(args)['in']

# read files and merge all alignments to one dataframe
for path in path_input:
    one_gene_table = file_read(path)
    if 'all_genes_table' not in globals():
        all_genes_table = one_gene_table
    else:
        all_genes_table = all_genes_table.merge(one_gene_table, left_on='taxon', right_on='taxon', how='outer')

# compute length of each gene and fill NA by '-'
length = all_genes_table.iloc[1].str.len()[1:].to_dict()
all_genes_table['all_genes'] = str()
for gene in length.keys():
    all_genes_table[gene].fillna(length[gene] * '-', inplace=True)
    all_genes_table['all_genes'] = all_genes_table['all_genes'] + all_genes_table[gene]

# write OUT and NEXUS files
final_table = all_genes_table[['taxon', 'all_genes']]
if vars(args)['out']:
    final_table.to_csv(vars(args)['out'], sep='\n', index=False, header=False)
else:
    print(final_table)

if vars(args)['nexus']:
    start = 0
    original_stdout = sys.stdout
    with open(vars(args)['nexus'], 'w') as f:
        sys.stdout = f
        print('#nexus\nbegin sets;')
        for gene in length.keys():
            print('\tCHARSET ', gene, ' = ', start + 1, '-', length[gene] + start, ';', sep='')
            start += length[gene]
        print('end;\n')
        sys.stdout = original_stdout
