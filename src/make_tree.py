import subprocess


def make_tree(input_fasta_fp, output_directory):
    #First need to do alignment
    mafft_output = f'{output_directory}/{input_fasta_fp.split("/")[-1].split(".")[0]}.clust'
    print(f'Alignment filename will be: {mafft_output}')
    subprocess.run(['mafft', '--thread', '30', input_fasta_fp , '>' , mafft_output], shell=True)
    print('Alignment complete!')
    #Then need to do trimming to take off the biggest gaps, so that the tree is not biased towards the most gappy sequences. This is done with clipkit
    clipkit_output = mafft_output + ".clipkit"
    print(f'Trimming filename will be: {clipkit_output}')
    subprocess.run(['clipkit', mafft_output, '-m', 'gappy', '-g', '0.9']) #This will output a file with the same name as the input file, but with .clipkit appended (mafft_output.clipkit)
    print('Trimming complete! If 90%/ of the sequences are gappy, they were removed.')
    #Now run iqtree to build the phylogenetic tree
    print('Building tree...')
    subprocess.run(['iqtree', '-s', clipkit_output, '-nt', '20', '-bb', '1000', '-bnni', '-o', 'outgroup'])
    print('Tree built!')
    return 