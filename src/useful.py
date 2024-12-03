#Function that takes in a fasta and outputs a pandas dataframe containing the protein name and the corresponding AA sequence
import pandas as pd

def fasta2df(filepath):    
    protein_names = []
    sequence = []
    with open(filepath, "r") as a_file:
        for line in a_file:
            if line.startswith('>'):
                #Get rid of the > and \n
                split_line = line.split('>')
                split_line = split_line[1].split('\n')
                prot_name = split_line[0]
                protein_names.append(prot_name)
                seq = next(a_file)
                split_seq = seq.split('\n')
                sequence.append(split_seq[0])
    #Once the list of protein names and sequences are complete, create a dataframe
    df = pd.DataFrame(list(zip(protein_names, sequence)),
               columns =['ProteinID', 'aa_sequence'])
    return df

def fasta2dict(filepath):
    protein_dict = {}
    with open(filepath, "r") as a_file:
        for line in a_file:
            if line.startswith('>'):
                #Get rid of the > and \n
                split_line = line.split('>')
                split_line = split_line[1].split('\n')
                prot_name = split_line[0]
                seq = next(a_file)
                split_seq = seq.split('\n')
                protein_dict[prot_name] = split_seq[0]
    return protein_dict


#Write the unique phyla to a fasta file
def write_fasta(filepath, df, id_col, seq_col):
    with open(filepath, 'w') as f:
        for seq in df[[id_col, seq_col]].values:
            print(seq[0])
            f.write('>' + str(seq[0]) + '\n')
            f.write(str(seq[1]) + '\n')
    f.close()
    return filepath

#Function that takes in a 3line file and outputs a pandas dataframe containing the protein names, the corresponding AA sequence, and prediction from TMBed
def read3line_tmbed(filepath):    
    protein_names = []
    sequence = []
    prediction = []
    with open(filepath, "r") as a_file:
        for line in a_file:
            if line.startswith('>'):
                #split_line = line.strip('\n')
                split_line = line.split('>')
                split_line = split_line[1].split('\n')
                prot_name = split_line[0]
                protein_names.append(prot_name)
                for a in range(2):
                    if a == 0:
                        seq = next(a_file)
                        split = seq.split('\n')
                        sequence.append(split[0])
                    else: 
                        pred = next(a_file)
                        pred_split = pred.split('\n')
                        prediction.append(pred_split[0])
    df = pd.DataFrame(list(zip(protein_names, sequence, prediction)),
               columns =['ProteinID', 'aa_sequence', 'TMBed_prediction'])
    return df

def readhmmer(filepath):
    columns = ['target_name', 'target_accession', 'query_name', 'query_accession', 'full_sequence_evalue', 'full_sequence_score', 'full_sequence_bias', 'domain_e_evalue', 'domain_score', 'domain_bias', 'domain#est_exp', 'domain#est_reg', 'domain#est_clu', 'domain#est_ov', 'domain#est_env', 'domain#est_dom', 'domain#est_rep', 'domain#est_inc', 'description_of_target']
    sep='\s+'
    df = pd.read_csv(filepath, sep=sep, skiprows=3, skipfooter=10, header=None, names=columns, index_col=False, engine='python')
    return df


def map_dive_temp(df):
    dive_temps = {'D4993_C5_LG':'35' , 'D4993_C5_H2': '72', 'D4993_C5_H1':'29', 'D4993_C5_H3':'33',
       'D4994_C39_H1':'89', 'D4994_C39_H2':'99', 'D4998_C1112_H1':'80', 'D4998_C2223_H2':'15',
       'D4998_C1112_H2':'100', 'D4991_C11_H1':'10', 'D4993_C5_H4': '39', 'D4998_C2223_H1':'4',
       'D4998_C1112_H3':'115'}
    df['temperature'] = df['dive'].map(dive_temps)
    #Convert temperature column to numeric
    df['temperature'] = pd.to_numeric(df['temperature'])
    return df

#Write out fasta files for each protein
#nylon_fp = '../../data/results/NylonHits_Guaymas2020_ALLBINS.fasta'
#pet_fp = '../../data/results/PETHits_Guaymas2020_ALLBINS.fasta'
#with open(pet_fp, 'w') as f:
#    for seq in pet_new[['target_name', 'aa_sequence']].values:
#        f.write('>' + seq[0] + '\n')
#        f.write(seq[1] + '\n')
#    f.close()