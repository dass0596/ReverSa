from argparse import ArgumentParser
import ConfigParser

def arguments():
    parser = ArgumentParser()
    config = ConfigParser.ConfigParser()
    
    #Input files
    parser.add_argument('--fasta_file')
    parser.add_argument('--config_file')
    args = parser.parse_args()
    #Read config_file
    config.read(args.config_file)
    #Create a dictionary to store values of the config_file
    config_d = {}
    config_d['win_length'] = int(config.get('Preprocessing','window_length'))
    config_d['win_step'] = int(config.get('Preprocessing','window_step'))
    config_d['score_adj'] = float(config.get('Similarity','score_adjacency'))
    config_d['kmer_len'] = int(config.get('Composition','kmer_length'))
    config_d['pca_comp'] = float(config.get('PCA','variance_components'))
    return args,config_d
    