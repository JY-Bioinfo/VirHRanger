import os
import pandas as pd
import argparse
import numpy as np
import pickle
import json


def np_relu(x):
    return (x > 0) * x


def violation_process(hierarchy_structure, probs):
    # skip if binary
    if hierarchy_structure == 'binary':
        return probs
    taxa_list = list(hierarchy_structure.keys())
    taxa_roots = range(min(taxa_list))

    output = np.zeros(probs.shape)

    output[:, taxa_roots] = probs[:, taxa_roots]

    for child in taxa_list:
        parent = hierarchy_structure[child]

        p_child = probs[:, child]
        p_parent = output[:, parent]
        output[:, child] = p_child - np_relu(p_child - p_parent)

    # sum for all the labels
    return output


def parse_args():

    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--ensemble_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--all_cates",
        type=str,
    )
    parser.add_argument(
        "--hierarchy_structure",
        type=str,
        default="binary"
    )
    parser.add_argument(
        "--cutoff",
        type=str,
        help="cutoff_path"
    )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    with open(args.all_cates, 'rb') as fp:
        all_cates = pickle.load(fp)
    with open(args.cutoff, 'rb') as fp:
        cutoffs = pickle.load(fp)

    if args.hierarchy_structure != 'binary':
        f = open(args.hierarchy_structure)
        taxa = json.load(f)
        hierarchy_structure = {int(k): int(v) for k, v in taxa.items()}
    else:
        hierarchy_structure = 'binary'
    '''check folder situations'''
    rec_exist = {}
    for mode in ['Genomic', 'Genomic_Protein', 'Genomic_Protein_PPI']:
        rec_exist[mode] = 0
        dirlen = 0
        for files in os.listdir(args.ensemble_path + f'{mode}'):
            if files.endswith('.pkl'):
                dirlen += 1

        if dirlen != 0:
            rec_exist[mode] = 1

    VIDS, thres, Probs, Preds = {}, {}, {}, {}

    for mode in ['Genomic', 'Genomic_Protein', 'Genomic_Protein_PPI']:
        if rec_exist[mode] == 1:
            '''load vid'''
            with open(args.ensemble_path + f'{mode}/vids.pkl', 'rb') as fp:
                VIDS[mode] = pickle.load(fp)

            '''load cutoff'''
            thres[mode] = cutoffs[mode]

            '''load probs'''
            Probs[mode] = pd.read_csv(
                args.ensemble_path + f'{mode}/pred_probs.csv')
    
    if rec_exist['Genomic_Protein'] == 1:
        mode = 'Genomic_Protein'
        
        '''special case for Genomic_Protein_PPI'''
        if rec_exist['Genomic_Protein_PPI'] == 1:
            '''update Genomic_Protein ensembled model with additional PPI information'''
            thres['Genomic_Protein_PPI_updated'] = thres['Genomic_Protein']
            # using stricter threshold will not violate hierarchy
            if thres['Genomic_Protein_PPI']['6_Homo sapiens'] > thres['Genomic_Protein']['6_Homo sapiens']:
                thres['Genomic_Protein']['6_Homo sapiens'] = thres['Genomic_Protein_PPI']['6_Homo sapiens']
            
            '''set index first'''
            Probs['Genomic_Protein'].insert(0, 'vid', VIDS['Genomic_Protein'])
            Probs['Genomic_Protein'] = Probs['Genomic_Protein'].set_index('vid')
            Probs['Genomic_Protein_PPI'].insert(0, 'vid', VIDS['Genomic_Protein_PPI'])
            Probs['Genomic_Protein_PPI'] = Probs['Genomic_Protein_PPI'].set_index('vid')
            temp_probs = Probs[mode]
            temp_probs.update(Probs['Genomic_Protein_PPI'])

            '''index removed'''
            temp_probs = violation_process(
                hierarchy_structure, temp_probs.values)
            Probs[mode] = pd.DataFrame(temp_probs, columns=all_cates)
            Preds[mode] = np.where(Probs[mode] > thres['Genomic_Protein_PPI_updated'], 1, 0)
            Preds[mode] = violation_process(hierarchy_structure, Preds[mode])

            Preds[mode] = pd.DataFrame(Preds[mode], columns=all_cates)
            Probs[mode].insert(0, 'vid', VIDS[mode])
            Preds[mode].insert(0, 'vid', VIDS[mode])
            Preds[mode].insert(1, 'ensemble', 'nuc_orf')
            Preds[mode].loc[Preds[mode].vid.isin(
                VIDS['Genomic_Protein_PPI']), 'ensemble'] = 'Genomic_Protein_PPI'

        else:
            Preds[mode] = np.where(Probs[mode] > thres[mode], 1, 0)
            Preds[mode] = violation_process(hierarchy_structure, Preds[mode])
            Preds[mode] = pd.DataFrame(Preds[mode], columns=all_cates)
            Probs[mode].insert(0, 'vid', VIDS[mode])
            Preds[mode].insert(0, 'vid', VIDS[mode])
            Preds[mode].insert(1, 'ensemble', 'Genomic_Protein')

    if rec_exist['Genomic'] == 1:
        mode = 'Genomic'
        Preds[mode] = np.where(Probs[mode] > thres[mode], 1, 0)
        Preds[mode] = violation_process(hierarchy_structure, Preds[mode])
        Preds[mode] = pd.DataFrame(Preds[mode], columns=all_cates)
        Probs[mode].insert(0, 'vid', VIDS[mode])
        Preds[mode].insert(0, 'vid', VIDS[mode])
        Preds[mode].insert(1, 'ensemble', 'Genomic')

    Probs_final, Preds_final = pd.DataFrame(), pd.DataFrame()
    for mode in ['Genomic', 'Genomic_Protein']:
        if rec_exist[mode] == 1:
            Probs_final = Probs_final.append(Probs[mode], ignore_index=True)
            Preds_final = Preds_final.append(Preds[mode], ignore_index=True)

    Probs_final = Probs_final.sort_values('vid')
    Preds_final = Preds_final.sort_values('vid')

    Probs_final.to_csv(args.output_path + 'whole_probs.csv',
                       index=False, sep=',')
    Preds_final.to_csv(args.output_path + 'whole_preds.csv',
                       index=False, sep=',')

    with open(args.output_path + "whole_results.txt", "w") as text_file:

        for i in range(len(Preds_final)):
            text_file.write(
                str(Preds_final.iloc[i, 0]) + ': '+str(Preds_final.iloc[i, 1]) + '\n')
            text_file.writelines(
                ', '.join(all_cates[np.where(Preds_final.iloc[i, 2:] == 1)].to_list()))
            text_file.write('\n\n')


if __name__ == "__main__":
    main()
