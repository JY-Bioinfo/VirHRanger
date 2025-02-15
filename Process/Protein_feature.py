import os
import numpy as np
import pandas as pd
import itertools
from Bio import SeqIO
import argparse
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn import preprocessing
import multiprocessing
from itertools import repeat

min_max_scaler = preprocessing.MinMaxScaler()
# define features names globally
AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
      "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
SC = ["1", "2", "3", "4", "5", "6", "7"]

di_pep = [''.join(i) for i in itertools.product(AA, repeat=2)]  # +400
tri_pep = [''.join(i) for i in itertools.product(AA, repeat=3)]  # +8000
di_sc = [''.join(i) for i in itertools.product(SC, repeat=2)]  # +49
tri_sc = [''.join(i) for i in itertools.product(SC, repeat=3)]  # +343
tetra_sc = [''.join(i) for i in itertools.product(SC, repeat=4)]  # +2401
myseq = "AILMVNQSTGPCHKRDEFWY"
trantab2 = myseq.maketrans("AILMVNQSTGPCHKRDEFWY", "11111222233455566777")


def parse_args():
    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--filepath",
        type=str,
        help="load in genome file, fasta or fna",
    )
    parser.add_argument(
        "--assembly_path",
        type=str,
        help="load in meta information",
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        required=True,
        help="output folder",
    )
    parser.add_argument(
        "--process",
        type=int,
        help="number of processors used",
    )

    args = parser.parse_args()

    return args


def produce_protein_feature(idx, rec_list, seq_dict):
    rec = rec_list[idx]
    ll = len(seq_dict[rec])
    # preprocessing
    seqq = seq_dict[rec].upper()
    seqqq = seqq.replace('X', 'A').replace('J', 'L').replace(
        '*', 'A').replace('Z', 'E').replace('B', 'D').replace('U', 'C').replace('O', 'C')
    X = ProteinAnalysis(seqqq)

    # +8 dim. 11201 feature dims in total
    tt = [X.isoelectric_point(), X.instability_index(), ll, X.aromaticity(),
          X.molar_extinction_coefficient(
    )[0], X.molar_extinction_coefficient()[1],
        X.gravy(), X.molecular_weight()]
    tt_n = np.asarray(tt, dtype=float)
    myseq = seqq.translate(trantab2)

    if ll > 2:
        # +400
        di_pep_count = [seqq.count(i) / (ll - 1) for i in di_pep]
        di_pep_count_n = np.asarray(di_pep_count, dtype=float)

        # +49
        di_sc_count = [myseq.count(i) / (ll - 1) for i in di_sc]
        di_sc_count_n = np.asarray(di_sc_count, dtype=float)
    else:
        di_pep_count_n = np.zeros(400)
        di_sc_count_n = np.zeros(49)

    if ll > 3:
        # +8000
        tri_pep_count = [seqq.count(i) / (ll - 2) for i in tri_pep]
        tri_pep_count_n = np.asarray(tri_pep_count, dtype=float)

        # +343
        tri_sc_count = [myseq.count(i) / (ll - 2) for i in tri_sc]
        tri_sc_count_n = np.asarray(tri_sc_count, dtype=float)
    else:
        tri_pep_count_n = np.zeros(8000)
        tri_sc_count_n = np.zeros(343)

    if ll > 4:
        # +2401
        tetra_sc_count = [myseq.count(i) / (ll - 3) for i in tetra_sc]
        tetra_sc_count_n = np.asarray(tetra_sc_count, dtype=float)
    else:
        tetra_sc_count_n = np.zeros(2401)

    # total 11201 feature dims
    cat_n = np.concatenate(
        (tt_n, di_pep_count_n, tri_pep_count_n, di_sc_count_n, tri_sc_count_n, tetra_sc_count_n))
    cat_n = cat_n.reshape((1, cat_n.shape[0]))

    return idx, cat_n


def main():
    args = parse_args()

    # load sequence into dict
    seq_dict = {}
    rec_list = []
    for file in os.listdir(args.filepath):
        if file[-3:] != 'faa' and file[-5:] != 'fasta':
            continue

        rec = file.split('.fa')[0]
        rec_list.append(rec)
        f = open(args.filepath + file, "r")
        seq_dict[rec] = f.read().replace('\n', '')

    len_Rec = len(rec_list)
    arr = np.empty((len_Rec, 11201), dtype=float)
    idx_range = range(len_Rec)

    # processing
    pool_obj = multiprocessing.Pool(args.process)  # need parallel
    idx_tuple, cat_n = zip(*pool_obj.starmap(produce_protein_feature,
                           zip(idx_range, repeat(rec_list), repeat(seq_dict))))

    for idx in idx_range:
        arr[idx_tuple[idx], :] = cat_n[idx]

    arr_min_max = min_max_scaler.fit_transform(arr)

    # feature name columns
    col_name = []
    col_name.extend(['isoelectric_point', 'instability_index', 'length', 'aromaticity',
                     'molar_extinction_coefficient_0', 'molar_extinction_coefficient_1', 
                     'gravy', 'molecular_weight'])
    # +400
    col_name.extend(di_pep)
    # +8000
    col_name.extend(tri_pep)
    # +49
    col_name.extend(['SC_'+x for x in di_sc])
    # +343
    col_name.extend(['SC_'+x for x in tri_sc])
    # +2401
    col_name.extend(['SC_'+x for x in tetra_sc])
    arr_min_max = pd.DataFrame(arr_min_max, columns=col_name)
    arr_min_max.insert(0, 'accession', rec_list)

    arr_min_max.to_hdf(args.out_folder +
                       'protein_features.h5', key='stage', mode='w')
    # arr_min_max.to_csv(args.out_folder + 'protein_features.csv', index=False)

    # metadata and dataframe
    header_agg, assembly_agg, ori_agg = [], [], []
    for filename in os.listdir(args.assembly_path):
        header_list = []
        if filename.split('.')[-1] not in ['faa', 'fasta', 'fna']:
            continue

        vid = '.'.join(filename.split('.')[:-1])
        with open(args.assembly_path + filename) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if "_prot_" in record.id:
                    header_list.append(record.id.split(
                        "_prot_", 1)[1].split('.')[0])
                elif "_protein" in record.id:
                    # New header format: >AAD12264.2_protein.
                    header_list.append(record.id.split("_protein", 1)[0])
                else:
                    # Fallback: use the entire record.id or handle as needed
                    header_list.append(record.id)

        vid_list = [vid] * len(header_list)
        assembly_agg.extend(vid_list)
        header_agg.extend(header_list)

    inputs_df = pd.DataFrame({'accession': header_agg,
                              'vid': assembly_agg,
                              'ori_seq': assembly_agg,
                              })
    list_vid_uniq = list(dict.fromkeys(assembly_agg))
    dict_temp = dict(zip(list_vid_uniq, range(len(list_vid_uniq))))
    inputs_df['ori_seq'] = inputs_df['ori_seq'].apply(lambda x: dict_temp[x])
    inputs_df.to_csv(args.out_folder + 'assembly_df.csv', index=False)

    # generate data df
    extract_df = pd.DataFrame()
    extract_df['accession'] = arr_min_max['accession']
    extract_df['proteomic_features'] = arr_min_max.iloc[:, 1:].apply(
        lambda x: '_'.join(round(x, 6).astype(str)), axis=1)
    inputs_df = inputs_df.merge(extract_df, on='accession')
    inputs_df.to_csv(args.out_folder + 'processed_input.csv', index=False)


if __name__ == "__main__":
    main()
