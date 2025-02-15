import pandas as pd
import argparse
import numpy as np
import os
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm


def produce_score(df, test_id):
    '''input'''

    support_df = pd.DataFrame(columns=df.filter(
        regex=("[0-9]_")).columns.insert(0, 'qseqid'))

    for tid in tqdm(test_id):

        temp_df = df[df['qseqid'] == tid]
        pair_score = temp_df['pident']
        pair_score = pair_score/pair_score.sum()

        # 1 keep label cols
        # 2 multiply support
        # 3 sum up
        support = temp_df.filter(
            regex=("[0-9]_")).mul(pair_score, axis=0).sum(axis=0)
        support_df = support_df.append(support, ignore_index=True)

    support_df['qseqid'] = test_id

    return support_df


def PN_saver(dictex):
    for key, val in dictex.items():
        val.to_csv("data_{}.csv".format(str(key)), index=False)

    with open("keys.txt", "w+") as f:  # saving keys to file
        f.write(str(list(dictex.keys())))


def PN_loader(PN_PATH):
    os.chdir(PN_PATH)
    """Reading data from keys"""
    with open("keys.txt", "r") as f:
        keys = eval(f.read())

    dictex = {}
    for key in keys:
        dictex[key] = pd.read_csv("data_{}.csv".format(str(key)))

    return dictex


def parse_args():
    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--filepath",
        required=True,
        type=str,
        help="load in Phylo features",
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        required=True,
        help="output folder",
    )
    parser.add_argument(
        "--reference_path",
        type=str,
        required=True,
        help="label reference",
    )
    parser.add_argument(
        "--query_path",
        type=str,
        required=True,
        help="get query id list",
    )

    parser.add_argument(
        "--support_num",
        default=5,
        type=int,
        help="number of processors used",
    )

    parser.add_argument(
        "--process",
        type=int,
        help="number of processors used",
    )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    Header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
              'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    blast_out = pd.read_csv(
        args.filepath + 'blast_out.txt', sep='\t', names=Header)
    blast_out = blast_out[blast_out.qseqid != blast_out.sseqid]
    blast_out = blast_out.dropna(axis=0)
    blast_out = blast_out.sort_values(by='bitscore').groupby(
        ["qseqid"]).tail(args.support_num)
    blast_out = blast_out.reset_index(drop=True)

    # import reference labels
    DF_labels = pd.read_csv(args.reference_path).drop(['label'], axis=1)
    DF_labels.rename({'VIV_accession': 'sseqid'}, axis='columns', inplace=True)

    # get query_ids
    query_id = []
    with open(args.query_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            query_id.append(record.id)

    DF_separate = pd.merge(blast_out, DF_labels, on='sseqid')
    support_df = produce_score(DF_separate, query_id)

    Path(args.out_folder).mkdir(parents=True, exist_ok=True)
    np.save(args.out_folder + 'pred_probs.npy', support_df.iloc[:, 1:].values)

    support_df.rename({'qseqid': 'vid'}, axis='columns', inplace=True)
    support_df.to_csv(args.out_folder + 'whole_probs.csv', index=False)


if __name__ == "__main__":
    main()
