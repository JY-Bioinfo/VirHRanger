import os
from Bio import SeqIO
import pandas as pd
import argparse
SEED = 42


def parse_args():
    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--filepath",
        type=str,
        help="load in genome file, fasta or fna",
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        required=True,
        help="output folder",
    )
    parser.add_argument(
        "--size",
        type=int,
        help="split sub sequence length"
    )
    parser.add_argument(
        "--overlap",
        type=int,
        help="split sub sequence overlap length"
    )
    parser.add_argument(
        "--mode",
        type=str,
        help="process mode for nc or protein"
    )
    args = parser.parse_args()

    return args


def gen_split_overlap(seq, size, overlap):
    if size < 1 or overlap < 0:
        raise ValueError('size must be >= 1 and overlap >= 0')

    if len(seq) < size:
        yield seq
    else:
        for i in range(0, len(seq) - overlap, size - overlap):
            yield seq[i:i + size]


def gen_split_overlap_breakpoint(seq, size, overlap):
    if size < 1 or overlap < 0:
        raise ValueError('size must be >= 1 and overlap >= 0')

    if len(seq) < size:
        yield '0...' + str(len(seq))
    else:
        for i in range(0, len(seq) - overlap, size - overlap):
            yield '{}...{}'.format(i, i + size)


def gen_genome_overlap(df, size, overlap):

    sum_df = pd.DataFrame(
        columns=['vid', 'sequence', 'ori_seq', 'break_point'])
    for i in range(len(df)):
        vid = df.iloc[i]['vid']
        sequence = df.iloc[i]['sequence']

        sub_seqs_list = list(gen_split_overlap(sequence, size, overlap))
        sub_break_point = list(
            gen_split_overlap_breakpoint(sequence, size, overlap))

        sub_ids_list = [vid]*len(sub_seqs_list)
        sub_split_list = [i]*len(sub_seqs_list)
        new_df = pd.DataFrame(list(zip(sub_ids_list, sub_ids_list, sub_seqs_list, sub_break_point)), columns=[
                              'vid', 'ori_seq', 'sequence', 'break_point'])

        # 循环添加
        sum_df = sum_df.append(new_df, ignore_index=True)

    vid_uniq = sum_df['ori_seq'].drop_duplicates().reset_index(drop=True)
    dict_temp = dict(zip(vid_uniq.values, vid_uniq.index))
    sum_df['ori_seq'] = sum_df['ori_seq'].apply(lambda x: dict_temp[x])

    return sum_df


def split_protein(df, size, overlap):
    df['sequence'] = df.sequence.apply(
        lambda x: list(gen_split_overlap(x, size, overlap)))
    df = df.explode('sequence')
    df['part'] = df.groupby('vid').cumcount()+1
    df['part'] = df.apply(lambda row: row['vid'] + '_part_' + str(row['part']) if
                          row['length'] >= size else row['vid'], axis=1)
    df['length'] = df['sequence'].apply(lambda seq: len(seq))
    df['sequence'] = df.sequence.apply(lambda x: ' '.join(list(x)))
    return df


def Process_Rec(inputs, save_path, args):
    header_list, seq_list = [], []

    with open(inputs) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header_list.append(record.id)
            seq_list.append(record.seq._data)

    length_seq = [len(x) for x in seq_list]
    inputs_df = pd.DataFrame({'vid': header_list,
                              'sequence': seq_list,
                              'length': length_seq,
                              })
    inputs_df.to_csv(save_path, index=False, sep=',')

    return inputs_df


def Process_Pro_Rec(file_path, save_path, args):

    header_agg, seq_agg, length_agg, assembly_agg, ori_agg = [], [], [], [], []
    for filename in os.listdir(file_path):
        header_list, seq_list = [], []
        if filename.split('.')[-1] not in ['faa', 'fasta', 'fna']:
            continue

        vid = '.'.join(filename.split('.')[:-1])

        with open(file_path + filename) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if "_prot_" in record.id:
                    # Existing header format: >lcl|OQ715684.1_prot_WPV62740.1_2 [...]
                    header = record.id.split("_prot_", 1)[1]
                elif "_protein." in record.id:
                    # New header format: >AAD12264.2_protein.
                    header = record.id.split("_protein", 1)[0]
                else:
                    # Fallback: use the entire record.id or handle as needed
                    header = record.id

                header_list.append(header)
                seq_list.append(record.seq._data)

        length_list = [len(x) for x in seq_list]
        vid_list = [vid] * len(length_list)

        assembly_agg.extend(vid_list)
        header_agg.extend(header_list)
        seq_agg.extend(seq_list)
        length_agg.extend(length_list)

    inputs_df = pd.DataFrame({'accession': header_agg,
                              'vid': assembly_agg,
                              'ori_seq': assembly_agg,
                              'sequence': seq_agg,
                              'length': length_agg,
                              })

    list_vid_uniq = list(dict.fromkeys(assembly_agg))
    dict_temp = dict(zip(list_vid_uniq, range(len(list_vid_uniq))))
    inputs_df['ori_seq'] = inputs_df['ori_seq'].apply(lambda x: dict_temp[x])

    inputs_df.to_csv(save_path, index=False, sep=',')

    return inputs_df


def main():
    args = parse_args()

    if args.mode == 'Genome':

        whole_path = args.out_folder + 'processed_whole.csv'
        processed_path = args.out_folder + 'processed_input.csv'
        whole_df = Process_Rec(args.filepath, whole_path, args)
        whole_df = gen_genome_overlap(
            whole_df, size=args.size, overlap=args.overlap)
        whole_df.to_csv(processed_path, index=False)

    elif args.mode == 'Protein':

        whole_path = args.out_folder + 'processed_whole.csv'
        processed_path = args.out_folder + 'processed_input.csv'
        whole_df = Process_Pro_Rec(args.filepath, whole_path, args)
        whole_df = split_protein(
            whole_df, size=args.size, overlap=args.overlap)
        whole_df.to_csv(processed_path, index=False)
    else:
        raise ValueError("process mode must be 'Genome' or 'Protein'! ")


if __name__ == "__main__":
    main()
