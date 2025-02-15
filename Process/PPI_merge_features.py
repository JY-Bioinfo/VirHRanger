import pandas as pd
import argparse
import pickle


def parse_args():
    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--out_folder",
        type=str,
        required=True,
        help="output folder",
    )
    parser.add_argument(
        "--Blast_out_file",
        type=str,
        required=True,
        help="Blast_out feature",
    )
    parser.add_argument(
        "--SLiM_folder",
        type=str,
        required=True,
        help="SLiM feature folder",
    )
    parser.add_argument(
        "--SLiM_all_motif",
        type=str,
        required=True
    )
    parser.add_argument(
        "--Human_protein_to_GO",
        type=str,
        required=True
    )
    parser.add_argument(
        "--HVPPI",
        type=str,
        required=True
    )
    parser.add_argument(
        "--Protein_rec_to_vid",
        type=str,
        required=True
    )

    args = parser.parse_args()

    return args

def add_columns(df, col_list):
    new_col = []
    for col in col_list:
        if col not in df.columns:
            new_col.append(col)

    df = df.join(pd.DataFrame(columns=new_col))
    return df[col_list].fillna(0).astype('float')

def main():
    args = parse_args()
    
    To_rec_DF = pd.read_csv(args.Protein_rec_to_vid)
    # part1. GO features
    mapping_to_go = pd.read_csv(args.Human_protein_to_GO)
    HVPPI = pd.read_csv(args.HVPPI)
    # blast 
    Header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore' ]

    PPI_blast = pd.read_csv(args.Blast_out_file, sep='\t', names =Header)
    PPI_blast['sseqid'] = PPI_blast['sseqid'].apply(
        lambda x: x.split('|')[1])
    PPI_blast['qseqid'] = PPI_blast['qseqid'].apply(
        lambda x: x.split("_prot_", 1)[1].split('.')[0])

    PPI_blast = PPI_blast.replace({'P19838-PRO_0000030311':'P19838'})
    PPI_blast = PPI_blast[['qseqid', 'sseqid', 'pident', 'bitscore']]
    PPI_blast['pident'] = PPI_blast['pident'].astype(int)
    PPI_blast['bitscore'] = PPI_blast['bitscore'].astype(int)
    # mapping
    PPI_blast = PPI_blast.merge(HVPPI).merge(mapping_to_go)
    blast_go = PPI_blast[[ 'qseqid', 'bitscore', 'GO_id']].groupby(['qseqid', 'GO_id']).max('bitscore').reset_index()  
    go_bit = blast_go.pivot(index = 'qseqid', columns = 'GO_id', values='bitscore')
    go_cols = mapping_to_go.GO_id.sort_values().drop_duplicates()
    go_bit = add_columns(go_bit, go_cols).reset_index()
    go_bit = go_bit.rename({'qseqid':'accession'}, axis=1)
    go_bit = To_rec_DF.merge(go_bit, how='left').fillna(0)
    
    
    # part2. SLiM features
    with open(args.SLiM_all_motif,"rb") as fp:   # Unpickling
         all_motif = pickle.load(fp)
            
    slim_df = pd.read_csv(args.SLiM_folder + '/slimprob.csv')
    slim_df = slim_df.reset_index(drop=True)
    slim_df = slim_df.drop(['SeqNum', 'RunID', 'Masking', 'RunTime', 'UPNum', 'N_UPC','E_UPC','p_UPC','pUnd_UPC','N_Seq','E_Seq','p_Seq','pUnd_Seq', 'Pattern'], axis=1)
    all_rec = slim_df.Dataset.drop_duplicates().to_list()
    slim_df['N_Occ_norm'] = slim_df['N_Occ']/slim_df['AANum']*1000
    slim_df['E_Occ_norm'] = slim_df['E_Occ']/slim_df['AANum']*1000
    slim_df['Occ_diff'] = slim_df['N_Occ'] - slim_df['E_Occ']
    slim_df['Occ_diff_norm'] = slim_df['N_Occ_norm'] - slim_df['E_Occ_norm']

    '''Occ features'''
    Occ_df = slim_df.pivot(index='Dataset', columns='Motif', values='N_Occ').fillna(0)
    Occ_df = Occ_df.reindex(index =all_rec, columns=all_motif, fill_value=0)
    new_columns = ['N_' + sub for sub in all_motif] 
    Occ_df = Occ_df.rename(columns=dict(zip(Occ_df.columns, new_columns))).reset_index().rename(columns={'Dataset':'accession'})


    Occ_df_norm = slim_df.pivot(index='Dataset', columns='Motif', values='N_Occ_norm').fillna(0)
    Occ_df_norm = Occ_df_norm.reindex(index =all_rec, columns=all_motif, fill_value=0)
    new_columns = ['N_' + sub +'norm' for sub in all_motif] 
    Occ_df_norm = Occ_df_norm.rename(columns=dict(zip(Occ_df_norm.columns, new_columns))).reset_index().rename(columns={'Dataset':'accession'})


    '''Diff features'''
    Diff_df = slim_df.pivot(index='Dataset', columns='Motif', values='Occ_diff').fillna(0)
    Diff_df = Diff_df.reindex(index =all_rec, columns=all_motif, fill_value=0)
    new_columns = ['diff_' + sub for sub in all_motif] 
    Diff_df = Diff_df.rename(columns=dict(zip(Diff_df.columns, new_columns))).reset_index().rename(columns={'Dataset':'accession'})


    Diff_df_norm = slim_df.pivot(index='Dataset', columns='Motif', values='Occ_diff_norm').fillna(0)
    Diff_df_norm = Diff_df_norm.reindex(index =all_rec, columns=all_motif, fill_value=0)
    new_columns = ['diff_' + sub + 'norm' for sub in all_motif] 
    Diff_df_norm = Diff_df_norm.rename(columns=dict(zip(Diff_df_norm.columns, new_columns))).reset_index().rename(columns={'Dataset':'accession'})

    '''length feature'''
    length_df = slim_df[['Dataset','AANum']].drop_duplicates().reset_index(drop=True).rename(columns={'Dataset':'accession'})
    SLiM_df = length_df.merge(Occ_df_norm).merge(Diff_df_norm)

    
    SLiM_df = To_rec_DF[['vid','ori_seq','accession']].merge(SLiM_df)
    
    ### merge
    merge_DF = SLiM_df.merge(go_bit).copy()
    
    merge_DF['PPI_features'] = merge_DF.iloc[:,3:].apply(lambda x: '_'.join(x.astype(int).astype(str)), axis =1)
    merge_DF = merge_DF[['vid','ori_seq','accession','PPI_features']].reset_index(drop=True)
    merge_DF.to_csv(args.out_folder  + 'processed_input.csv', index=False)
    
if __name__ == "__main__":
    main()