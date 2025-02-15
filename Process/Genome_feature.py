import argparse
import pandas
import os
import pandas as pd
import warnings
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from subprocess import run
from tqdm import tqdm
from sklearn import preprocessing
import pickle
import re


def fuzzy_extract(dictionary, key, returnType='string'):
    """Extract a value from a dictionary which may not contain the key.
    Instead of failing when key is not present, returns either an empty
    string or 0, depending on type requested
    """
    try:
        result = dictionary[key]
    except KeyError:
        if returnType == 'string':
            result = ''
        elif returnType == 'intList':
            result = [0]
        else:
            raise NotImplementedError('Invalid returnType')
    return result


def check_cds(feature, extractedSeq, name=None, problems_file=None):
    """Check extracted CDS sequences for validity, returning an informative warning 
    messages before skipping/discarding invalid ones.
    Data on skipped sequences will be appended to the open file handle 'problems_file',
    if supplied (in csv format with columns: name, accession, protein_name, protein_id, length, problem_description).
    """

    accession = extractedSeq.id
    proteinName = ", ".join(fuzzy_extract(feature.qualifiers, "product"))
    proteinID = ", ".join(fuzzy_extract(feature.qualifiers, "protein_id"))
    baseMessage = "Sequence {0}: {1} (protein id {2})".format(
        accession, proteinName, proteinID)

    # Check if for causes of concern:
    wrongLength = len(extractedSeq) % 3 != 0
    truncatedStart = type(feature.location.start) == BeforePosition
    truncatedEnd = type(feature.location.end) == AfterPosition

    startPos = fuzzy_extract(feature.qualifiers, "codon_start", "intList")[0]
    shiftedStart = int(startPos) > 1

    # Respond to these
    if truncatedStart or truncatedEnd:
        problem = "truncation"
        message = "CDS skipped - truncation: {}".format(baseMessage)

    elif shiftedStart:
        problem = "translation start occurs before sequence start (truncation?)"
        message = "CDS skipped - shifted translation start found (codon {0}), assuming truncation: {1}".format(
            startPos, baseMessage)

    elif wrongLength:
        problem = "length not divisible by 3 (cause not identified)"
        message = "CDS skipped - length not a multiple of 3, but cause not identified: {}".format(
            baseMessage)

    else:
        problem = ""
        message = ""

    # Return:
    if len(message) != 0:
        warnings.warn(message, Warning)

        if problems_file is not None:
            outline = ",".join(
                [name, accession, proteinName, proteinID, str(len(extractedSeq)), problem])
            outline = "{0}\n".format(outline)
            _ = problems_file.write(outline)

        return None

    return extractedSeq


def rename_sequence(seqrecord, name):
    """Rename a seqrecord object so the specified name gets used when 
    writing to fasta
    """
    seqrecord.name = name
    seqrecord.id = name
    seqrecord.description = ''

    return seqrecord


def extract_and_name_cds(seqrecord, name, problems_file=None):
    """Extract noseg (valid) CDS sequences from a parsed genbank record and rename
    them to a common name.
    If given, 'problems_file' should be an open file handle to which data on problematic 
    sequences can be recorded.
    """
    codingSeqs = []

    for feature in seqrecord.features:
        if feature.type.lower() == "cds":
            seq = feature.extract(seqrecord)
            validSeq = check_cds(feature, seq, name, problems_file)

            if validSeq is not None:
                renamedSeq = rename_sequence(validSeq, name)
                codingSeqs.append(renamedSeq)

    return codingSeqs


def parse_args():
    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--filepath",
        type=str,
        help="load in genbank file, gb or gbff",
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
        "--CPB_path",
        type=str,
        help="path of CPB machine",
    )
    parser.add_argument(
        "--ref_cols",
        type=str,
        help="path of ref cols",
    )

    args = parser.parse_args()

    return args


# Step 1: Define a function to adjust df_b based on df_a
def adjust_accession(row, df_a):
    accession_b = row['accession']

    # Remove version number (if exists) from df_b
    base_accession_b = accession_b.split('.')[0]

    # Check if df_a has the base accession
    if base_accession_b in df_a['accession'].values:
        # If df_a has the base accession (without version), return the base accession (without version)
        return base_accession_b
    else:
        # If df_a has the full accession (with version), return the original accession from df_b
        return accession_b


def main():
    args = parse_args()

    # 1.load sequences to dict
    seqData = {}
    for filename in tqdm(os.listdir(args.filepath)):
        if filename[-2:] != 'gb' and filename[-4:] != 'gbff':
            continue
        acc = filename.split('.gb')[0]
        seqData[acc] = SeqIO.read(args.filepath + filename, format="genbank")

    processedIDs = []
    # 2.output file & problematic file
    outFastaName = args.out_folder + 'processed.fasta'
    outFasta = open(outFastaName, 'w')
    outWriter = SeqIO.FastaIO.FastaWriter(outFasta, wrap=None)
    outWriter.write_header()

    # 3. extract Fasta
    problemFile_name = args.out_folder + 'problematic.csv'
    problemFile = open(problemFile_name, 'w')
    problemFile.write(
        "name,accession,protein_name,protein_id,length,problem_description\n")

    for accession in tqdm(seqData.keys()):

        sequences = extract_and_name_cds(
            seqData[accession], accession, problems_file=problemFile)
        processedIDs.append(accession)

        if sequences is not None:
            outWriter.write_records(sequences)

    if problemFile is not None:
        problemFile.close()

    outWriter.write_footer()
    outFasta.close()

    # 4.run & create output file
    outData_name = args.out_folder + 'GF_accession.tsv'
    outData = open(outData_name, 'w')
    outDF = pandas.DataFrame()
    outDF["Name"] = seqData.keys()
    outDF["TaxID"] = "NA"
    outDF["Species"] = "NA"
    outDF["Family"] = "NA"

    outDF.to_csv(outData, sep='\t', index=False, header=False)
    outData.close()
    run(["java", "-jar", args.CPB_path, outFasta.name, outData.name], check=True)

    outputFile = "Genomic_features.tsv"
    os.rename(args.out_folder + 'processed_dat.txt',
              args.out_folder + outputFile)

    # 5. metadata and dataframe
    header_agg, assembly_agg, ori_agg = [], [], []
    for filename in os.listdir(args.assembly_path):
        header_list = []
        if filename.split('.')[-1] not in ['faa', 'fasta', 'fna']:
            continue

        vid = '.'.join(filename.split('.')[:-1])
        with open(args.assembly_path + filename) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header_list.append(record.id.split('_genome')[0])

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

    # 6. generate data df

    with open(args.ref_cols, 'rb') as fp:
        ref_cols = pickle.load(fp)

    min_max_scaler = preprocessing.MinMaxScaler()
    extract_df = pd.read_csv(args.out_folder + outputFile, sep='\t')
    # extract_df = extract_df.drop(['TaxID','Species','Good','Complete','Seqs','SeqLength','Codons','BadCodons','CodonPairs','Stops',], axis =1)
    extract_df = extract_df.replace(-9999.0, -10)

    extract_df = extract_df.rename(columns=lambda x: re.sub(
        '[()-]', '.', x))  # rename column names
    extract_df = extract_df.loc[:, ['SeqName'] + ref_cols]

    extract_df = extract_df.rename(columns={'SeqName': 'accession'})

    # input to follow the rule of extract

    inputs_df['accession'] = inputs_df.apply(
        adjust_accession, axis=1, df_a=extract_df)
    inputs_df = inputs_df.merge(extract_df, on='accession')
    inputs_df = inputs_df.drop(['accession'], axis=1)

    inputs_df.iloc[:, 2:] = min_max_scaler.fit_transform(inputs_df.iloc[:, 2:])
    inputs_df['genomic_features'] = inputs_df.iloc[:, 2:].apply(
        lambda x: '_'.join(x.astype(str)), axis=1)

    inputs_df = inputs_df[['vid', 'ori_seq',
                           'genomic_features']].reset_index(drop=True)
    inputs_df.to_csv(args.out_folder + 'processed_input.csv', index=False)


if __name__ == "__main__":
    main()
