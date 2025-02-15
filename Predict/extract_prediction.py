import pandas as pd
import argparse
import numpy as np
import pickle
import torch
import sys
import json

PERC = 0.4


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
        "--predict_folder",
        type=str,
        required=True,
        help="working folder",
    )
    parser.add_argument(
        "--feature_folder",
        type=str,
        required=True,
        help="working folder",
    )
    parser.add_argument(
        "--categories",
        type=str,
        help="pred class labels"
    )
    parser.add_argument(
        "--feature",
        type=str,
        help="feature model types"
    )
    parser.add_argument(
        "--cutoff",
        type=str,
        help="cutoff_path"
    )
    parser.add_argument(
        "--hierarchy_structure",
        type=str,
        default="binary"
    )
    parser.add_argument(
        "--agg_method",
        type=str,
        help="chunks prediction aggregation methods"
    )

    args = parser.parse_args()

    return args


def mean_result_probs(probs, oriseq, THRESHOLD, hierarchy_structure, cates):
    upper, lower = 1, 0
    mean_results = []
    mean_probs = []

    for i in range(oriseq.max()+1):
        sub_probs = probs[oriseq == i]
        if len(sub_probs) == 0:
            print("##################")
            print("There is no record for sequence {}".format(i))
            continue

        mean_prob = sub_probs.mean(axis=0)
        mean_result = np.where(pd.DataFrame(
            [mean_prob], columns=cates) > THRESHOLD, upper, lower)

        mean_results.append(np.squeeze(mean_result))
        mean_probs.append(mean_prob)

    mean_probs = violation_process(hierarchy_structure, np.stack(mean_probs))
    mean_results = violation_process(
        hierarchy_structure, np.stack(mean_results))

    return mean_results,  mean_probs


def max_result_probs(probs, oriseq, THRESHOLD, hierarchy_structure, cates):
    upper, lower = 1, 0
    max_results = []
    max_probs = []

    for i in range(oriseq.max()+1):
        sub_probs = probs[oriseq == i]
        if len(sub_probs) == 0:
            print("##################")
            print("There is no record for sequence {}".format(i))
            continue

        max_prob = sub_probs.max(axis=0)
        max_result = np.where(pd.DataFrame(
            [max_prob], columns=cates) > THRESHOLD, upper, lower)

        max_results.append(np.squeeze(max_result))
        max_probs.append(max_prob)
    max_probs = violation_process(hierarchy_structure, np.stack(max_probs))
    max_results = violation_process(hierarchy_structure, np.stack(max_results))

    return max_results,  max_probs


def topk_result_probs(probs, oriseq, THRESHOLD, hierarchy_structure, cates, perc):
    upper, lower = 1, 0
    mean_results = []
    mean_probs = []

    for i in range(oriseq.max()+1):
        sub_probs = probs[oriseq == i]
        if len(sub_probs) == 0:
            print("##################")
            print("There is no record for sequence {}".format(i))
            continue

        val, _ = torch.topk(torch.Tensor(sub_probs), int(
            np.ceil(perc*len(sub_probs))), axis=0)
        mean_prob = val.mean(axis=0).numpy()
        mean_result = np.where(pd.DataFrame(
            [mean_prob], columns=cates) > THRESHOLD, upper, lower)

        mean_results.append(np.squeeze(mean_result))
        mean_probs.append(mean_prob)

    mean_probs = violation_process(hierarchy_structure, np.stack(mean_probs))
    mean_results = violation_process(
        hierarchy_structure, np.stack(mean_results))

    return mean_results,  mean_probs


def main():
    args = parse_args()
    # save_file = args.predict_folder + 'pred_output.csv'
    probs = np.load(args.predict_folder + 'pred_probs.npy')
    
    if args.hierarchy_structure != 'binary':
        f = open(args.hierarchy_structure)
        taxa = json.load(f)
        hierarchy_structure = {int(k): int(v) for k, v in taxa.items()}
    else:
        hierarchy_structure = 'binary'

    probs = violation_process(hierarchy_structure, probs)

    with open(args.categories, 'rb') as fp:
        cates = pickle.load(fp)

    with open(args.cutoff, 'rb') as fp:
        cutoffs = pickle.load(fp)

    cutoff = cutoffs[args.feature]

    # special case for blast
    if args.feature == 'Homology':
        probs_df = pd.read_csv(args.predict_folder + 'whole_probs.csv')
        probs = violation_process(
            hierarchy_structure, probs_df.iloc[:, 1:].values)
        probs = pd.DataFrame(probs, columns=cates)
        result = np.where(probs > cutoff, 1, 0)
        result = violation_process(hierarchy_structure, result)
        pred_df = probs_df.vid
        names = probs_df.iloc[:, 1:].columns
        pred_df = pd.concat(
            [pred_df, pd.DataFrame(result, columns=names)], axis=1)
        pred_df.to_csv(args.predict_folder +
                       'whole_preds.csv', index=False, sep=',')

        agg_preds = pred_df.iloc[:, 1:]

        with open(args.predict_folder + "whole_results.txt", "w") as text_file:
            for i in range(len(pred_df)):
                text_file.write(str(pred_df.iloc[i, 0]) + '\n')
                text_file.writelines(
                    ', '.join(cates[np.where(agg_preds.iloc[i, :] == 1)].to_list()))
                text_file.write('\n\n')

        sys.exit(0)

    # sub sequences
    sub_df = pd.read_csv(args.feature_folder + 'processed_input.csv', sep=',')
    if 'sequence' in sub_df.columns:
        sub_df = sub_df.drop(['sequence'], axis=1)
    probs_df = pd.DataFrame(probs, columns=cates)

    sub_df = pd.concat([sub_df, probs_df], axis=1)
    sub_df.to_csv(args.predict_folder +
                  'probs_for_sub_sequences.csv', index=False, sep=',')

    if args.agg_method == 'mean':
        agg_preds, agg_probs = mean_result_probs(
            probs, sub_df.ori_seq, cutoff, hierarchy_structure, cates)
    elif args.agg_method == 'max':
        agg_preds, agg_probs = max_result_probs(
            probs, sub_df.ori_seq, cutoff, hierarchy_structure, cates)
    elif args.agg_method == 'topk':
        agg_preds, agg_probs = topk_result_probs(
            probs, sub_df.ori_seq, cutoff, hierarchy_structure, cates,  PERC)

    else:
        raise ValueError("agg_method must be 'mean', 'max' or 'topk'! ")

    agg_preds = pd.DataFrame(agg_preds, columns=cates)
    agg_probs = pd.DataFrame(agg_probs, columns=cates)

    sum_df = sub_df[['vid', 'ori_seq']
                    ].drop_duplicates().reset_index(drop=True)

    pd.concat([sum_df, agg_preds], axis=1).to_csv(
        args.predict_folder + 'whole_preds.csv', index=False, sep=',')
    pd.concat([sum_df, agg_probs], axis=1).to_csv(
        args.predict_folder + 'whole_probs.csv', index=False, sep=',')

    if args.feature != 'PPI':
        with open(args.predict_folder + "whole_results.txt", "w") as text_file:

            for i in range(len(sum_df)):
                text_file.write(str(sum_df.iloc[i, 0]) + '\n')
                text_file.writelines(
                    ', '.join(cates[np.where(agg_preds.iloc[i, :] == 1)].to_list()))
                text_file.write('\n\n')


if __name__ == "__main__":
    main()
