import os
import pandas as pd
import argparse
import numpy as np
import pickle
from Bio import SeqIO


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
        help="feature folder",
    )
    parser.add_argument(
        "--assembly_path",
        type=str,
        required=True,
        help="assembly_path to extract metadata",
    )
    parser.add_argument(
        "--PPI_evidence",
        type=str,
        required=True,
        help="PPI metadata served as evidence",
    )
    parser.add_argument(
        "--PPI_evidence_col",
        type=str,
        required=True,
        help="PPI metadata served as evidence",
    )
    parser.add_argument(
        "--model_name_or_path",
        type=str,
        required=True,
        help="saved model path",
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

    args = parser.parse_args()

    return args
