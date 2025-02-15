# coding=utf-8
# Copyright 2021 The HuggingFace Inc. team. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Finetuning a 🤗 Transformers model for sequence classification on viruses."""
from transformers import bioinfo_compute_metrics as compute_metrics
from transformers import (
    # modified bert models
    BigBirdConfig,
    BigBirdTokenizer,
    BigBirdForSequenceClassification,
    BigBirdForMultiLabel_SC,
    BigBirdForMHC,
    BertForMHC,
    BertConfig,
    BertMLPCLS,
    BertForSequenceClassification,
    BertTokenizer,
    AlbertForMHC,
    AlbertConfig,
    AlbertTokenizer,
    PreTrainedTokenizerFast,
    PreTrainedTokenizer,
    AdamW,
    DataCollatorWithPadding,
    default_data_collator,
    SchedulerType,
    get_scheduler,
    set_seed,
    AutoConfig
)
import argparse
import logging
import math
import os
import random
import numpy as np
import datasets
import torch
from datasets import load_dataset, load_metric
from torch.utils.data.dataloader import DataLoader
from tqdm.auto import tqdm

# checkpoint watermark
import glob
import re
import shutil
from datetime import datetime
import json

# import tensorboard
from torch.utils.tensorboard import SummaryWriter
from accelerate import Accelerator, DeepSpeedPlugin

import transformers
print(transformers.__version__)


# added metrics for downstream tasks
output_modes = {
    "PRIOR": "boolean_cls",
    "VHP": "multi_cls",
    "Corona": "multi_cls",
    "VIV_mhc": "mhc_label",
    "VIV_mhc_pos": "mhc_label",
    "VIV_Pro": "mhc_label",
    "VIV_tropism_tr": "multi_label",
    "VIV_tropism_tp": "multi_label",
    "VIV_root": "boolean_cls",
    "VIV_binary": "boolean_cls",
}

MODEL_CLASSES = {
    "BigBird_cls": (BigBirdConfig, BigBirdForSequenceClassification, BigBirdTokenizer),
    "BigBirdMulti": (BigBirdConfig, BigBirdForMultiLabel_SC, BigBirdTokenizer),
    "BigBirdMHC": (BigBirdConfig, BigBirdForMHC, BigBirdTokenizer),
    "BertMHC": (BertConfig, BertForMHC, BertTokenizer),
    "Bert_cls": (BertConfig, BertForSequenceClassification, BertTokenizer),
    "Bert_MLP": (BertConfig, BertMLPCLS, BertTokenizer),
    "Bert_Flat": (BertConfig, BertMLPCLS, BertTokenizer),
    "AlbertMHC": (AlbertConfig, AlbertForMHC, BertTokenizer),
}

task_to_keys = {
    "PRIOR": ("sequence", None),
    "VHDB_cls": ("sequence", None),
    "VHDB_multi": ("sequence", None),
    "VHP": ("sequence", None),
    "Corona": ("sequence", None),
    "VIV_mhc": ("sequence", None),
    "VIV_mhc_pos": ("sequence", None),
    "VIV_Pro": ("sequence", None),
    "VIV_tropism_tr": ("sequence", None),
    "VIV_tropism_tp": ("sequence", None),
    "VIV_root": ("sequence", None),
    "VIV_binary": ("sequence", None),
}

CUDA_VISIBLE_DEVICES = 0, 1, 2, 3

# write logging info to path


def get_logger(args, verbosity=1, name=None):
    filename = args.logger_path
    level_dict = {0: logging.DEBUG, 1: logging.INFO, 2: logging.WARNING}
    formatter = logging.Formatter(
        "[%(asctime)s][%(filename)s][line:%(lineno)d][%(levelname)s] %(message)s"
    )
    logger = logging.getLogger(name)
    logger.setLevel(level_dict[verbosity])

    # output to file
    fh = logging.FileHandler(filename, "w")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # output to stream
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger


def parse_args():
    parser = argparse.ArgumentParser(
        description="Finetune a transformers model on a text classification task")

# 1. data & files
    parser.add_argument(
        "--task_name",
        type=str,
        default=None,
        help="The name of the virus prediction task to train on.",
        choices=list(task_to_keys.keys()),
    )
    parser.add_argument(
        "--train_file", type=str, default=None, help="A csv or a json file containing the training data."
    )
    parser.add_argument(
        "--validation_file", type=str, default=None, help="A csv or a json file containing the validation data."
    )
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Where to store the final model.")
    parser.add_argument(
        "--max_length",
        type=int,
        default=512,
        help=(
            "The maximum total input sequence length after tokenization. Sequences longer than this will be truncated,"
            " sequences shorter will be padded if `--pad_to_max_lengh` is passed."
        ),
    )
    parser.add_argument(
        "--pad_to_max_length",
        action="store_true",
        help="If passed, pad all samples to `max_length`. Otherwise, dynamic padding is used.",
    )
    parser.add_argument(
        "--model_name_or_path",
        type=str,
        help="Path to pretrained model or model identifier from huggingface.co/models.",
        required=True,
    )
    parser.add_argument(
        "--use_slow_tokenizer",
        action="store_true",
        help="If passed, will use a slow tokenizer (not backed by the 🤗 Tokenizers library).",
    )
    parser.add_argument(
        "--logger_path",
        default='Logger/default_log_file.log',
        type=str,
        help="Save the logger output to files",
    )
    parser.add_argument("--logging_steps", type=int,
                        default=0, help="Log every X updates steps.")
    parser.add_argument("--logging_per_epoch", type=int,
                        default=0, help="Log per epoch.")
    parser.add_argument("--save_every_log", type=int,
                        default=1, help="Save checkpoint every X logs.")
    parser.add_argument(
        "--save_total_limit",
        type=int,
        default=None,
        help="Limit the total amount of checkpoints, delete the older checkpoints in the output_dir, does not delete by default",)
    parser.add_argument("--cv_num", default=0, type=str,
                        help="cross validation fold")
    parser.add_argument(
        "--hierarchy_structure",
        type=str,
        help="Give hierarchy structure path"
    )


# 2. training hyperparameters
    parser.add_argument(
        "--per_device_train_batch_size",
        type=int,
        default=8,
        help="Batch size (per device) for the training dataloader.",
    )
    parser.add_argument(
        "--per_device_eval_batch_size",
        type=int,
        default=8,
        help="Batch size (per device) for the evaluation dataloader.",
    )
    parser.add_argument(
        "--learning_rate",
        type=float,
        default=5e-5,
        help="Initial learning rate (after the potential warmup period) to use.",
    )
    parser.add_argument(
        "--learning_rate_classifier",
        type=float,
        default=1e-3,
        help="Initial learning rate for classifier.",
    )
    parser.add_argument("--weight_decay", type=float,
                        default=0.0, help="Weight decay to use.")
    parser.add_argument("--num_train_epochs", type=int, default=3,
                        help="Total number of training epochs to perform.")
    parser.add_argument(
        "--max_train_steps",
        type=int,
        default=None,
        help="Total number of training steps to perform. If provided, overrides num_train_epochs.",
    )
    parser.add_argument(
        "--gradient_accumulation_steps",
        type=int,
        default=1,
        help="Number of updates steps to accumulate before performing a backward/update pass.",
    )
    parser.add_argument(
        "--lr_scheduler_type",
        type=SchedulerType,
        default="linear",
        help="The scheduler type to use.",
        choices=["linear", "cosine", "cosine_with_restarts",
                 "polynomial", "constant", "constant_with_warmup"],
    )
    parser.add_argument(
        "--num_warmup_steps", type=int, default=0, help="Number of steps for the warmup in the lr scheduler."
    )
    parser.add_argument("--warmup_percent", default=0, type=float,
                        help="Linear warmup over warmup_percent*total_steps.")
    parser.add_argument(
        "--early_stop", default=-1, type=int, help="set this to a positive integet if you want to perfrom early stop. The model will stop \
                                                    if the val loss keep increasing early_stop times",
    )
    parser.add_argument(
        "--early_stop_epoch_manual", default=10000, type=int, help="set this to a defined epoch number to perform early stop without validation. we get this number first as a hyperparameter and use that for the whole training datas"
    )
    parser.add_argument("--beta1", default=0.9, type=float,
                        help="Beta1 for Adam optimizer.")
    parser.add_argument("--beta2", default=0.999, type=float,
                        help="Beta2 for Adam optimizer.")
    # added for class balance
    parser.add_argument("--weighting", default=None,
                        type=str, help="Class Balancing methods.")
    parser.add_argument("--CLS_dropout", default=0.0,
                        type=float, help="Dropout rate of classifier layer.")
    parser.add_argument(
        "--multi_threshold",
        type=float,
        default=0.5,
        help="Give threshold for sigmoid"
    )
    parser.add_argument(
        "--coe_l",
        type=float,
        default=0.1,
        help="penalty for hierarchical loss"
    )
    parser.add_argument(
        "--MLP_width",
        type=int,
        default=768,
        help="MLP layer width"
    )

# 3. model and config
    parser.add_argument(
        "--model_type_limit",
        type=str,
        required=True,
        help="Model type to use if training from scratch.",
    )
    parser.add_argument(
        "--config_name",
        type=str,
        default=None,
        help="Pretrained config name or path if not the same as model_name",
    )
    parser.add_argument(
        "--tokenizer_name",
        type=str,
        default=None,
        help="Pretrained tokenizer name or path if not the same as model_name",
    )
    parser.add_argument(
        "--customize_tokenizer", action="store_true", help="Overwrite the tokenizer used"
    )
    parser.add_argument(
        "--unfrozen_layers",
        type=int,
        default=-1,
        help="use para to freeze BERT layers, default -1 not freeze"
    )

    # added for focal loss/ penalty loss
    parser.add_argument("--loss_type", default="CELoss",
                        type=str, help="Type of loss used")

# 4. working modes
    parser.add_argument("--do_train", action="store_true",
                        help="Whether to run training.")
    parser.add_argument("--do_predict", action="store_true",
                        help="Whether to do prediction on the given dataset.")
    parser.add_argument("--do_predict_no_label", action="store_true",
                        help="Whether to do prediction on the given dataset with no label information.")
    parser.add_argument("--do_visualize", action="store_true",
                        help="Whether to calculate attention score.")
    parser.add_argument("--do_embedding", action="store_true",
                        help="Whether to output embedding.")


# 5. computational parameters
    parser.add_argument("--seed", type=int, default=None,
                        help="A seed for reproducible training.")
    parser.add_argument(
        "--preprocessing_num_workers",
        type=int,
        default=None,
        help="The number of processes to use for the preprocessing.",
    )
    parser.add_argument("--adam_epsilon", default=1e-8,
                        type=float, help="Epsilon for Adam optimizer.")
    parser.add_argument("--fp16", default=False,
                        type=bool, help="FP 16 accelerate")

    args = parser.parse_args()

    # Sanity checks
    if args.task_name is None and args.train_file is None and args.validation_file is None:
        raise ValueError(
            "Need either a task name or a training/validation file.")
    else:
        if args.train_file is not None:
            extension = args.train_file.split(".")[-1]
            assert extension in [
                "csv", "json"], "`train_file` should be a csv or a json file."
        if args.validation_file is not None:
            extension = args.validation_file.split(".")[-1]
            assert extension in [
                "csv", "json"], "`validation_file` should be a csv or a json file."

    if args.output_dir is not None:
        os.makedirs(args.output_dir, exist_ok=True)

    return args


def _sorted_checkpoints(args, checkpoint_prefix="checkpoint", use_mtime=False):
    ordering_and_checkpoint_path = []

    glob_checkpoints = glob.glob(os.path.join(
        args.output_dir, "{}-*".format(checkpoint_prefix)))

    for path in glob_checkpoints:
        if use_mtime:
            ordering_and_checkpoint_path.append((os.path.getmtime(path), path))
        else:
            regex_match = re.match(
                ".*{}-([0-9]+)".format(checkpoint_prefix), path)
            if regex_match and regex_match.groups():
                ordering_and_checkpoint_path.append(
                    (int(regex_match.groups()[0]), path))

    checkpoints_sorted = sorted(ordering_and_checkpoint_path)
    checkpoints_sorted = [checkpoint[1] for checkpoint in checkpoints_sorted]
    return checkpoints_sorted

# save model


def _rotate_checkpoints(args, checkpoint_prefix="checkpoint", use_mtime=False):
    if not args.save_total_limit:
        return
    if args.save_total_limit <= 0:
        return

    # Check if we should delete older checkpoint(s)
    checkpoints_sorted = _sorted_checkpoints(
        args, checkpoint_prefix, use_mtime)
    if len(checkpoints_sorted) <= args.save_total_limit:
        return

    number_of_checkpoints_to_delete = max(
        0, len(checkpoints_sorted) - args.save_total_limit)
    checkpoints_to_be_deleted = checkpoints_sorted[:number_of_checkpoints_to_delete]
    for checkpoint in checkpoints_to_be_deleted:
        logger.info(
            "Deleting older checkpoint [{}] due to args.save_total_limit".format(checkpoint))
        shutil.rmtree(checkpoint)


def np_relu(x):
    return (x > 0) * x


def violation_process(hierarchy_structure, probs):

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


def evaluate(args, model, eval_dataloader, completed_steps, tb_time, log_info):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running evaluation {} *****")
    logger.info("  Num examples = %d", len(eval_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    eval_loss, avg_steps = 0, 0
    logits = None
    model.eval()

    # sanity check, avoid nan loss
    output_eval_nan_file = os.path.join(args.output_dir, "eval_nan_loss.txt")

    for eval_step, eval_batch in enumerate(tqdm(eval_dataloader)):
        # save GPU memory in evaluate
        with torch.no_grad():
            outputs = model(**eval_batch)
        if logits is None:
            logits = outputs.logits.detach().cpu().numpy()
            reference = eval_batch["labels"].detach().cpu().numpy()
        else:

            logits = np.append(
                logits, outputs.logits.detach().cpu().numpy(), axis=0)
            reference = np.append(
                reference, eval_batch["labels"].detach().cpu().numpy(), axis=0)

        # sanity check, avoid nan for logits
        loss_add = outputs.loss.item()
        if loss_add == loss_add:
            eval_loss += loss_add
            avg_steps += 1

    # avg
    eval_loss = eval_loss / avg_steps

    # add eval_loss_result
    output_eval_loss_file = os.path.join(args.output_dir, "eval_losses.txt")

    d3 = tb_time.strftime("%h%d_%H.%M")
    board_writer = SummaryWriter(log_dir='runs/{}/{}_{}_LR_{}_clsLR_{}_Batch_{}_accum_{}_CV_{}_{}_seed{}'.format(args.task_name, d3, args.max_length, args.learning_rate,
                                 args.learning_rate_classifier, args.per_device_train_batch_size, args.gradient_accumulation_steps, args.cv_num, args.weighting, args.seed))
    board_writer.add_scalar('Loss/eval', eval_loss, completed_steps)

    with open(output_eval_loss_file, "a") as writer:
        writer.write(str(eval_loss) + "\n")

    if args.output_mode == "boolean_cls":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    elif args.output_mode == "multi_cls":
        probs = softmax(torch.tensor(logits, dtype=torch.float32)).numpy()
        preds = np.argmax(logits, axis=1)

    elif args.output_mode == "mhc_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()

        # post_processing of hierarchical violation
        f = open(args.hierarchy_structure, )
        taxa = json.load(f)
        taxa = {int(k): int(v) for k, v in taxa.items()}
        probs = violation_process(taxa, probs)

        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)
    elif args.output_mode == "multi_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    else:
        print("Wrong output mode!")

    if args.model_type_limit == 'BertMHC' or args.model_type_limit == 'Bert_Flat' or args.model_type_limit == 'Bert_MLP':
        result = compute_metrics('mhc', preds, reference, probs)
    else:
        result = compute_metrics('binary', preds, reference, probs)
    # define eval output file name
    output_eval_file = os.path.join(args.output_dir, "eval_results.csv")
    output_eval_file_cm = os.path.join(args.output_dir, "eval_cm.npy")
    output_eval_file_probs = os.path.join(
        args.output_dir, "eval_probs_{}.npy".format(log_info))
    with open(output_eval_file, "a") as writer:
        logger.info("***** Eval results*****")
        for key in sorted(result.keys()):
            if key != 'confusion_matrix':
                logger.info("  %s = %s", key, str(result[key]))
                if args.output_mode == "mhc_label" or args.output_mode == "multi_label":
                    writer.write(str(result[key]) + ",")
                else:
                    writer.write(str(round(result[key], 3)) + ",")
        writer.write("\n")

    logger.info("Confusion Matrix : \n {}".format(result['confusion_matrix']))
    np.save(output_eval_file_cm, result['confusion_matrix'])
    np.save(output_eval_file_probs, probs)
    # add new values to dict

    # tensorboard output

    if args.model_type_limit == 'BertMHC' or args.model_type_limit == 'Bert_Flat' or args.model_type_limit == 'Bert_MLP':
        board_writer.add_pr_curve(
            'pr_curve', reference, probs, completed_steps)
        board_writer.add_scalar('Metrics/precision_micro',
                                result['precision_mi'], completed_steps)
        board_writer.add_scalar('Metrics/recall_micro',
                                result['recall_mi'], completed_steps)
        board_writer.add_scalar('Metrics/precision_macro',
                                result['precision_ma'], completed_steps)
        board_writer.add_scalar('Metrics/recall_macro',
                                result['recall_ma'], completed_steps)

        board_writer.add_scalar(
            'Metrics/f1_micro', result['f1_micro'], completed_steps)
        board_writer.add_scalar(
            'Metrics/f1_macro', result['f1_macro'], completed_steps)

        board_writer.add_scalar('Metrics/auc_micro',
                                result['auc_micro'], completed_steps)
        board_writer.add_scalar('Metrics/auc_macro',
                                result['auc_macro'], completed_steps)
        board_writer.add_scalar('Metrics/aupr_micro',
                                result['aupr_micro'], completed_steps)
        board_writer.add_scalar('Metrics/aupr_macro',
                                result['aupr_macro'], completed_steps)
        board_writer.add_scalar('Metrics/hamming_loss',
                                result['hamming_loss'], completed_steps)
        board_writer.add_scalar('Metrics/precision_micro',
                                result['precision_mi'], completed_steps)

    else:
        board_writer.add_scalar('Metrics/acc', result['acc'], completed_steps)
        board_writer.add_scalar('Metrics/auc', result['auc'], completed_steps)
        board_writer.add_scalar('Metrics/precision',
                                result['precision'], completed_steps)
        board_writer.add_scalar(
            'Metrics/recall', result['recall'], completed_steps)
        board_writer.add_scalar('Metrics/f1', result['f1'], completed_steps)

    board_writer.close()

    return result, eval_loss


def predict(args, model, pred_dataloader):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running predict {} *****")
    logger.info("  Num examples = %d", len(pred_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    pred_loss, avg_steps = 0, 0
    logits = None
    embeddings = None

    model.eval()

    # sanity check, avoid nan loss
    output_pred_nan_file = os.path.join(args.output_dir, "pred_nan_loss.txt")

    for pred_step, pred_batch in enumerate(tqdm(pred_dataloader)):
        # save GPU memory in predict
        with torch.no_grad():
            outputs = model(**pred_batch)
        if logits is None:
            logits = outputs.logits.detach().cpu().numpy()
            reference = pred_batch["labels"].detach().cpu().numpy()
            embeddings = outputs.embedding.detach().cpu().numpy()
        else:

            logits = np.append(
                logits, outputs.logits.detach().cpu().numpy(), axis=0)
            reference = np.append(
                reference, pred_batch["labels"].detach().cpu().numpy(), axis=0)
            embeddings = np.append(
                embeddings, outputs.embedding.detach().cpu().numpy(), axis=0)

        # sanity check, avoid nan for logits
        loss_add = outputs.loss.item()
        if loss_add == loss_add:
            pred_loss += loss_add
            avg_steps += 1

    # avg
    pred_loss = pred_loss / avg_steps

    # add pred_loss_result
    output_pred_loss_file = os.path.join(args.output_dir, "pred_losses.txt")

    with open(output_pred_loss_file, "a") as writer:
        writer.write(str(pred_loss) + "\n")

    if args.output_mode == "boolean_cls":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    elif args.output_mode == "multi_cls":
        probs = softmax(torch.tensor(logits, dtype=torch.float32)).numpy()
        preds = np.argmax(logits, axis=1)

    elif args.output_mode == "mhc_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()

        # post_processing of hierarchical violation
        f = open(args.hierarchy_structure, )
        taxa = json.load(f)
        taxa = {int(k): int(v) for k, v in taxa.items()}
        probs = violation_process(taxa, probs)

        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)
    elif args.output_mode == "multi_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    else:
        print("Wrong output mode!")

    # post-processing of probs

    if args.model_type_limit == 'BertMHC' or args.model_type_limit == 'Bert_Flat' or args.model_type_limit == 'Bert_MLP':
        result = compute_metrics('mhc', preds, reference, probs)
    else:
        result = compute_metrics('binary', preds, reference, probs)

    output_embedding = os.path.join(args.output_dir, "embedding.npy")
    output_reference = os.path.join(args.output_dir, "refs.npy")
    # define pred output file name
    output_pred_file = os.path.join(args.output_dir, "pred_results.csv")
    output_pred_file_cm = os.path.join(args.output_dir, "pred_cm.npy")
    output_pred_file_probs = os.path.join(args.output_dir, "pred_probs.npy")
    with open(output_pred_file, "a") as writer:
        logger.info("***** Pred results*****")
        for key in sorted(result.keys()):
            if key != 'confusion_matrix':
                logger.info("  %s = %s", key, str(result[key]))
                if args.output_mode == "mhc_label" or args.output_mode == "multi_label":
                    writer.write(str(result[key]) + ",")
                else:
                    writer.write(str(round(result[key], 3)) + ",")
        writer.write("\n")

    logger.info("Confusion Matrix : \n {}".format(result['confusion_matrix']))
    np.save(output_pred_file_cm, result['confusion_matrix'])
    np.save(output_pred_file_probs, probs)
    # save embeddings

    np.save(output_embedding, embeddings)
    np.save(output_reference, reference)

    return result


def predict_no_label(args, model, pred_dataloader):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running predict {} *****")
    logger.info("  Num examples = %d", len(pred_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    logits = None
    embeddings = None

    model.eval()

    for pred_step, pred_batch in enumerate(tqdm(pred_dataloader)):
        # save GPU memory
        with torch.no_grad():
            outputs = model(**pred_batch)
        if logits is None:
            logits = outputs.logits.detach().cpu().numpy()
            embeddings = outputs.embedding.detach().cpu().numpy()
        else:
            logits = np.append(
                logits, outputs.logits.detach().cpu().numpy(), axis=0)
            embeddings = np.append(
                embeddings, outputs.embedding.detach().cpu().numpy(), axis=0)

    if args.output_mode == "boolean_cls":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    elif args.output_mode == "multi_cls":
        probs = softmax(torch.tensor(logits, dtype=torch.float32)).numpy()
        preds = np.argmax(logits, axis=1)

    elif args.output_mode == "mhc_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        # post_processing of hierarchical violation
        f = open(args.hierarchy_structure, )
        taxa = json.load(f)
        taxa = {int(k): int(v) for k, v in taxa.items()}
        probs = violation_process(taxa, probs)
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)
    elif args.output_mode == "multi_label":
        probs = torch.sigmoid(torch.tensor(
            logits, dtype=torch.float32)).numpy()
        upper, lower = 1, 0
        THRESHOLD = args.multi_threshold
        preds = np.where(probs > THRESHOLD, upper, lower)

    else:
        print("Wrong output mode!")

    output_embedding = os.path.join(args.output_dir, "embedding.npy")
    # define pred output file name
    output_pred_file_probs = os.path.join(args.output_dir, "pred_probs.npy")
    np.save(output_pred_file_probs, probs)

    # save embedding
    np.save(output_embedding, embeddings)

    return None


def give_embedding(args, model, emb_dataloader):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running embedding {} *****")
    logger.info("  Num examples = %d", len(emb_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    logits = None
    embeddings = None

    model.eval()

    for emb_step, emb_batch in enumerate(tqdm(emb_dataloader)):
        # save GPU memory
        with torch.no_grad():
            outputs = model(**emb_batch)

        if logits is None:
            logits = outputs.logits.detach().cpu().numpy()
            embeddings = outputs.embedding.detach().cpu().numpy()
        else:
            logits = np.append(
                logits, outputs.logits.detach().cpu().numpy(), axis=0)
            embeddings = np.append(
                embeddings, outputs.embedding.detach().cpu().numpy(), axis=0)

    output_embedding = os.path.join(args.output_dir, "embedding.npy")
    np.save(output_embedding, embeddings)

    return None


def visualize(args, model, pred_dataloader):

    return None


def main():
    args = parse_args()
    tb_time = datetime.now()
    args.output_mode = output_modes[args.task_name]
    global logger
    logger = get_logger(args)

    # Initialize the accelerator. We will let the accelerator handle device placement for us in this example.
    # add deepspeed
    deepspeed_plugin = DeepSpeedPlugin(
        zero_stage=2, gradient_accumulation_steps=args.gradient_accumulation_steps)

    accelerator = Accelerator(
        fp16=args.fp16, deepspeed_plugin=deepspeed_plugin)

    # Make one log on every process with the configuration for debugging.
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
        level=logging.INFO,
    )
    logger.info(accelerator.state)

    # Setup logging, we only want one process per machine to log things on the screen.
    # accelerator.is_local_main_process is only True for one process per machine.
    logger.setLevel(
        logging.INFO if accelerator.is_local_main_process else logging.ERROR)
    if accelerator.is_local_main_process:
        datasets.utils.logging.set_verbosity_warning()
        transformers.utils.logging.set_verbosity_info()
    else:
        datasets.utils.logging.set_verbosity_error()
        transformers.utils.logging.set_verbosity_error()

    # If passed along, set the training seed now.
    if args.seed is not None:
        set_seed(args.seed)

    # define tasks
    if args.task_name is not None:

        data_files = {}
        if args.train_file is not None:
            data_files["train"] = args.train_file
        if args.validation_file is not None:
            data_files["validation"] = args.validation_file
        extension = (
            args.train_file if args.train_file is not None else args.validation_file).split(".")[-1]
        raw_datasets = load_dataset(extension, data_files=data_files)

    _, model_class, tokenizer_class = MODEL_CLASSES[args.model_type_limit]

    # use customized tokenizer
    if args.customize_tokenizer:
        tokenizer = PreTrainedTokenizerFast(tokenizer_file=args.tokenizer_name)
        tokenizer.add_special_tokens({"unk_token": '<unk>'})
        tokenizer.add_special_tokens({'pad_token': '<pad>'})
        tokenizer.add_special_tokens({"bos_token": '<s>'})
        tokenizer.add_special_tokens({"eos_token": '</s>'})
        tokenizer.add_special_tokens({"sep_token": '[SEP]'})
        tokenizer.add_special_tokens({"cls_token": '[CLS]'})
        tokenizer.add_special_tokens({'mask_token': '[MASK]'})
        print(tokenizer.vocab_size)

    else:
        if args.tokenizer_name:
            tokenizer = tokenizer_class.from_pretrained(
                args.tokenizer_name, use_fast=True)
        elif args.model_name_or_path:
            tokenizer = tokenizer_class.from_pretrained(
                args.model_name_or_path, use_fast=True)
        else:
            raise ValueError(
                "You are instantiating a new tokenizer from scratch. This is not supported by this script."
                "You can do it from another script, save it, and load it from here, using --tokenizer_name."
            )

    # if False, use dynamic padding; if True, use static
    padding = "max_length" if args.pad_to_max_length else True
    # padding = "max_length" if args.pad_to_max_length else False

    # Preprocessing the datasets
    if args.task_name is not None:
        sentence1_key, sentence2_key = task_to_keys[args.task_name]

    # added to translate mhc label structure and sample-wised weight
    def Translate_multi_label(l):
        label = l.split('_')
        label = [int(i) for i in label]
        return label

    def Translate_multi_weight(w):
        weight = w.split('_')
        weight = [float(i) for i in weight]
        return weight

    def preprocess_function(examples):
        # Tokenize the texts
        texts = (
            (examples['sequence'],)
        )
        result = tokenizer(*texts, padding=padding,
                           max_length=args.max_length, truncation=True)

        if "label" in examples:
            if label_to_id is not None:
                # Map labels to IDs (not necessary for GLUE tasks)
                result["labels"] = [label_to_id[l] for l in examples["label"]]
            else:
                # In all cases, rename the column to labels because the model will expect that.
                result["labels"] = examples["label"]

        if "weight" in examples:
            # translate labels to list
            result["weights"] = examples["weight"]

        return result

    def preprocess_multilabel_function(examples):
        # Tokenize the texts
        texts = (
            (examples['sequence'],)
        )
        result = tokenizer(*texts, padding=padding,
                           max_length=args.max_length, truncation=True)

        if "label" in examples:
            # translate labels to list
            result["labels"] = list(
                map(lambda l: Translate_multi_label(l), examples['label']))

        if "weight" in examples:
            # translate labels to list
            result["weights"] = list(
                map(lambda w: Translate_multi_weight(w), examples['weight']))
        return result

    # Part 1: Three Load data modes

    if args.do_train:
        if args.output_mode == "boolean_cls":
            num_labels = 1
            label_to_id = {0: 0, 1: 1}
            processed_datasets = raw_datasets.map(
                preprocess_function, batched=True, remove_columns=raw_datasets["train"].column_names,
            )
        elif args.output_mode == "multi_cls":

            label_list = raw_datasets["train"].unique("label")
            label_list.sort()  # Let's sort it for determinism
            num_labels = len(label_list)

            label_to_id = {v: i for i, v in enumerate(label_list)}

            processed_datasets = raw_datasets.map(
                preprocess_function, batched=True, remove_columns=raw_datasets["train"].column_names,
            )
        else:
            new_labels = list(map(lambda l: Translate_multi_label(
                l), raw_datasets['train']['label']))
            num_labels = len(new_labels[0])
            processed_datasets = raw_datasets.map(
                preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["train"].column_names,
            )

        train_dataset = processed_datasets["train"]
        eval_dataset = processed_datasets["validation"]

    elif args.do_predict:
        if args.output_mode == "boolean_cls":
            num_labels = 1
            label_to_id = {0: 0, 1: 1}
            processed_datasets = raw_datasets.map(
                preprocess_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
            )
        elif args.output_mode == "multi_cls":

            label_list = raw_datasets["validation"].unique("label")
            label_list.sort()  # Let's sort it for determinism
            num_labels = len(label_list)

            label_to_id = {v: i for i, v in enumerate(label_list)}

            processed_datasets = raw_datasets.map(
                preprocess_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
            )
        else:
            new_labels = list(map(lambda l: Translate_multi_label(
                l), raw_datasets['validation']['label']))
            num_labels = len(new_labels[0])
            processed_datasets = raw_datasets.map(
                preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
            )

        eval_dataset = processed_datasets["validation"]

    else:
        # manually set num of labels
        num_labels = 79
        if args.output_mode != 'multi_label':

            processed_datasets = raw_datasets.map(
                preprocess_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
            )
        else:
            processed_datasets = raw_datasets.map(
                preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
            )

        eval_dataset = processed_datasets["validation"]

    # All config
    config = BertConfig.from_pretrained(
        args.model_name_or_path, num_labels=num_labels, finetuning_task=args.task_name)
    config.loss_type = args.loss_type

    if args.model_type_limit == 'Bert_Flat':
        config.cls = 'six_layers'

    if args.model_type_limit == 'Bert_CLS':
        config.cls = 'simple_head'

    if args.task_name == 'VIV_mhc_pos' or args.task_name == 'VIV_Pro':

        MHC_hidden, G2L = 768, 384
        config.g_size = MHC_hidden
        config.G2L = G2L
        config.hierarchy_classes = [3, 5, 14, 21, 20, 16]  # magic number
        config.hierarchical_depth = len(config.hierarchy_classes)

        # for hierarchy violation penalty
        fj = open(args.hierarchy_structure)
        hierarchy_structure = json.load(fj)
        hierarchy_structure = {int(k): int(v)
                               for k, v in hierarchy_structure.items()}
        config.hierarchy_structure = hierarchy_structure
        # penalty coefficient
        config.coe_l = args.coe_l

    else:
        config.g_size = args.MLP_width

    # save args in configs
    config.batch_size = args.per_device_train_batch_size * \
        args.gradient_accumulation_steps
    config.lr = args.learning_rate
    config.cls_lr = args.learning_rate_classifier
    config.cls_dropout_prob = args.CLS_dropout

    if args.model_name_or_path:
        model = model_class.from_pretrained(
            args.model_name_or_path,
            from_tf=bool(".ckpt" in args.model_name_or_path),
            config=config,
        )
    else:
        logger.info("Training new model from scratch")
        model = model_class(config)

    model.resize_token_embeddings(len(tokenizer))

    # DataLoaders creation:
    if args.pad_to_max_length:
        # If padding was already done ot max length, we use the default data collator that will just convert everything
        # to tensors.
        data_collator = default_data_collator
    else:
        # Otherwise, `DataCollatorWithPadding` will apply dynamic padding for us (by padding to the maximum length of
        # the samples passed). When using mixed precision, we add `pad_to_multiple_of=8` to pad all tensors to multiple
        # of 8s, which will enable the use of Tensor Cores on NVIDIA hardware with compute capability >= 7.5 (Volta).
        data_collator = DataCollatorWithPadding(
            tokenizer, pad_to_multiple_of=(8 if accelerator.use_fp16 else None))
        # added drop last to avoid unstable training

    # Part 2: Dataloader preparation
    if args.do_train:
        # Log a few random samples from the training set:
        for index in random.sample(range(len(train_dataset)), 2):
            logger.info(
                f"Sample {index} of the training set: {train_dataset[index]}.")

    # Freeze function: freeze several layers of model
        if args.unfrozen_layers != -1:
            for param in model.bert.embeddings.parameters():
                param.requires_grad = False

            # if unfrozen_layers == -1, we won't do the layer freeze
            if args.unfrozen_layers == 0:
                for param in model.bert.parameters():
                    param.requires_grad = False
            else:
                for layer in model.bert.encoder.layer[:-args.unfrozen_layers]:
                    for param in layer.parameters():
                        param.requires_grad = False

        train_dataloader = DataLoader(train_dataset, shuffle=True, drop_last=True,
                                      collate_fn=data_collator, batch_size=args.per_device_train_batch_size)
        eval_dataloader = DataLoader(
            eval_dataset, collate_fn=data_collator, batch_size=args.per_device_eval_batch_size)

        # Optimizer
        # Split weights in two groups, one with weight decay and the other not.
        no_decay = ["bias", "LayerNorm.weight"]

        # use different learning rate for classifier and feature-extractor
        # contributed by tianyi
        if args.learning_rate_classifier != args.learning_rate:
            optimizer_grouped_parameters = [
                {
                    "params": [p for n, p in model.bert.named_parameters() if not any(nd in n for nd in no_decay)],
                    "weight_decay": args.weight_decay,
                },
                {
                    "params": [p for n, p in model.bert.named_parameters() if any(nd in n for nd in no_decay)],
                    "weight_decay": 0.0,
                },
                {
                    "params": [p for n, p in model.classifier.named_parameters() if not any(nd in n for nd in no_decay)],
                    "weight_decay": args.weight_decay,
                    "lr": args.learning_rate_classifier,
                },
                {
                    "params": [p for n, p in model.classifier.named_parameters() if any(nd in n for nd in no_decay)],
                    "weight_decay": 0.0,
                    "lr": args.learning_rate_classifier,
                },
            ]
        else:
            optimizer_grouped_parameters = [
                {
                    "params": [p for n, p in model.named_parameters() if not any(nd in n for nd in no_decay)],
                    "weight_decay": args.weight_decay,
                },
                {
                    "params": [p for n, p in model.named_parameters() if any(nd in n for nd in no_decay)],
                    "weight_decay": 0.0,
                },
            ]

        optimizer = AdamW(optimizer_grouped_parameters, lr=args.learning_rate)

        # Prepare everything with our `accelerator`.
        model, optimizer, train_dataloader, eval_dataloader = accelerator.prepare(
            model, optimizer, train_dataloader, eval_dataloader
        )

        # Note -> the training dataloader needs to be prepared before we grab his length below (cause its length will be
        # shorter in multiprocess)

        # Scheduler and math around the number of training steps.
        num_update_steps_per_epoch = math.ceil(
            len(train_dataloader) / args.gradient_accumulation_steps)
        if args.max_train_steps is None:
            args.max_train_steps = args.num_train_epochs * num_update_steps_per_epoch
        else:
            args.num_train_epochs = math.ceil(
                args.max_train_steps / num_update_steps_per_epoch)

        # added
        step_total = args.num_train_epochs * num_update_steps_per_epoch
        warmup_steps = args.num_warmup_steps if args.warmup_percent == 0 else int(
            args.warmup_percent*step_total)

        lr_scheduler = get_scheduler(
            name=args.lr_scheduler_type,
            optimizer=optimizer,
            num_warmup_steps=args.num_warmup_steps,
            num_training_steps=args.max_train_steps,
        )

        # Train!

        d3 = tb_time.strftime("%h%d_%H.%M")
        board_writer = SummaryWriter(log_dir='runs/{}/{}_{}_LR_{}_clsLR_{}_Batch_{}_accum_{}_CV_{}_{}_seed{}'.format(args.task_name, d3, args.max_length, args.learning_rate,
                                     args.learning_rate_classifier, args.per_device_train_batch_size, args.gradient_accumulation_steps, args.cv_num, args.weighting, args.seed))

    else:
        eval_dataloader = DataLoader(
            eval_dataset, collate_fn=data_collator, batch_size=args.per_device_eval_batch_size)
        # Prepare everything with our `accelerator`.
        model, eval_dataloader = accelerator.prepare(
            model, eval_dataloader
        )

    # Part 3: Start Training
    if args.do_train:
        total_batch_size = args.per_device_train_batch_size * \
            accelerator.num_processes * args.gradient_accumulation_steps

        # dynamic logging steps
        logging_step = args.logging_steps if args.logging_per_epoch == 0 else int(
            1/args.logging_per_epoch * num_update_steps_per_epoch)

        logger.info("***** Running training *****")
        logger.info(f"  Num examples = {len(train_dataset)}")
        logger.info(f"  Num Epochs = {args.num_train_epochs}")
        logger.info(
            f"  Instantaneous batch size per device = {args.per_device_train_batch_size}")
        logger.info(
            f"  Total train batch size (w. parallel, distributed & accumulation) = {total_batch_size}")
        logger.info(
            f"  Gradient Accumulation steps = {args.gradient_accumulation_steps}")
        logger.info(f"  Total optimization steps = {args.max_train_steps}")
        logger.info(f"  logging every steps = {logging_step}")
        if args.unfrozen_layers != -1:
            logger.info(
                f"  frozen bert layers = {len(model.bert.encoder.layer) - args.unfrozen_layers}")
        # Only show the progress bar once on each machine.
        progress_bar = tqdm(range(args.max_train_steps),
                            disable=not accelerator.is_local_main_process)
        completed_steps = 0
        logging_times = 0

        # added for reproducibility
        if args.seed is not None:
            set_seed(args.seed)
        model.zero_grad()

        # logging loss file
        output_tr_loss_file = os.path.join(
            args.output_dir, "trainning_losses.tsv")

        # we use val loss perform an early stopping
        best_record = 1000
        last_record = 1000
        stop_counter = 0

        for epoch in range(args.num_train_epochs):
            # flag for eval log info
            log_time = 0
            model.train()

            # reset for each epoch
            total_train_loss, logging_loss, saving_loss = 0.0, 0.0, 0.0

            for step, batch in enumerate(train_dataloader):
                outputs = model(**batch)
                loss = outputs.loss

                # Scale the loss to the mean of the accumulated batch size
                loss = loss / args.gradient_accumulation_steps
                # backward to destory the compute map to save memory
                accelerator.backward(loss)

                # add total training loss
                total_train_loss += loss.item()

                # skip the first step
                if step == 0:
                    continue

                    # gradient_accumulation equals to huge batch
                if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:

                    # update
                    optimizer.step()
                    lr_scheduler.step()

                    # Reset gradients, for the next accumulated batches
                    optimizer.zero_grad()
                    progress_bar.update(1)
                    completed_steps += 1

                    with open(output_tr_loss_file, "a") as writer:

                        # save losses accumulated every step
                        loss_print = total_train_loss - saving_loss
                        saving_loss = total_train_loss

                        board_writer.add_scalar(
                            'Loss/train', loss_print, completed_steps)

                        tr_loss_result = "steps:" + "\t" + \
                            str(completed_steps)+"\t" + str(loss_print)
                        writer.write(tr_loss_result + "\n")

                        # add logging
                        if completed_steps > 0 and logging_step > 0 and completed_steps % logging_step == 0:

                            # logging only for logging steps
                            loss_scalar = (total_train_loss -
                                           logging_loss) / logging_step
                            logging_loss = total_train_loss
                            logger.info(
                                "completed_steps = %s, avg_logging_loss = %s ", completed_steps, loss_scalar)

                    # for eval, same as logging train
                    results = {}
                    # no evaluation
                    if args.logging_per_epoch == 0:
                        if completed_steps >= args.max_train_steps:
                            break
                        else:
                            continue

                if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:

                    # logging
                    if completed_steps > 0 and logging_step > 0 and completed_steps % logging_step == 0:

                        # evaluation
                        log_time += 1
                        log_info = 'epoch{}-{}'.format(epoch+1, log_time)
                        result, eval_loss = evaluate(
                            args, model, eval_dataloader, completed_steps, tb_time, log_info)

                        logging_times += 1
                        results.update(result)

                        # early stop
                        if args.early_stop != -1:

                            if eval_loss <= best_record:
                                # update best, reset counter
                                best_record = eval_loss
                                stop_counter = 0

                                # Save model checkpoint
                                checkpoint_prefix = "best_checkpoint"
                                output_check_dir = os.path.join(
                                    args.output_dir, "best_checkpoint-{}".format(log_info))
                                if not os.path.exists(output_check_dir):
                                    os.makedirs(output_check_dir)

                                accelerator.wait_for_everyone()
                                unwrapped_model = accelerator.unwrap_model(
                                    model)
                                unwrapped_model.save_pretrained(
                                    output_check_dir, save_function=accelerator.save)

                                logger.info(
                                    "Saving model best checkpoint to %s", output_check_dir)
                                _rotate_checkpoints(args, checkpoint_prefix)

                                torch.save(args, os.path.join(
                                    output_check_dir, "training_args.bin"))
                                # torch.save(optimizer.state_dict(), os.path.join(output_check_dir, "optimizer.pt"))
                                torch.save(lr_scheduler.state_dict(), os.path.join(
                                    output_check_dir, "scheduler.pt"))
                                logger.info(
                                    "Saving best optimizer and scheduler states to %s", output_check_dir)

                            else:
                                stop_counter += 1

                            # record the last val loss

                            last_record = eval_loss

                            if stop_counter == args.early_stop:
                                logger.info("Early stop!")
                                logger.info("Best val loss: {}\n Last val loss: {}".format(
                                    best_record, last_record))
                                break

                        model.train()
                # eval end

                if completed_steps >= args.max_train_steps:
                    break
            # epoch end

            if stop_counter == args.early_stop:
                logger.info("Early stop")
                logger.info("Best val loss: {}\n Last val loss: {}".format(
                    best_record, last_record))
                break

            # save checkpoint， use completed steps
            else:
                checkpoint_prefix = "checkpoint"
                # Save model checkpoint
                output_check_dir = os.path.join(
                    args.output_dir, "checkpoint-epoch-{}-end".format(epoch+1))
                if not os.path.exists(output_check_dir):
                    os.makedirs(output_check_dir)

                accelerator.wait_for_everyone()
                unwrapped_model = accelerator.unwrap_model(model)
                unwrapped_model.save_pretrained(
                    output_check_dir, save_function=accelerator.save)
                # tokenizer.save_pretrained(output_check_dir)

                logger.info("Saving model checkpoint to %s", output_check_dir)

                _rotate_checkpoints(args, checkpoint_prefix)

                torch.save(args, os.path.join(
                    output_check_dir, "training_args.bin"))
                # torch.save(optimizer.state_dict(), os.path.join(output_check_dir, "optimizer.pt"))
                torch.save(lr_scheduler.state_dict(), os.path.join(
                    output_check_dir, "scheduler.pt"))
                logger.info(
                    "Saving optimizer and scheduler states to %s", output_check_dir)

                logger.info(
                    f"  Average training losses for this epoch = {total_train_loss * args.gradient_accumulation_steps / len(train_dataloader)}")

                # stop training early
                if epoch+1 == args.early_stop_epoch_manual:
                    logger.info("Early stop at epoch-{}-end".format(epoch+1))
                    break

            if args.output_dir is not None:
                accelerator.wait_for_everyone()
                unwrapped_model = accelerator.unwrap_model(model)
                unwrapped_model.save_pretrained(
                    args.output_dir, save_function=accelerator.save)

    # Part 4: Predicting
    prediction = {}
    if args.do_predict:
        logger.info("Predict using the saved model")
        prediction = predict(args, model, eval_dataloader)

    # Predict no label!
    if args.do_predict_no_label:
        logger.info("Inference  using the saved model")
        predict_no_label(args, model, eval_dataloader)

    # Visualize!
    if args.do_visualize:
        logger.info(
            "Calculate attention score and visualize using the saved model")
        scores = None
        attention_scores, probs = visualize(args, model, eval_dataloader)
        if scores is not None:
            scores = attention_scores
        else:
            scores = deepcopy(attention_scores)
        np.save(os.path.join(args.output_dir, "atten.npy"), scores)
    # Embedding!
    if args.do_embedding:
        logger.info("Output Embedding for sequence using the saved model")
        give_embedding(args, model, eval_dataloader)

    # tensorboard
    if args.do_train:
        board_writer.flush()
        board_writer.close()


if __name__ == "__main__":
    main()
