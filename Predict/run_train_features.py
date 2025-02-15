""" HCMN-F training script
"""
from mhc_model import MHC_Net
from transformers import (
    AdamW,
    DataCollatorWithPadding,
    default_data_collator,
    SchedulerType,
    get_scheduler,
    set_seed,
)
from transformers import bioinfo_compute_metrics as compute_metrics
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
from datetime import datetime
import argparse
import json

import configparser

#  tensorboard
from torch.utils.tensorboard import SummaryWriter
from accelerate import Accelerator, DeepSpeedPlugin

import transformers
print(transformers.__version__)

#  added metrics


# import Network
CUDA_VISIBLE_DEVICES = 0, 1, 2, 3

MODEL_CLASSES = {
    "Genomic_Traits": (MHC_Net, 4229, 768),
    "Protein_Traits": (MHC_Net, 11201, 1024),
    "PPI": (MHC_Net, 484, 768),
}


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
        description="Training a HCMN model on a multi-label classification task")

# 1. data & files
    parser.add_argument(
        "--feature_type",
        type=str,
        default=None,
        help="The name of feature used for the virus prediction to train on.",
    )
    parser.add_argument(
        "--train_file", type=str, default=None, help="A csv or a json file containing the training data."
    )
    parser.add_argument(
        "--validation_file", type=str, default=None, help="A csv or a json file containing the validation data."
    )
    parser.add_argument(
        "--model_name_or_path",
        type=str,
        help="Path to pretrained model or model identifier from huggingface.co/models.",
    )
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Where to store the final model.")
    parser.add_argument("--cv_num", default=0, type=str,
                        help="cross validation fold")
    parser.add_argument(
        "--logger_path",
        default='Logger/default_log_file.log',
        type=str,
        help="Save the logger output to files",
    )
    parser.add_argument(
        "--record_path",
        type=str,
        help="Save the model records",
    )
    parser.add_argument("--logging_steps", type=int,
                        default=0, help="Log every X updates steps.")
    parser.add_argument("--logging_per_epoch", type=int,
                        default=0, help="Log per epoch.")
    parser.add_argument("--save_every_log", type=int,
                        default=1, help="Save checkpoint every X logs.")
    parser.add_argument("--save_every_epoch", action="store_true",
                        help="Whether to save every epoch.")
    parser.add_argument(
        "--save_total_limit",
        type=int,
        default=None,
        help="Limit the total amount of checkpoints, delete the older checkpoints in the output_dir, does not delete by default",)

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
    parser.add_argument("--weight_decay", type=float,
                        default=0.0, help="Weight decay to use.")
    parser.add_argument(
        "--learning_rate",
        type=float,
        default=1e-3,
        help="learning rate.",
    )
    parser.add_argument(
        "--max_train_steps",
        type=int,
        default=None,
        help="Total number of training steps to perform. If provided, overrides num_train_epochs.",
    )
    parser.add_argument("--num_train_epochs", type=int, default=3,
                        help="Total number of training epochs to perform.")
    parser.add_argument(
        "--num_warmup_steps", type=int, default=0, help="Number of steps for the warmup in the lr scheduler."
    )
    parser.add_argument("--warmup_percent", default=0, type=float,
                        help="Linear warmup over warmup_percent*total_steps.")
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
    parser.add_argument("--beta1", default=0.9, type=float,
                        help="Beta1 for Adam optimizer.")
    parser.add_argument("--beta2", default=0.999, type=float,
                        help="Beta2 for Adam optimizer.")
    # added for class balance
    parser.add_argument("--weighting", default=None,
                        type=str, help="Class Balancing methods.")
    # added for training and logging
    parser.add_argument("--dropout_prob", type=float,
                        help="Dropout rate of the HMCN net.")
    parser.add_argument(
        "--multi_threshold",
        type=float,
        default=0.5,
        help="Give threshold for sigmoid"
    )
    parser.add_argument(
        "--early_stop", default=-1, type=int, help="set this to a positive integet if you want to perfrom early stop. The model will stop \
                                                    if the auc keep decreasing early_stop times",
    )

# 3. model and config

    # added for penalty loss
    parser.add_argument(
        "--hierarchy_structure",
        type=str,
        help="Give hierarchy structure path"
    )
    parser.add_argument("--loss_type", default="CELoss",
                        type=str, help="Type of loss used")
    parser.add_argument("--G2L", type=int, required=True,
                        help="Local layer size, Global to Local layer projection.")
    parser.add_argument("--g_size", type=int, required=True,
                        help="Global layer size.")
#     parser.add_argument("--hidden_size", type=int, required=True, help="Input layer size of the HMCN net.")
    parser.add_argument("--label_num", type=int,
                        help="Predict mode label numbers.")

# 4. working modes
    parser.add_argument("--do_train", action="store_true",
                        help="Whether to run training.")
    parser.add_argument("--do_predict", action="store_true",
                        help="Whether to do prediction on the given dataset.")
    parser.add_argument("--do_predict_no_label", action="store_true",
                        help="Whether to do prediction on the given dataset with no label information.")

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
    if args.feature_type is None and args.train_file is None and args.validation_file is None:
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

#  save model


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

    logger.info("***** Running evaluation {} *****")
    logger.info("  Num examples = %d", len(eval_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    eval_loss, avg_steps = 0, 0
    logits = None

    model.eval()
    # torch.no_grad()
    # need to eval or all the examples

    # add eval_prob_false_result
    output_eval_nan_file = os.path.join(args.output_dir, "eval_nan_loss.txt")

    for eval_step, eval_batch in enumerate(tqdm(eval_dataloader)):
        #  save GPU memory
        with torch.no_grad():
            outputs = model(**eval_batch)
        if logits is None:
            logits = outputs['logits'].detach().cpu().numpy()
            reference = eval_batch["labels"].detach().cpu().numpy()
        else:
            logits = np.append(
                logits, outputs['logits'].detach().cpu().numpy(), axis=0)
            reference = np.append(
                reference, eval_batch["labels"].detach().cpu().numpy(), axis=0)

        # sanity check, avoid nan for logits
        loss_add = outputs['loss'].item()
        if loss_add == loss_add:
            eval_loss += loss_add
            avg_steps += 1

    # avg
    eval_loss = eval_loss / avg_steps

    # add eval_loss_result
    output_eval_loss_file = os.path.join(args.output_dir, "eval_losses.txt")

    d3 = tb_time.strftime("%h%d_%H.%M")
    board_writer = SummaryWriter(log_dir='runs/{}/{}_LR{}_Batch{}_({},{},{},{})_CV{}_{}_seed{}'.format(args.feature_type, d3, args.learning_rate,
                                 args.per_device_train_batch_size, args.hidden_size, args.g_size, args.G2L, args.dropout_prob, args.cv_num, args.weighting, args.seed))
    board_writer.add_scalar('Loss/eval', eval_loss, completed_steps)

    with open(output_eval_loss_file, "a") as writer:
        writer.write(str(eval_loss) + "\n")

    probs = torch.sigmoid(torch.tensor(logits, dtype=torch.float32)).numpy()

    # post_processing of hierarchical violation
    f = open(args.hierarchy_structure, )
    taxa = json.load(f)
    taxa = {int(k): int(v) for k, v in taxa.items()}
    probs = violation_process(taxa, probs)

    upper, lower = 1, 0
    THRESHOLD = args.multi_threshold
    preds = np.where(probs > THRESHOLD, upper, lower)

    result = compute_metrics("mhc", preds, reference, probs)

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

                writer.write(str(result[key]) + ",")
        writer.write("\n")

    logger.info("Confusion Matrix : \n {}".format(result['confusion_matrix']))
    np.save(output_eval_file_cm, result['confusion_matrix'])
    np.save(output_eval_file_probs, probs)
    #  add new values to dict

    # tensorboard
    board_writer.add_pr_curve('pr_curve', reference, probs, completed_steps)

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
    board_writer.close()

    return result, eval_loss


def predict(args, model, pred_dataloader):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running predict {} *****")
    logger.info("  Num examples = %d", len(pred_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)

    logits = None

    model.eval()
    # sanity check, avoid nan loss
    output_pred_nan_file = os.path.join(args.output_dir, "pred_nan_loss.txt")

    for pred_step, pred_batch in enumerate(tqdm(pred_dataloader)):
        # save GPU memory in predict
        with torch.no_grad():
            outputs = model(**pred_batch)
        if logits is None:
            logits = outputs['logits'].detach().cpu().numpy()
            reference = pred_batch["labels"].detach().cpu().numpy()
        else:
            logits = np.append(
                logits, outputs['logits'].detach().cpu().numpy(), axis=0)
            reference = np.append(
                reference, pred_batch["labels"].detach().cpu().numpy(), axis=0)

        # sanity check, avoid nan for logits
        loss_add = outputs['loss'].item()

    probs = torch.sigmoid(torch.tensor(logits, dtype=torch.float32)).numpy()
    # post_processing of hierarchical violation
    f = open(args.hierarchy_structure, )
    taxa = json.load(f)
    taxa = {int(k): int(v) for k, v in taxa.items()}
    probs = violation_process(taxa, probs)

    upper, lower = 1, 0
    THRESHOLD = args.multi_threshold
    preds = np.where(probs > THRESHOLD, upper, lower)
    result = compute_metrics("mhc", preds, reference, probs)
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
                writer.write(str(result[key]) + ",")
        writer.write("\n")

    logger.info("Confusion Matrix : \n {}".format(result['confusion_matrix']))
    np.save(output_pred_file_cm, result['confusion_matrix'])
    np.save(output_pred_file_probs, probs)
    np.save(output_reference, reference)

    return result


def predict_no_label(args, model, pred_dataloader):
    softmax = torch.nn.Softmax(dim=1)

    logger.info("***** Running predict {} *****")
    logger.info("  Num examples = %d", len(pred_dataloader))
    logger.info("  Batch size = %d", args.per_device_eval_batch_size)
    logits = None

    model.eval()

    for pred_step, pred_batch in enumerate(tqdm(pred_dataloader)):
        # save GPU memory
        with torch.no_grad():
            outputs = model(**pred_batch)
        if logits is None:
            logits = outputs['logits'].detach().cpu().numpy()
        else:
            logits = np.append(
                logits, outputs['logits'].detach().cpu().numpy(), axis=0)

    probs = torch.sigmoid(torch.tensor(logits, dtype=torch.float32)).numpy()
    # post_processing of hierarchical violation
    f = open(args.hierarchy_structure, )
    taxa = json.load(f)
    taxa = {int(k): int(v) for k, v in taxa.items()}
    probs = violation_process(taxa, probs)

    upper, lower = 1, 0
    THRESHOLD = args.multi_threshold
    preds = np.where(probs > THRESHOLD, upper, lower)

    # define pred output file name
    output_pred_file_probs = os.path.join(args.output_dir, "pred_probs.npy")
    np.save(output_pred_file_probs, probs)

    return None


def main():
    args = parse_args()
    tb_time = datetime.now()
    config = configparser.ConfigParser()
    #  write logging info to path
    global logger
    logger = get_logger(args)

    #  add deepspeed
    deepspeed_plugin = DeepSpeedPlugin(
        zero_stage=2, gradient_accumulation_steps=args.gradient_accumulation_steps)

    accelerator = Accelerator(
        fp16=args.fp16, deepspeed_plugin=deepspeed_plugin)

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

    data_files = {}
    if args.train_file is not None:
        data_files["train"] = args.train_file
    if args.validation_file is not None:
        data_files["validation"] = args.validation_file
    extension = (
        args.train_file if args.train_file is not None else args.validation_file).split(".")[-1]
    raw_datasets = load_dataset(extension, data_files=data_files)

    def Translate_feature(f):
        feature = f.split('_')
        feature = [float(i) for i in feature]
        return feature

    def Translate_multi_label(l):
        label = l.split('_')
        label = [int(i) for i in label]
        return label

    def Translate_multi_weight(w):
        weight = w.split('_')
        weight = [float(i) for i in weight]
        return weight

    def preprocess_multilabel_function(examples):

        result = {}

        if "genomic_features" in examples:
            result["features"] = list(
                map(lambda l: Translate_feature(l), examples['genomic_features']))
        if "proteomic_features" in examples:
            result["features"] = list(
                map(lambda l: Translate_feature(l), examples['proteomic_features']))
        if "structure_features" in examples:
            result["features"] = list(
                map(lambda l: Translate_feature(l), examples['structure_features']))

        if "embedding" in examples:
            result["features"] = list(
                map(lambda l: Translate_feature(l), examples['embedding']))

        if "label" in examples:
            # translate labels to list
            result["labels"] = list(
                map(lambda l: Translate_multi_label(l), examples['label']))

        if "weight" in examples:
            result["weights"] = list(
                map(lambda l: Translate_multi_weight(l), examples['weight']))

        return result

    # Part 1: Three Load data modes

    if args.do_train:

        new_labels = list(map(lambda l: Translate_multi_label(
            l), raw_datasets['train']['label']))
        num_labels = len(new_labels[0])

        processed_datasets = raw_datasets.map(
            preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["train"].column_names, )

        train_dataset = processed_datasets["train"]
        eval_dataset = processed_datasets["validation"]

    elif args.do_predict:
        new_labels = list(map(lambda l: Translate_multi_label(
            l), raw_datasets['validation']['label']))
        num_labels = len(new_labels[0])

        processed_datasets = raw_datasets.map(
            preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
        )

        eval_dataset = processed_datasets["validation"]

    elif args.do_predict_no_label:
        num_labels = args.label_num  # magic number
        processed_datasets = raw_datasets.map(
            preprocess_multilabel_function, batched=True, remove_columns=raw_datasets["validation"].column_names,
        )

        eval_dataset = processed_datasets["validation"]

    # args

    # MHC args

    config.g_size = args.g_size
    config.G2L = args.G2L

    if args.feature_type == 'Genomic_Traits':
        config.hierarchy_classes = [3, 5, 14, 21, 20, 16]
    elif args.feature_type == 'Protein_Traits':
        config.hierarchy_classes = [3, 5, 14, 21, 20, 16]

    elif args.feature_type == 'PPI':
        config.hierarchy_classes = [1, 1, 1, 1, 1, 1]

    config.hierarchical_depth = len(config.hierarchy_classes)

    def jsonKeys2int(x):
        if isinstance(x, dict):
            return {int(k): v for k, v in x.items()}
        return x
    fj = open(args.hierarchy_structure)
    hierarchy_structure = json.load(fj)
    hierarchy_structure = jsonKeys2int(hierarchy_structure)

    config.hierarchy_structure = hierarchy_structure
    # penalty coefficient
    config.coe_l = 0.1
    config.num_labels = num_labels
    config.cls_dropout_prob = args.dropout_prob

    config.loss_type = args.loss_type

    model_class, args.feature_dim, args.hidden_size = MODEL_CLASSES[args.feature_type]
    config.hidden_size = args.hidden_size
    config.feature_dim = args.feature_dim

    model = model_class(config)

    data_collator = default_data_collator

    # Part 2: Dataloader preparation
    if args.do_train:
        #         # Log a few random samples from the training set:
        #         for index in random.sample(range(len(train_dataset)), 2):
        #             logger.info(f"Sample {index} of the training set: {train_dataset[index]}.")

        train_dataloader = DataLoader(
            train_dataset, shuffle=True, drop_last=True, collate_fn=data_collator, batch_size=args.per_device_train_batch_size
        )
        eval_dataloader = DataLoader(
            eval_dataset, collate_fn=data_collator, batch_size=args.per_device_eval_batch_size)

        # Optimizer
        # Split weights in two groups, one with weight decay and the other not.
        no_decay = ["bias", "LayerNorm.weight"]

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

        num_update_steps_per_epoch = math.ceil(
            len(train_dataloader) / args.gradient_accumulation_steps)
        if args.max_train_steps is None:
            args.max_train_steps = args.num_train_epochs * num_update_steps_per_epoch
        else:
            args.num_train_epochs = math.ceil(
                args.max_train_steps / num_update_steps_per_epoch)

        #  added
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
        board_writer = SummaryWriter(log_dir='runs/{}/{}_LR{}_Batch{}_({},{},{},{})_CV{}_{}_seed{}'.format(args.feature_type, d3, args.learning_rate,
                                     args.per_device_train_batch_size, args.hidden_size, args.g_size, args.G2L, args.dropout_prob, args.cv_num, args.weighting, args.seed))

    # predict mode
    else:

        model.load_state_dict(torch.load(
            args.model_name_or_path+'/best_model.pt'))
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

        #  dynamic logging steps
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
        # Only show the progress bar once on each machine.
        progress_bar = tqdm(range(args.max_train_steps),
                            disable=not accelerator.is_local_main_process)
        completed_steps = 0
        logging_times = 0

        #  added for reproducibility
        if args.seed is not None:
            set_seed(args.seed)
        model.zero_grad()

        #  logging loss file
        output_tr_loss_file = os.path.join(
            args.output_dir, "trainning_losses.tsv")

        best_record = 1000
        last_record = 1000
        stop_counter = 0
        for epoch in range(args.num_train_epochs):
            model.train()
            log_time = 0
            #  reset for each epoch
            total_train_loss, logging_loss, saving_loss = 0.0, 0.0, 0.0

            for step, batch in enumerate(train_dataloader):
                outputs = model(**batch)
                loss = outputs['loss']

                #  Scale the loss to the mean of the accumulated batch size
                loss = loss / args.gradient_accumulation_steps
                #  backward to destory the compute map to save memory
                accelerator.backward(loss)

                #  add total training loss
                total_train_loss += loss.item()

                #  skip the first step
                if step == 0:
                    continue

                #  gradient_accumulation equals to huge batch
                if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:

                    #  update
                    optimizer.step()
                    lr_scheduler.step()

                    #   Reset gradients, for the next accumulated batches
                    optimizer.zero_grad()
                    progress_bar.update(1)
                    completed_steps += 1

                if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:
                    with open(output_tr_loss_file, "a") as writer:

                        #  save losses accumulated every step
                        loss_print = total_train_loss - saving_loss
                        saving_loss = total_train_loss

                        board_writer.add_scalar(
                            'Loss/train', loss_print, completed_steps)

                        tr_loss_result = "steps:" + "\t" + \
                            str(completed_steps)+"\t" + str(loss_print)
                        writer.write(tr_loss_result + "\n")

                        #  add logging
                        if completed_steps > 0 and logging_step > 0 and completed_steps % logging_step == 0:

                            #  logging only for logging steps
                            loss_scalar = (total_train_loss -
                                           logging_loss) / logging_step
                            logging_loss = total_train_loss
                            logger.info("completed_steps = %s, avg_logging_loss = %s ", completed_steps, loss_scalar)

                #  for eval, same as logging train
                results = {}
                if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:

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
                                opt_auc = result['auc_micro']
                                opt_aupr = result['aupr_micro']
                                best_valid_loss = eval_loss
                                log_info = 'best_epoch{}-{}'.format(
                                    epoch+1, log_time)
                                # result, eval_loss = evaluate(args, model, eval_dataloader, completed_steps, tb_time, log_info)

                                accelerator.wait_for_everyone()
                                unwrapped_model = accelerator.unwrap_model(
                                    model)
                                torch.save(unwrapped_model.state_dict(), os.path.join(
                                    args.output_dir, 'best_model.pt'))
                                logger.info(
                                    "Saving model checkpoint to %s", args.output_dir)
                                torch.save(args, os.path.join(
                                    args.output_dir, "training_args.bin"))
                                torch.save(optimizer.state_dict(), os.path.join(
                                    args.output_dir, "optimizer.pt"))
                                torch.save(lr_scheduler.state_dict(), os.path.join(
                                    args.output_dir, "scheduler.pt"))

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
                #  eval end

                #  save checkpoint， use completed steps
                if args.logging_per_epoch != 0:
                    if completed_steps > 0 and args.save_every_log > 0 and logging_times % args.save_every_log == 0 and completed_steps % logging_step == 0:
                        #  save only once
                        if step % args.gradient_accumulation_steps == 0 or step == len(train_dataloader) - 1:

                            checkpoint_prefix = "checkpoint"
                            # Save model checkpoint
                            output_check_dir = os.path.join(
                                args.output_dir, "checkpoint-{}".format(log_info))
                            if not os.path.exists(output_check_dir):
                                os.makedirs(output_check_dir)

                if completed_steps >= args.max_train_steps:
                    break

            # epoch end
            if stop_counter == args.early_stop:
                logger.info("Early stop")
                logger.info("Best val loss: {}\n Last val loss: {}".format(
                    best_record, last_record))
                break

            # save every epoch

            if args.save_every_epoch:
                accelerator.wait_for_everyone()
                unwrapped_model = accelerator.unwrap_model(model)
                torch.save(unwrapped_model.state_dict(), os.path.join(
                    args.output_dir, f'model_epoch_{epoch+1}.pt'))
                logger.info("Saving model checkpoint to %s", args.output_dir)
                torch.save(args, os.path.join(
                    args.output_dir, "training_args.bin"))
                torch.save(optimizer.state_dict(), os.path.join(
                    args.output_dir, "optimizer.pt"))
                torch.save(lr_scheduler.state_dict(), os.path.join(
                    args.output_dir, "scheduler.pt"))

            logger.info(
                f"  Average training losses for this epoch = {total_train_loss * args.gradient_accumulation_steps / len(train_dataloader)}")

        # save final results of early stop
        outputFilename = args.record_path + '/results.txt'

        # save best performance
        if os.path.exists(outputFilename) == False:
            header = "model\tAUC\tAUPR\tLoss\n"
            with open(outputFilename, 'a') as file_writer:
                file_writer.write(header)
                file_writer.close()

        f = open(outputFilename, 'a')
        f.write(
            f'Model_{args.cv_num}_{args.per_device_train_batch_size}_{args.learning_rate}_{args.dropout_prob}:\t')
        f.write(f'{opt_auc:.4f}\t')
        f.write(f'{opt_aupr:.4f}\t')
        f.write(f'{best_valid_loss:.4f}\n')
        f.close()

    # Part 4: Predicting
    # Predict!
    prediction = {}
    if args.do_predict:
        logger.info("Predict using the saved model")
        prediction = predict(args, model, eval_dataloader)

    # Predict no label!
    if args.do_predict_no_label:
        logger.info("Inference using the saved model")
        predict_no_label(args, model, eval_dataloader)

    # tensorboard
    if args.do_train:
        board_writer.flush()
        board_writer.close()


if __name__ == "__main__":
    main()
