#!/bin/bash
export NUM=$1
export FEATURE=GO

###### search hyper-parameters
export LR=$2
export BS=$3
export Dropout=$4

###### data and result directory
export DATA=/data1/jyguo/VBird/BERT/Scripts/SLim/clustered_go/
export TRAIN_FILE=${DATA}/${NUM}/train.csv
export TEST_FILE=${DATA}/${NUM}/test.csv
export OUTPUT_FOLDER=/data1/jyguo/VBird/BERT/trained_model/GO/APR/
export OUTPUT_PATH=${OUTPUT_FOLDER}/LR_${LR}_BS_${BS}_DR_${Dropout}/${NUM}
export LOG_PATH=/data1/jyguo/VBird/BERT/logger/VIV_slim_LR_$(date "+%m.%d").log
export SEED_TRAIN=42

mkdir $OUTPUT_FOLDER
mkdir $OUTPUT_PATH


CUDA_VISIBLE_DEVICES=$5 python run_slim.py \
    --cv_num $NUM \
    --feature_type $FEATURE \
    --train_file $TRAIN_FILE \
    --validation_file $TEST_FILE \
    --output_dir $OUTPUT_PATH \
    --save_path $OUTPUT_FOLDER \
    --do_train \
    --dropout_prob $Dropout \
    --seed $SEED_TRAIN \
    --gradient_accumulation_steps 1 \
    --per_device_train_batch_size $BS \
    --per_device_eval_batch_size 64 \
    --num_train_epochs 100 \
    --warmup_percent 0.1 \
    --learning_rate $LR \
    --adam_epsilon 1e-6 \
    --weight_decay 0.01 \
    --beta1 0.9 \
    --beta2 0.98 \
    --preprocessing_num_workers 128 \
    --logger_path $LOG_PATH \
    --logging_per_epoch 1 \
    --save_total_limit 1 \
    --save_every_log 2000 \
    --fp16 True \
    --early_stop 20 \
    --loss_type penalty \