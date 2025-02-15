#!/bin/bash

# Prompt the user for the project name
echo "Please enter the project name:"
read project_name

# Set the output directory
export o=Out/

# Run the Python script with the specified project name and output directory
python Predict/extract_prediction.py \
    --output_dir ${o}/${project_name} \
    --input_dir  # You can specify the input directory here if needed