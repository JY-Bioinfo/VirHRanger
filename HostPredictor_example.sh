#!/bin/bash

# Your script logic starts here
printf "VirHRanger Version-0.1"

################### 1. Experimental Variables ###################
export TRANSFORMERS_OFFLINE=1
export MAX_PROC=$(nproc --all) # get cores
export g=3
export b=8 # batch size
export c=$((MAX_PROC / 2)) # use half of max core
export m=Model/
export Tool_Path=${PWD}


# default values
export i=Example_input/Bats
export o=Out/
export e=Genomic_Protein
export f=Genomic_FM # default value
export w=whole
export p=individual




################### 2. Options ###################

usage() { printf -k "Usage: $0
[-h "help"]
[-i "input_path" <string>]
[-o "output_path" <string>]  
[-m "model_path" <string>]
[-c "cores_used" <int>]
[-e "ensemble modes: choose from {Genomic, Genomic_Protein, Genomic_Protein_PPI}" <string>]
[-f "feature extraction modes: choose from {Genomic_FM, Protein_FM, Genomic_Traits, Protein_Traits, Homology, PPI}" <string>]
[-p "processing input by individual file or batch: choose from {individual, batch}" <string>]
[-w "whether to run whole process or skip the feature extraction and continue from last intermediate features: choose from {whole, only_predict}" <string>]
[-g "which gpu = 0,1,2,3 or cpu = -1" <int>]" 1>&2; exit 0; }


while getopts ":i:o:m:c:e:f:p:w:g:" para; do
    case "${para}" in
        i)                           
            i=${OPTARG} 
            ;;
        o) 
            o=${OPTARG}
            ;;
        m) 
            m=${OPTARG}
            ;;
        c) 
            c=${OPTARG}
            ;;
        e) 
            e=${OPTARG}
            ;;
        f) 
            f=${OPTARG}
            ;;
        p) 
            p=${OPTARG}
            ;;
        w) 
            w=${OPTARG}
            ;;
        g) 
            g=${OPTARG}
            ;;
        *) 
            usage
            ;;
            
    esac
done
shift $((OPTIND-1))

### check for the valid input options
if [[ $e != @(None|Genomic|Genomic_Protein|Genomic_Protein_PPI) ]]; then
    printf " -e \"${e}\" is not an valid option! Please choose from {Genomic, Genomic_Protein, Genomic_Protein_PPI, None}.\n"
    exit 0; 
fi

if [[ $f != @(Genomic_FM|Protein_FM|Genomic_Traits|Protein_Traits|Homology|PPI) ]]; then
    printf " -f \"${f}\" is not an valid option! Please choose from {Genomic_FM, Protein_FM, Genomic_Traits, Protein_Traits, Homology, PPI}.\n"
    exit 0; 
fi

if [[ $p != @(individual|batch) ]]; then
    printf " -p \"${p}\" is not an valid option! Please choose from {individual, batch}.\n"
    exit 0; 
fi

if [[ $w != @(whole|only_predict) ]]; then
    printf " -w \"${w}\" is not an valid option! Please choose from {whole, only_predict}.\n"
    exit 0; 
fi


export i=${i}/
### create output folder
mkdir -p ${o}
export o=${o}/Host_prediction_"$(basename -- ${i%%.*})"/
mkdir -p ${o}
mkdir -p ${o}/sequence/
export LOG_PATH=${o}/logging_info.txt


### check if input valid

if [[ ! -f ${i} ]] && [[ ! -d ${i} ]]
then
    printf "${i} does not exist on your filesystem."
    exit 0; 
fi


### input types of files by mode
export genome=FALSE
export genbank=FALSE
export protein=FALSE


if [[ $e == "Genomic" ]]; then
    export genome=TRUE
    
    
elif [[ $e == "Genomic_Protein" ]] || [[ $e == "Genomic_Protein_PPI" ]] ; then
    export genome=TRUE
    export genbank=TRUE
    export protein=TRUE
fi

### input types of files by given
if [[ $f == "Genomic_FM" ]] || [[ $f == "Homology" ]]; then
    export genome=TRUE  

elif [[ $f == "Protein_FM" ]] || [[ $f == "Protein_Traits" ]] || [[ $f == "PPI" ]]; then
    export protein=TRUE 

elif [[ $f == "Genomic_Traits" ]]; then
    export genome=TRUE
    export genbank=TRUE
fi

### check features needed
printf "Need genome: $genome\n"
printf "Need genbank: $genbank\n"
printf "Need protein: $protein\n\n"

if [[ $genome == 'FALSE' ]] && [[ $genbank == 'FALSE' ]] && [[ $protein == 'FALSE' ]];
then
    printf "No feature selected! need correct -f or -e parameters"
    exit 0; 
fi

################### 3. Preprocess sequences ###################
### only reprocess when in whole mode
if [[ $w == 'only_predict' ]]; then
    printf "Skip preprocessing steps! \n"
fi

if [[ $w == 'whole' ]]; then
    mkdir -p ${o}/sequence/intermediate

    if [[ $genome == "TRUE" ]]; then
        mkdir -p ${o}/sequence/genome
    fi
    if [[ $genbank == "TRUE" ]]; then
        mkdir -p ${o}/sequence/genbank 
    fi
    if [[ $protein == "TRUE" ]]; then
        mkdir -p ${o}/sequence/protein
        mkdir -p ${o}/sequence/all_protein
    fi

    ## remove old

    printf "Remove old sequences cached\n\n"
    rm -f ${o}/sequence/genome/*.f*a
    rm -f ${o}/sequence/genbank/*.gb*
    rm -f ${o}/sequence/protein/*.f*a
    rm -f ${o}/sequence/all_protein/*.f*a

    ### processing input modes ######

    if [[ $p == "individual" ]]; then
        printf "Individual Mode.\n\n"
    elif [[ $p == "batch" ]]; then
        printf "Batch Mode.\n\n"
    fi

    ### 3.1.processing genome files
    if [[ $genome == "TRUE" ]]; then
        printf "Processing genome sequences..."
        # copy files
        for dir in ${i}*/
        do

            dir=${dir%*/}      # remove the trailing "/"
            cp ${dir}/genome/*.f*a -t ${o}/sequence/intermediate/  # copy to intermediate dir 
            # merge all fasta into one and rename with assembly name
            cat ${o}/sequence/intermediate/*.f*a > ${o}/sequence/genome/$(basename ${dir}).fna

            rm -f ${o}/sequence/intermediate/*    # remove cached files in intermediate
        done



        ###### preprocessing files #######
        cd ${o}/sequence/genome/
        for file in *.f*a; do mv -n "$file" "${file// /_}" 2>/dev/null; done    ## remove space in filename
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a                 ## remove newline char
        sed -i '/^$/d' *.f*a ## remove empty lines
        gawk -i inplace 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' *.f*a    ## uppercase
        sed -i '/^>/!s/[RHNYWBMSVKD]/N/g' *.f*a   ## avoid out-of-vocabulary
        cd ${Tool_Path}  ## change working dir
    fi

    ### 3.2.processing protein files
    if [[ $protein == "TRUE" ]]; then
        printf "Processing protein sequences..."
        # copy files
        for dir in ${i}*/
        do

            dir=${dir%*/}      # remove the trailing "/"
            # Check if any files match the pattern
            if ls ${dir}/protein/*.f*a 1> /dev/null 2>&1; then
                cp ${dir}/protein/*.f*a -t ${o}/sequence/intermediate/  ## copy to intermediate dir 
                # merge all fasta into one and rename with assembly name
                cat ${o}/sequence/intermediate/*.f*a > ${o}/sequence/protein/$(basename ${dir}).faa
                rm -f ${o}/sequence/intermediate/*    ## remove cached files in intermediate
            fi
        done

        ###### preprocessing files #######

        cd ${o}/sequence/protein/
        for file in *.f*a; do mv -n "$file" "${file// /_}" 2>/dev/null; done    ## remove space in filename
        #sed -i '/>/s/\ /_/g' *.f*a    ## process header split with '_'
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a  ## remove newline char
        sed -i '/^$/d' *.f*a ## remove empty lines
        gawk -i inplace 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' *.f*a ## uppercase
        sed -i '/^>/!s/J/X/g' *.f*a  ## avoid out-of-vocabulary
        cd ${Tool_Path}
    fi

    ### 3.3.processing protein files separately
    if [[ $protein == "TRUE" ]]; then
        printf "Processing protein sequences separately..."
        # copy files
        for dir in ${i}*/
        do
            dir=${dir%*/}      # remove the trailing "/"
            if ls ${dir}/protein/*.f*a 1> /dev/null 2>&1; then
                cp ${dir}/protein/*.f*a -t ${o}/sequence/all_protein/  ## copy separately
                cp ${dir}/protein/*.f*a -t ${o}/sequence/all_protein_with_header/  ## copy separately
            fi
        done

        ###### preprocessing files #######
        cd ${o}/sequence/all_protein/
        for file in *.f*a; do mv -n "$file" "${file// /_}" 2>/dev/null; done    ## remove space in filename
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a  ## remove newline char
        sed -i '/^>/d' *.f*a ## remove headers
        sed -i '/^$/d' *.f*a ## remove empty lines
        gawk -i inplace 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' *.f*a ## uppercase
        sed -i '/^>/!s/J/X/g' *.f*a  ## avoid out-of-vocabulary
        cd ${Tool_Path}

        ###### preprocessing files #######
        cd ${o}/sequence/all_protein_with_header/
        for file in *.f*a; do mv -n "$file" "${file// /_}" 2>/dev/null; done    ## remove space in filename
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a  ## remove newline char
        sed -i '/^$/d' *.f*a ## remove empty lines
        gawk -i inplace 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' *.f*a ## uppercase
        sed -i '/^>/!s/J/X/g' *.f*a  ## avoid out-of-vocabulary
        cd ${Tool_Path}
    fi

    ### 3.4.processing genbank files
    if [[ $genbank == "TRUE" ]]; then
        printf "Processing genbank files...\n\n"
        # copy files
        for dir in ${i}*/
        do
            dir=${dir%*/}      # remove the trailing "/"
            cp ${dir}/genbank/*.gb* ${o}/sequence/genbank/
        done

    fi
fi


################### 4. Extract Features ###################
### feature extraction mode
export Genomic_FM=FALSE
export Homology=FALSE
export Genomic_Traits=FALSE
export Protein_FM=FALSE
export Protein_Traits=FALSE
export PPI=FALSE

######## default parameters for feature extraction
export Genome_MAX_LEN=512
export Genome_SIZE=3000
export Genome_OVERLAP=1500

export Protein_MAX_LEN=510
export Protein_SIZE=510
export Protein_OVERLAP=205

##########################################################

if [[ $e == "Genomic" ]]; then
    export Genomic_FM=TRUE
    export Homology=TRUE
    
elif [[ $e == "Genomic_Protein" ]]; then
    export Genomic_FM=TRUE
    export Homology=TRUE
    export Genomic_Traits=TRUE
    export Protein_FM=TRUE
    export Protein_Traits=TRUE
    
elif [[ $e == "Genomic_Protein_PPI" ]]; then
    export Genomic_FM=TRUE
    export Homology=TRUE
    export Genomic_Traits=TRUE
    export Protein_FM=TRUE
    export Protein_Traits=TRUE
    export PPI=TRUE
fi

### input types of files by given
if   [[ $f == "Genomic_FM" ]]; then
    export Genomic_FM=TRUE
elif [[ $f == "Protein_FM" ]]; then
    export Protein_FM=TRUE
elif [[ $f == "Genomic_Traits" ]]; then
    export Genomic_Traits=TRUE
elif [[ $f == "Protein_Traits" ]]; then
    export Protein_Traits=TRUE
elif [[ $f == "Homology" ]]; then
    export Homology=TRUE
elif [[ $f == "PPI" ]]; then
    export PPI=TRUE
fi

printf "Feature extraction methods: \n"
printf "By Genomic_FM: $Genomic_FM \n"
printf "By Protein_FM: $Protein_FM \n"
printf "By Homology: $Homology \n"
printf "By Genomic_Traits: $Genomic_Traits \n"
printf "By Protein_Traits: $Protein_Traits \n"
printf "By PPI: $PPI \n\n"

if [[ $Genomic_FM == 'FALSE' ]] && [[ $Protein_FM == 'FALSE' ]] && [[ $Homology == 'FALSE' ]] && [[ $Protein_Traits == 'FALSE' ]] && [[ $Genomic_Traits == 'FALSE' ]] && [[ $PPI == 'FALSE' ]];
then
    printf "No feature selected! need correct -f or -e parameters"
    exit 0; 
    
fi

###### process features
if [[ $w == 'only_predict' ]]; then
    printf "Skip feature extraction steps! \n"
fi

# ##### processing all features
if [[ $w == 'whole' ]]; then

    mkdir -p ${o}/features/
    for FEAT in Genomic_FM Protein_FM Homology Protein_Traits Genomic_Traits PPI
    do
        mkdir -p ${o}/features/${FEAT}
        rm -f ${o}/features/${FEAT}/* #remove old results
    done

    #### P1.Language Model Genome
    if   [[ ${Genomic_FM} == "TRUE" ]] ; then
        export FEAT=Genomic_FM

        cp ${o}/sequence/genome/*.f*a  ${o}/features/${FEAT}/
        cd ${o}/features/${FEAT}/
        sed -i '/>/d' *.f*a     ## remove header
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a  ## remove newline char
        for file in *.f*a; do sed -i "s/^/>${file%.fna}\n/g" "${file}"; done  ## add header
        for temp_f in *.f*a; do (cat "$temp_f" ; echo) >> all_sequences.fna  && rm "$temp_f" || break ; done ## merge into one

        cd ${Tool_Path}

        ### process features for Language model
        python Process/data_process.py \
            --mode Genome \
            --filepath ${o}/features/${FEAT}/all_sequences.fna \
            --out_folder ${o}/features/${FEAT}/ \
            --size $Genome_SIZE \
            --overlap $Genome_OVERLAP \

        printf "Finish processing ${FEAT} features"
        
    fi

    #### P2.Homology
    if   [[ ${Homology} == "TRUE" ]] ; then
        export FEAT=Homology

        cp ${o}/sequence/genome/*.f*a -t ${o}/features/${FEAT}/
        cd ${o}/features/${FEAT}/
        sed -i '/>/d' *.f*a     ## remove header
        gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } 
        /^>/ { print n $0; n = "" }
        END { printf "%s", n }
        ' *.f*a  ## remove newline char
        for file in *.f*a; do sed -i "s/^/>${file%.fna}\n/g" "${file}"; done  ## add header
        for temp_f in *.f*a; do (cat "$temp_f" ; echo) >> all_sequences.fna  && rm "$temp_f" || break ; done ## merge into one
        cd ${Tool_Path}

        ### run blast
        Query=${o}/features/${FEAT}/all_sequences.fna
        DB=$m/${FEAT}/viv_all.db
        Out=${o}/features/${FEAT}/blast_out.txt
        ### run
        blastn -task blastn -db $DB -query $Query -out $Out -outfmt 6  -num_threads $c  -max_hsps 1  -reward 2  -word_size 8  -gapopen 2 -gapextend 2 -evalue 0.001 -perc_identity 30 -qcov_hsp_perc 60

        printf "Finish processing ${FEAT} features \n"
        
    fi

    #### P3.Genome features
    if   [[ ${Genomic_Traits} == "TRUE" ]] ; then
        export FEAT=Genomic_Traits  
        

        python Process/Genome_feature.py \
            --assembly_path ${o}/sequence/genome/ \
            --filepath ${o}/sequence/genbank/ \
            --ref_cols Process/ref_cols.pickle \
            --out_folder ${o}/features/${FEAT}/ \
            --CPB_path External_Utils/CPB_Machine.jar \
            
        printf "Finish processing ${FEAT} features \n"
    fi


    #### P4.Language Model Protein
    if   [[ ${Protein_FM} == "TRUE" ]] ; then
        export FEAT=Protein_FM

        ### process features for Language model
        python Process/data_process.py \
            --mode Protein \
            --filepath ${o}/sequence/protein/ \
            --out_folder ${o}/features/${FEAT}/ \
            --size $Protein_SIZE \
            --overlap $Protein_OVERLAP \

        printf "Finish processing ${FEAT} features \n"

    fi

    #### P5.Protein features
    if   [[ ${Protein_Traits} == "TRUE" ]] ; then
        export FEAT=Protein_Traits

        python Process/Protein_feature.py \
            --assembly_path ${o}/sequence/protein/ \
            --filepath ${o}/sequence/all_protein/ \
            --out_folder ${o}/features/${FEAT}/ \
            --process $c \

        printf "Finish processing ${FEAT} features \n"
        
    fi

    #### P6.Protein PPI
    if   [[ ${PPI} == "TRUE" ]] ; then
        export FEAT=PPI

    ### part.1 HVIDB mapping
        printf "Start BLASTp for PPI features \n"
        cp ${o}/sequence/protein/*.f*a -t ${o}/features/${FEAT}/
        cd ${o}/features/${FEAT}/

        for temp_f in *.f*a; do (cat "$temp_f" ; echo) >> all_sequences.fna  && rm "$temp_f" || break ; done ## merge into one
        sed -i '/^$/d' all_sequences.fna ## remove empty lines
        cd ${Tool_Path}

        ### run blast
        Query=${o}/features/${FEAT}/all_sequences.fna
        DB=PPI_sources/HVIDB/hvidb_struct.db
        Out=${o}/features/${FEAT}/blast_out.txt

        blastp -task blastp -db $DB -query $Query -out $Out -outfmt 6  -num_threads $c  -word_size 6  -evalue 0.001 -qcov_hsp_perc 60

    #### part.2 extract SLiM
        printf "Extract SLiMs for PPI features \n"
        n | python2 External_Utils/slimsuite-SLiMSuite/tools/slimprob.py \
                    motifs=PPI_sources/elm2022_human.motifs \
                    batch=$o/sequence/all_protein_with_header/*.fasta \
                    dismask=T \
                    LogMask=F \
                    Pickle=F \
                    SaveSpace=1 \
                    disorder=iupred \
                    iucut=0.2 \
                    maxsize=10000000 \
                    maxseq=10000000 \
                    iupath=External_Utils/iupred/iupred \
                    iuchdir=T
        mv slimprob* -t $o/features/${FEAT}/
        cp ${o}/features/Protein_Traits/assembly_df.csv -t ${o}/features/${FEAT}/
        rm -r SLiMProb
        
    #### part.3 merge features  
        python Process/PPI_merge_features.py \
            --Blast_out_file ${o}/features/${FEAT}/blast_out.txt \
            --SLiM_folder ${o}/features/${FEAT}/ \
            --SLiM_all_motif PPI_sources/all_motif.pickle \
            --out_folder ${o}/features/${FEAT}/ \
            --Human_protein_to_GO PPI_sources/Human_protein_to_GO.csv \
            --HVPPI PPI_sources/HVPPI.csv \
            --Protein_rec_to_vid ${o}/features/${FEAT}/assembly_df.csv \
        
        printf "Finish processing ${FEAT} features \n"
        
    fi
    
fi 

################### 5. Run predictions ###################
export TASK=VIV_mhc_pos
export Genomic_FM_Tokenizer=Tokenizer/Genomic_FM/viro_bpe_32000.json
export Protein_FM_Tokenizer=Tokenizer/Protein_FM/ProtTrans
export Cols=Cates/categories.cols
export hierarchy=Cates/taxa_struct.json
export G_size=768
export G2L=384
export Support_NUM=5
export Num_host_categories=79
export Homology_REF=Predict/Reference_all.csv
export Cut_off=Predict/cut_off.pickle



mkdir -p ${o}/prediction/
for FEAT in Genomic_FM Protein_FM Homology Protein_Traits Genomic_Traits PPI
do
    mkdir -p ${o}/prediction/${FEAT}
    rm -f ${o}/prediction/${FEAT}/* #remove old results
done


#### Pre1.Language Model Genome
if   [[ ${Genomic_FM} == "TRUE" ]] ; then
    export FEAT=Genomic_FM
    CUDA_VISIBLE_DEVICES=$g python Predict/run_finetune_vbird.py \
        --output_dir ${o}/prediction/${FEAT} \
        --model_type_limit BertMHC \
        --model_name_or_path $m/${FEAT}/ \
        --max_length $Genome_MAX_LEN \
        --do_predict_no_label \
        --validation_file ${o}/features/${FEAT}/processed_input.csv \
        --task_name $TASK \
        --customize_tokenizer \
        --tokenizer_name $Genomic_FM_Tokenizer \
        --per_device_eval_batch_size $b \
        --preprocessing_num_workers $c \
        --logger_path Logger/${FEAT}.txt \
        --hierarchy_structure $hierarchy \
        --CLS_dropout 0.5 \
        --fp16 True \
        
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --hierarchy_structure $hierarchy \
        --categories ${Cols} \
        --cutoff ${Cut_off} \
        --agg_method topk \
        
fi
#### Pre2.Homology
if   [[ ${Homology} == "TRUE" ]] ; then
    export FEAT=Homology
    python Predict/PN_cal.py \
        --filepath ${o}/features/${FEAT}/ \
        --out_folder ${o}/prediction/${FEAT}/ \
        --reference_path ${Homology_REF} \
        --query_path ${o}/features/${FEAT}/all_sequences.fna \
        --support_num $Support_NUM \
        --process $c \
    
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --hierarchy_structure $hierarchy \
        --categories ${Cols} \
        --cutoff ${Cut_off} \
        --agg_method topk \
   
fi

#### Pre3.Genome features
if   [[ ${Genomic_Traits} == "TRUE" ]] ; then
    export FEAT=Genomic_Traits
    CUDA_VISIBLE_DEVICES=$g python Predict/run_train_features.py \
        --output_dir ${o}/prediction/${FEAT} \
        --feature_type Genomic_Traits \
        --label_num $Num_host_categories \
        --model_name_or_path $m/${FEAT}/ \
        --validation_file ${o}/features/${FEAT}/processed_input.csv \
        --do_predict_no_label \
        --hierarchy_structure $hierarchy \
        --g_size $G_size \
        --G2L $G2L \
        --per_device_eval_batch_size $b \
        --preprocessing_num_workers $c \
        --dropout_prob 0.1 \
        --fp16 True \
        
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --hierarchy_structure $hierarchy \
        --categories ${Cols} \
        --cutoff ${Cut_off} \
        --agg_method topk \
        
fi
    
#### Pre4.Language Model Protein
if   [[ ${Protein_FM} == "TRUE" ]] ; then
    export FEAT=Protein_FM
    CUDA_VISIBLE_DEVICES=$g python Predict/run_finetune_vbird.py \
        --output_dir ${o}/prediction/${FEAT} \
        --model_type_limit BertMHC \
        --model_name_or_path $m/${FEAT}/ \
        --max_length $Protein_MAX_LEN \
        --do_predict_no_label \
        --validation_file ${o}/features/${FEAT}/processed_input.csv \
        --task_name $TASK \
        --tokenizer_name $Protein_FM_Tokenizer \
        --per_device_eval_batch_size $b \
        --preprocessing_num_workers $c \
        --logger_path Logger/${FEAT}.txt \
        --hierarchy_structure $hierarchy \
        --CLS_dropout 0 \
        --fp16 True \
        
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --hierarchy_structure $hierarchy \
        --categories ${Cols} \
        --cutoff ${Cut_off} \
        --agg_method topk \
        
fi

### Pre5.Protein features
if   [[ ${Protein_Traits} == "TRUE" ]] ; then
    export FEAT=Protein_Traits
    CUDA_VISIBLE_DEVICES=$g python Predict/run_train_features.py \
        --output_dir ${o}/prediction/${FEAT}/ \
        --feature_type Protein_Traits \
        --label_num $Num_host_categories \
        --model_name_or_path $m/${FEAT}/ \
        --validation_file ${o}/features/${FEAT}/processed_input.csv \
        --do_predict_no_label \
        --hierarchy_structure $hierarchy \
        --g_size $G_size \
        --G2L $G2L \
        --per_device_eval_batch_size $b \
        --preprocessing_num_workers $c \
        --dropout_prob 0.1 \
        --fp16 True \
        
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --hierarchy_structure $hierarchy \
        --categories ${Cols} \
        --cutoff ${Cut_off} \
        --agg_method topk \
        
fi

# #### Pre6.Human Virus PPI
if   [[ ${PPI} == "TRUE" ]] ; then
    export FEAT=PPI
    CUDA_VISIBLE_DEVICES=$g python Predict/run_PPI.py \
        --feature_type PPI \
        --output_dir ${o}/prediction/${FEAT} \
        --model_name_or_path $m/${FEAT}/ \
        --do_predict_no_label \
        --validation_file ${o}/features/${FEAT}/processed_input.csv \
        --per_device_eval_batch_size $b \
        --preprocessing_num_workers $c \
        --logger_path Logger/${FEAT}.txt \
        --dropout_prob 0 \
        --fp16 True \
        
        
    python Predict/extract_prediction.py \
        --predict_folder ${o}/prediction/${FEAT}/ \
        --feature_folder ${o}/features/${FEAT}/ \
        --feature ${FEAT} \
        --categories Cates/only_human.cols \
        --cutoff ${Cut_off} \
        --agg_method topk \
        
fi



# ################### 6. Ensemble predictions ###################
### check
if [[ $e != @(Genomic|Genomic_Protein|Genomic_Protein_PPI) ]]; then
    printf "Finish prediction without using stacking models. \n"
    exit 0;
    
fi 


mkdir -p ${o}/ensemble
for Ensemble_mode in Genomic_Protein_PPI Genomic_Protein Genomic final
do
    mkdir -p ${o}/ensemble/$Ensemble_mode/
    rm -f ${o}/ensemble/$Ensemble_mode/* #remove old results
done



python Predict/aggregate_for_stacking.py \
    --predict_folder ${o}/prediction/ \
    --save_folder ${o}/ensemble/ \
    --ensemble_mode $e \
    
printf "Finish aggregating model probs. \n"
       

### predict


for Ensemble_mode in Genomic_Protein_PPI  Genomic_Protein  Genomic
do
    export Save_path=${o}/ensemble/$Ensemble_mode/
    export model_path=${Tool_Path}/$m/ensemble/$Ensemble_mode/LR_model.pkl
    if [[ -f ${o}/ensemble/$Ensemble_mode/test_prob_agg.pkl ]]; then
        if [[ $Ensemble_mode == "Genomic_Protein_PPI" ]]; then
            python Predict/StackingEnsemble.py \
                --trained_model $model_path \
                --save_path $Save_path \
                --individual_model_list MODEL_list.pkl \
                --valid_prob test_prob_agg.pkl \
                --LABELS LABELS.pkl \

        fi



        if [[ $Ensemble_mode == @(Genomic_Protein|Genomic) ]]; then
            python Predict/StackingEnsemble.py \
                --trained_model $model_path \
                --save_path $Save_path \
                --hierarchy_structure $hierarchy \
                --individual_model_list MODEL_list.pkl \
                --valid_prob test_prob_agg.pkl \
                --LABELS LABELS.pkl \

        fi
    fi
done

printf "Finish prediction of stacking models. \n"

export ensemble_path=${o}/ensemble/
export output_path=${o}/ensemble/final/
export all_cates=${Tool_Path}/Cates/categories.cols

python Predict/ensemble_extract_prediction.py \
    --ensemble_path $ensemble_path \
    --output_path $output_path \
    --all_cates $all_cates \
    --hierarchy_structure $hierarchy \
    --cutoff ${Cut_off} \
        
printf "Prediction completed! Final results stored in ${output_path}. \n"