#!/bin/bash
source ~/.bashrc

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"

# define task
task=rest

# define directories
BIDS_dir="$path"/MAGMOT_BIDS
deriv_dir=$path/derivatives
ROI_dir="$deriv_dir"/ROI_masks/output
anal_dir=$deriv_dir/analysis/"$task"_paper
FC_dir=$anal_dir/FC
code_dir=$path/code/"$task"_paper/FC
data_dir=$anal_dir/data

# change directory to BIDS folder
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002)
#subjects=(sub-control001)
# sort array
subjects=($(echo ${subjects[*]}| tr " " "\n" | sort -n))

preproc=(s0 s4 s6 s8)
#preproc=(s4 s6 s8)

phases=(pre post online offline)
phases=(post online offline)

for smooth in "${preproc[@]}"; do

    # define ROI masks
    VTA_aHPC=$FC_dir/VTA_RSFC_"$smooth"_aHPC.nii.gz
    aHPC=$ROI_dir/MNI_res-epi_label-aHPC_mask.nii.gz

    # create file to save labels
    out_labels=labels.txt

    # for each subject in the subjects array
    for subject in "${subjects[@]}"; do

	    echo $subject

	    # create label: substitute group name with c / e respectively
	    if [[ "$subject" == *"cont"* ]]; then
            label="${subject/"sub-control0"/"c"}"
        else
			label="${subject/"sub-experimental0"/"e"}"
        fi

	    # print label into file: creates c01-c49 & e2-e50 (one value per row)
		if [[ "$subject" == *"sub-control001"* ]]; then
		   printf "$label" > "$code_dir"/$out_labels
		else
		    printf "\n$label" >> "$code_dir"/$out_labels
		fi

	    # go to directory with concat files
	    concat_dir=$data_dir/$subject
	    cd $concat_dir

	    # PRE-PROCESSED AND CONCATENATED RESTING FILES
        inputs=("$subject"_task-rest_run-1_desc-"$smooth"preproc_bold.nii.gz)
        inputs=("$subject"_task-rest_run-1_desc-"$smooth"preproc_bold.nii.gz "$subject"_task-rest_run-2_desc-"$smooth"preproc_bold.nii.gz "$subject"_task-online_desc-"$smooth"concat_bold.nii.gz "$subject"_task-offline_desc-"$smooth"concat_bold.nii.gz)
        inputs=("$subject"_task-rest_run-2_desc-"$smooth"preproc_bold.nii.gz "$subject"_task-online_desc-"$smooth"concat_bold.nii.gz "$subject"_task-offline_desc-"$smooth"concat_bold.nii.gz)


	    # loop over the preprocessed resting state files and extract average time course
	    for (( f=0; f<${#inputs[@]}; f++)); do

	    	# define file_id
	        file=${inputs[$f]}

            echo $file

            # define phase
            phase=${phases[$f]}

	        # create file to save correlation values
            out_gcor=gcor_"$phase"_"$smooth".txt

            # open file for correlation output
		    if [[ "$subject" == *"sub-control001"* ]]; then
                printf "BIDS\tgcor" > "$code_dir"/$out_gcor #file header
		    fi

		    # compute clobal correlation within GM mask
		    gcor=$(@compute_gcor -no_demean -verb 0 -input $file -mask $deriv_dir/$subject/func/"$subject"_task-rest_space-MNI152NLin2009cAsym_label-epiGM_mask.nii.gz)

            # add gcor to table
			printf "\n$label\t$gcor" >> "$code_dir"/$out_gcor # gcor for each run


	    done


    done

done
