#!/bin/bash
source ~/.bashrc

################################################################################
# concatenating task BOLD series after pre-processing
################################################################################


# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# define derivatves dir
deriv_dir="$DIR"/derivatives

# change directory to BIDS folder
BIDS_dir="$DIR"/rawdata
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-control001)
#subjects=(sub-experimental016 sub-control045)

# define TSV file to read in
### input=$DIR/"code/task_paper/concat/MRC_inputForFinalConcatenation.tsv"
#online=$DIR"/code/rest_paper/preproc/MIC_inputForOnlineConcatenation.tsv"
#offline=$DIR"/code/rest_paper/preproc/MIC_inputForOfflineConcatenation.tsv"

phases=(online offline)

# define pre-proc
preproc=(s0 s4 s6 s8)


# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	# define concat directory
	outdir=$deriv_dir/analysis/rest_paper/data/$subject

	# load in pre-processed task file
	cd $outdir

	for smooth in "${preproc[@]}"; do

		############### concatenate TASK files ###############

		# define task
		task=magictrickwatching


		# loop over online and offline data
		for phase in "${phases[@]}"; do

		    # determine input file for 3dTcat
		    if [[ ${phase} = "online" ]]; then
			    input=$DIR"/code/rest_paper/preproc/MIC_inputForOnlineConcatenation.tsv"
			elif [[ ${phase} = "offline" ]]; then
			    input=$DIR"/code/rest_paper/preproc/MIC_inputForOfflineConcatenation.tsv"
			fi

		    # define searchstring and replace string
		    searchstring=task-"$task"_"desc-"$smooth"preproc_bold.nii.gz"
		    replacestring=task-"$phase"_"desc-"$smooth"concat_bold.nii.gz"

		    # define niifile
		    taskfile="$subject"_"$searchstring"

		    # set variable to track appearance of first magic trick
		    i=0
		    first=1

		    # read in file with information about scan duration
		    {
		    IGNORE_FL=T
		    while read ID stim_file start_vol end_vol
			    do

			    # if the line is header set IGNORE_FL to F and skip line
			    if [[ ${ID} = "ID" ]]; then
				    IGNORE_FL=F
					    continue
			    fi

			    # Ignore all lines until actual columns are found
			    if [[ ${IGNORE_FL} = T  ]]; then
				    continue
			    fi

			    # when the right information are available for the subject
			    if [[ "$ID" == *"$subject"* ]]; then

				    # update variable that tracks appearance of magictricks
				    ((i=i+1))

				    # when we're dealing with the first magic trick
				    if [[ $i -eq $first ]]; then

					    echo "processing concatenation"
					    # determine start and end vol for the first magic trick and save them in select
					    select=`echo $start_vol..$end_vol`

				    else

					    # determine start and end vol for the next magic trick and save them in select_temp
					    select_temp=`echo $start_vol..$end_vol`
					    #update select to include start and end of all magic tricks processed so far
					    select=`echo $select,$select_temp`

				    fi


			    fi

			    done < "$input"
			    }



		    # define prefix
		    prefix="${taskfile/$searchstring/$replacestring}"

		    # create task concat
		    if [ ! -f "$prefix" ]; then

			    echo $prefix

			    # once input file is processed
			    echo "selected vols: $select"

			    # 3dTcat to concatenate data using selected volumes
			    3dTcat $taskfile[$select] -prefix $prefix -session $outdir

			    # extract final file length
			    3dinfo $prefix > info.txt
			    numvol=$(awk 'NR==18 {print $6}' info.txt)
			    echo "final number of vols: $numvol"
			    rm info.txt

		    fi


		done






		############### concatenate REST files ###############

		# define task
		task=rest


		# define searchstring and replace string
		searchstring="desc-"$smooth"preproc_bold.nii.gz"
		replacestring1="run-1_desc-"$smooth"preproc_bold.nii.gz"
		replacestring2="run-2_desc-"$smooth"preproc_bold.nii.gz"

		# define niifile
		restfile="$subject"_task-"$task"_"$searchstring"

		# define tcat prefix
		prefix1="${restfile/$searchstring/$replacestring1}"
		prefix2="${restfile/$searchstring/$replacestring2}"

		# create run-1 (first 300 vols)
		if [ ! -f "$prefix1" ]; then
			echo $prefix1
			3dTcat -prefix $prefix1 $restfile[0..299]
		fi

		# create run-2 (second 300 vols)
		if [ ! -f "$prefix2" ]; then
			echo $prefix2
			3dTcat -prefix $prefix2 $restfile[300..599]
		fi

	done

done
