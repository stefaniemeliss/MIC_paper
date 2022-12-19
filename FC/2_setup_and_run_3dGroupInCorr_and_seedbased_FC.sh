#!/bin/bash
source ~/.bashrc

# this code runs a whole brain seed-based functional connectivity analysis using the HPC as seed to determine areas that exhibit functional connectivity

# the whole brain map is then intersecte with a VTA map to determine which voxel in the VTA show a significant FC with the HPC at pre learning rest


# define path
#path="/storage/shared/research/cinn/2018/MAGMOT"
path="/Users/nt807202/Dropbox/Reading/PhD/Magictricks/fmri_study/MMC"

# define task
task=rest

# define directories
deriv_dir=$path/derivatives
ROI_dir="$deriv_dir"/ROI_masks/output
anal_dir=$deriv_dir/analysis/"$task"_paper
CS_dir=$anal_dir/3dClustSim
FC_dir=$anal_dir/FC
mkdir $FC_dir
code_dir=$path/code/"$task"_paper/3_FC
data_dir=$anal_dir/data

############### specify seed masks and p value ###############

# specify mask (created in previous script)
epi_mask=sample_label-dilatedGM_mask.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define template
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`
sleep 10s

# define aHPC and VTA/SN
VTASN=MNI_res-epi_label-VTASN_mask.nii.gz
aHPC=MNI_res-epi_label-aHPC_mask.nii.gz

# define seeds
seeds=(aHPC VTASN)

############### other variables needed ###############

# define flow variables
run_setup=true
run_corr=true
run_data=true
run_diff=true
clustsim=true
#clustsim=false # debug
#run_setup=false
#run_corr=false
#run_data=false
#run_diff=false


# define different smoothing kernels
preproc=(s0 s4 s6 s8)
#preproc=(s4) # debug
#preproc=(s0) # debug

# define prefix
task="rest"

# define phases
phases=(pre post online offline)
#phases=(pre post) # debug

# define groups
group1=control
group2=experimental

# define behavioural covariates
behav_covs=(highConf_abs CDMB)
num_covs="${#behav_covs[@]}"


############### create loop for FC analysis for different smoothing kernels ###############

for smooth in "${preproc[@]}"; do

  # create folder for each kernel
  mkdir $FC_dir/$smooth

  # define input file pattern
  inputs=(*_task-rest_run-1_desc-"$smooth"preproc_bold.nii.gz *_task-rest_run-2_desc-"$smooth"preproc_bold.nii.gz *_task-online_desc-"$smooth"concat_bold.nii.gz *_task-offline_desc-"$smooth"concat_bold.nii.gz)
  inputs=(*_task-rest_run-1_desc-"$smooth"preproc_bold.nii.gz *_task-rest_run-2_desc-"$smooth"preproc_bold.nii.gz) # debug

  # loop over the preprocessed resting state files and extract average time course
  for (( f=0; f<${#inputs[@]}; f++)); do

    cd "$data_dir" # change directory to where pre-processed files are

    # define file pattern
    file=${inputs[$f]}
    echo $file

    # define phase
    phase=${phases[$f]}

    if ($run_setup == true) then

      ############### set up the file used to compute correlations (pearson) ###############

      # 3dSetupGroupInCorr: Summarise data for 3dGroupInCorr in an object
      out_GroupInCorr=FC_"$phase"_"$smooth"
      3dSetupGroupInCorr -mask $ROI_dir/$gm_mask -prefix "$FC_dir"/$preproc/$out_GroupInCorr -byte -labels "$code_dir"/labels.txt ./sub-*/$file

    fi

    if ($run_corr == true) then

      ############### compute functional connectivity using HPC as seed ###############

      cd "$FC_dir"/$smooth

      # copy clustsim output
      cp $CS_dir/ClustSim_"$task"_"$smooth"*  .

      # define phase
      phase=${phases[$f]}
      echo $phase

      # create seeds_FC file
      printf "seedFC_aHPC_"$phase"_"$smooth" "$ROI_dir"/""$aHPC" > "$code_dir"/seeds_FC # aHPC mask
      printf "\nseedFC_VTASN_"$phase"_"$smooth" "$ROI_dir"/""$VTASN" >> "$code_dir"/seeds_FC # VTA/SN mask

      # 3dGroupInCorr: Use data object created using 3dSetupGroupInCorr to  Summarise data for 3dGroupInCorr in an object
      out_GroupInCorr=FC_"$phase"_"$smooth"
      3dGroupInCorr -setA "$FC_dir"/$preproc/$out_GroupInCorr.grpincorr.niml -labelA "$phase""_all" -verb -clust ClustSim_"$task"_"$smooth" -covariates "$code_dir"/gcor_"$phase"_"$smooth".txt -donocov -sendall -batch MASKAVE "$code_dir"/seeds_FC

      # remove seeds_FC file
      rm "$code_dir"/seeds_FC

      # define aHPC seed map
      aHPC_map=seedFC_aHPC_"$phase"_"$smooth"+tlrc

      # extract ROI average time course: for each volume, the average of voxel is computed in saved in txt file
      3dmaskave -quiet -mask $ROI_dir"/"$VTASN $aHPC_map > maskave_VTASN_"$phase"_"$smooth" # FC values in VTA (determined with anterior HPC)

      # the output file will contain the average value for the sample in the first 6 bricks (mean, z_mean, gcor, z_cor, mean_NC, mean_NC_z)
      # as well as as the average values for sub-bricks 7:56 (data for each subject)

      # remove the Clustsim output from directory
      rm $FC_dir/$smooth/ClustSim_"$task"_"$smooth"*

    fi


    if ($run_data == true) then

      ############### extract pre and post data for each subject ###############

      if [[ "$phase" == *"p"* ]]; then

        cd "$FC_dir"/$smooth

        for seed in "${seeds[@]}"; do

          # define aHPC and VTA seed map
          seed_map=seedFC_"$seed"_"$phase"_"$smooth"+tlrc

          end=$(3dinfo -nvi $seed_map) # total number of volumes

          for i in $(seq 0 $end); do

            # determine name of sub-brick label
            subbrick_label=$(3dinfo -label "$seed_map".[$i])
            echo $subbrick_label

            # check whether subbrick belongs to a single subject (those start with A_)
            if [[ "$subbrick_label" == *"A_"* ]]; then

              # create prefix
              searchstr="A"
              replacestr="$seed""_""$phase""_""$smooth"
              replacestr="$seed""_""$phase"
              subbrick_new="${subbrick_label/$searchstr/$replacestr}"

              # extract sub-brick and save as new file
              3dbucket -prefix "$subbrick_new".nii.gz "$seed_map".[$i]

            fi

          done

        done

      fi

    fi

  done

  if ($run_diff == true) then

    cd "$FC_dir"/$smooth

    ############### compute difference image for each seed and smoothing kernel ###############

    for seed in "${seeds[@]}"; do

      # define pre and post files
      pre_files=($(ls "$seed"_pre_"$smooth"_*_zcorr.nii.gz))
      post_files=($(ls "$seed"_post_"$smooth"_*_zcorr.nii.gz))
      pre_files=($(ls "$seed"_pre_*_zcorr.nii.gz))
      post_files=($(ls "$seed"_post_*_zcorr.nii.gz))

      # sort arrays
      pre_files=($(echo ${pre_files[*]}| tr " " "\n" | sort -n))
      post_files=($(echo ${post_files[*]}| tr " " "\n" | sort -n))


      num_files="${#post_files[@]}"
      echo "${#pre_files[@]}"
      echo "${#post_files[@]}"

      for ((j=0; j<$num_files; j++)); do

        # determine pre- and post file
        pre_file="${pre_files[j]}"
        post_file="${post_files[j]}"

        echo $pre_file
        echo $post_file

        # create prefix
        searchstr="post"
        replacestr="diff"
        diff_prefix="${post_file/$searchstr/$replacestr}"

        # create diff file
        if [ ! -f "$diff_prefix" ]; then

          echo $diff_prefix

          # calculate difference
          3dcalc -a $post_file -b $pre_file -expr 'a-b' -prefix $diff_prefix

        fi

      done

      ############### T tests / ANCOVA on difference maps ###############

      # run 3dttest++ with Clustsim for thresholding
      if ($clustsim); then

        # two sample t-test as well as one sample within each group
        #3dttest++ -setA "$seed"_diff_"$smooth"_c*_zcorr.nii.gz  -setB "$seed"_diff_"$smooth"_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_"$seed"_diff_"$smooth" -ClustSim
        3dttest++ -setA "$seed"_diff_c*_zcorr.nii.gz  -setB "$seed"_diff_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_"$seed"_diff_"$smooth" -ClustSim

        # one sample t-test across whole sample
        #3dttest++ -setA "$seed"_diff_"$smooth"_*_zcorr.nii.gz -labelA "whole_sample" -mask $ROI_dir/$gm_mask -prefix 1sample_"$seed"_diff_"$smooth" -ClustSim
        # 3dttest++ -setA "$seed"_diff_*_zcorr.nii.gz -labelA "whole_sample" -mask $ROI_dir/$gm_mask -prefix 1sample_"$seed"_diff_"$smooth" -ClustSim

        for ((c=0; c<$num_covs; c++)); do

          # determine covariate and column
          cov="${behav_covs[c]}"
          echo ""
          echo $cov
          echo ""

          # run ANCOVA
          3dttest++ -setA "$seed"_diff_c*_zcorr.nii.gz  -setB "$seed"_diff_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -mask $ROI_dir/$gm_mask -prefix 2sample_"$cov"_"$seed"_diff_"$smooth" -covariates "$code_dir"/"$cov"_"$seed".txt -ClustSim

        done

        # remove CSim num_files
        #rm *CSim*

      else

        # two sample t-test as well as one sample within each group
        #3dttest++ -setA "$seed"_diff_"$smooth"_c*_zcorr.nii.gz  -setB "$seed"_diff_"$smooth"_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_"$seed"_diff_"$smooth"
        3dttest++ -setA "$seed"_diff_c*_zcorr.nii.gz  -setB "$seed"_diff_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_"$seed"_diff_"$smooth"

        # one sample t-test across whole sample
        #3dttest++ -setA "$seed"_diff_"$smooth"_*_zcorr.nii.gz -labelA "whole_sample" -mask $ROI_dir/$gm_mask -prefix 1sample_"$seed"_diff_"$smooth"
        # 3dttest++ -setA "$seed"_diff_*_zcorr.nii.gz -labelA "whole_sample" -mask $ROI_dir/$gm_mask -prefix 1sample_"$seed"_diff_"$smooth"

        for ((c=0; c<${#num_covs[@]}; c++)); do

          # determine covariate and column
          cov="${behav_covs[c]}"
          echo ""
          echo $cov
          echo ""

          # run ANCOVA
          3dttest++ -setA "$seed"_diff_c*_zcorr.nii.gz  -setB "$seed"_diff_e*_zcorr.nii.gz -labelA $group1 -labelB $group2 -mask $ROI_dir/$gm_mask -prefix 2sample_"$cov"_"$seed"_diff_"$smooth" -covariates "$code_dir"/"$cov"_"$seed".txt -center SAME

        done

      fi

      # delete all single subject files
      rm "$seed"_pre_*_zcorr.nii.gz
      rm "$seed"_post_*_zcorr.nii.gz
      rm "$seed"_diff_*_zcorr.nii.gz

    done

  fi

  cp $template_path/$template ./$template # copy template to directory


done

# code to create nii images for NeuroVault #
cd "$FC_dir"/s4

3dAFNItoNIFTI -prefix seedFC_aHPC_pre.nii.gz seedFC_aHPC_pre_s4+tlrc.
3dAFNItoNIFTI -prefix seedFC_aHPC_post.nii.gz seedFC_aHPC_post_s4+tlrc.
3dAFNItoNIFTI -prefix 2sample_aHPC_diff.nii.gz 2sample_aHPC_diff_s4+tlrc.
3dAFNItoNIFTI -prefix ANCOVA_aHPC_CDMB.nii.gz 2sample_CDMB_aHPC_diff_s4+tlrc.BRIK
3dAFNItoNIFTI -prefix ANCOVA_aHPC_total.nii.gz 2sample_highConf_abs_aHPC_diff_s4+tlrc.BRIK

3dAFNItoNIFTI -prefix seedFC_VTASN_pre.nii.gz seedFC_VTASN_pre_s4+tlrc.
3dAFNItoNIFTI -prefix seedFC_VTASN_post.nii.gz seedFC_VTASN_post_s4+tlrc.
3dAFNItoNIFTI -prefix 2sample_VTASN_diff.nii.gz 2sample_VTASN_diff_s4+tlrc.
3dAFNItoNIFTI -prefix ANCOVA_VTASN_CDMB.nii.gz 2sample_CDMB_VTASN_diff_s4+tlrc.BRIK
3dAFNItoNIFTI -prefix ANCOVA_VTASN_total.nii.gz 2sample_highConf_abs_VTASN_diff_s4+tlrc.BRIK
