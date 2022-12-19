#!/bin/bash
source ~/.bashrc

# change afni environment vars
@AfniEnv -set AFNI_ATLAS_COLORS MNI_Glasser_HCP_v1.0
@AfniEnv -set AFNI_IMAGE_LABEL_MODE 5
@AfniEnv -set AFNI_VALUE_LABEL YES
@AfniEnv -set AFNI_graph_ggap 0
@AfniEnv -set AFNI_ATLAS_LIST "CA_ML_18_MNI,MNI_Glasser_HCP_v1.0,Brainnetome_1.0,CA_MPM_22_MNI"


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

code_dir=$path/code/"$task"_paper/FC
data_dir=$anal_dir/data


############### specify masks and p value ###############

# specify mask (created in previous script)
epi_mask=sample_label-dilatedGM_mask.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define template
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`
sleep 10s

# inclusively mask the cluster map and the anatomical VTA ROI
VTASN=MNI_res-epi_label-VTASN_mask.nii.gz
aHPC=MNI_res-epi_label-aHPC_mask.nii.gz

# define seeds
seeds=(aHPC VTASN)


# define different smoothing kernels
preproc=(s0 s4 s6 s8)
#preproc=(s0)

# define groups
group1=control
group2=experimental

# define template
template=MNI152_2009_template_SSW.nii.gz


# define behavioural covariates
behav_covs=(highConf_abs CDMB)


############### create loop for FC analysis for different smoothing kernels ###############

for smooth in "${preproc[@]}"; do


cd "$FC_dir"/$smooth


  for seed in "${seeds[@]}"; do

    # define prefixes for ouput files

    #zval_all=report_FC_"$seed"_diff_"$smooth"_sample-all_type-zval.txt
    zval_diff=report_FC_"$seed"_diff_"$smooth"_sample-diff_type-zval.txt
    #zval_cont=report_FC_"$seed"_diff_"$smooth"_sample-cont_type-zval.txt
    #zval_exp=report_FC_"$seed"_diff_"$smooth"_sample-exp_type-zval.txt

    #ml_caez_all=report_FC_"$seed"_diff_"$smooth"_sample-all_type-CAEZ_ML.txt # ML = Macro Label
    ml_caez_diff=report_FC_"$seed"_diff_"$smooth"_sample-diff_type-CAEZ_ML.txt # ML = Macro Label
    #ml_caez_cont=report_FC_"$seed"_diff_"$smooth"_sample-cont_type-CAEZ_ML.txt # ML = Macro Label
    #ml_caez_exp=report_FC_"$seed"_diff_"$smooth"_sample-exp_type-CAEZ_ML.txt # ML = Macro Label

    #ml_hcp_all=report_FC_"$seed"_diff_"$smooth"_sample-all_type-HCP_ML.txt # ML = Macro Label
    ml_hcp_diff=report_FC_"$seed"_diff_"$smooth"_sample-diff_type-HCP_ML.txt # ML = Macro Label
    #ml_hcp_cont=report_FC_"$seed"_diff_"$smooth"_sample-cont_type-HCP_ML.txt # ML = Macro Label
    #ml_hcp_exp=report_FC_"$seed"_diff_"$smooth"_sample-exp_type-HCP_ML.txt # ML = Macro Label


    # define p and alpha value for thresholding
    p_val=0.001
    alpha_val=0.05

    # lenient threshold for exploratoty whole-brain analysis
    p_val=0.005
    k=5


    ### ONE SAMPLE t-TEST ACROSS WHOLE SAMPLE ###

    # in=1sample_"$seed"_diff_"$smooth"

    # extract corresponding cluster size k
    # k=$(1d_tool.py -infile "$in".CSimA.NN1_bisided.1D -csim_show_clustsize -csim_pthr $p_val -csim_alpha $alpha_val -verb 0)

    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates mask
    # 3dClusterize -inset "$in"+tlrc -ithr 1 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map mask.nii.gz > $zval_all # extract peak coordinate and binary mask

    # if any cluster survive thresholding
    # if [ -f "mask.nii.gz" ]; then

      # determine location of cluster: run report through atlas
      # whereami -omask mask.nii.gz -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_all # ML = macro label
      # whereami -omask mask.nii.gz -atlas CA_ML_18_MNI > $ml_caez_all # ML = macro label

      # rm mask.nii.gz

    # fi


    ### TWO SAMPLE t-TEST ###

    in=2sample_"$seed"_diff_"$smooth"

    # extract corresponding cluster size k
    # k=$(1d_tool.py -infile "$in".CSimA.NN1_bisided.1D -csim_show_clustsize -csim_pthr $p_val -csim_alpha $alpha_val -verb 0)

    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates mask
    3dClusterize -inset "$in"+tlrc -ithr 1 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map mask.nii.gz > $zval_diff # extract peak coordinate and binary mask

    # if any cluster survive thresholding
    if [ -f "mask.nii.gz" ]; then

      # determine location of cluster: run report through atlas
      whereami -omask mask.nii.gz -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_diff # ML = macro label
      whereami -omask mask.nii.gz -atlas CA_ML_18_MNI > $ml_caez_diff # ML = macro label

      # extract thresholded effect size map
      3dClusterize -inset "$in"+tlrc -ithr 3 -idat 2 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_dat effsize.nii.gz -pref_map bin.nii.gz -binary

      # create plot
      func_range=0.5
      fig=lenient_diff_"$seed"
      #fig=lenient_CDMB_aHPC.axi.jpg # debug
      op=8

      n=$(3dinfo -dmax mask.nii.gz)
      n=`expr $n + $n`


      @chauffeur_afni -ulay $template -olay effsize.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx $n -monty 1 -box_focus_slices bin.nii.gz -no_cor -no_sag

      rm mask.nii.gz
      rm effsize.nii.gz
      rm bin.nii.gz

    fi

    ### ANCOVA ###

    for ((c=0; c<${#behav_covs[@]}; c++)); do

      # determine covariate and column
      cov="${behav_covs[c]}"
      echo ""
      echo $cov
      echo ""

      in=2sample_"$cov"_"$seed"_diff_"$smooth"

      # debug
      # in=2sample_CDMB_aHPC_diff_"$smooth"

      zval_diff=report_FC_"$seed"_diff_"$smooth"_sample-"$cov"_type-zval.txt
      ml_caez_diff=report_FC_"$seed"_diff_"$smooth"_sample-"$cov"_type-CAEZ_ML.txt # ML = Macro Label
      ml_hcp_diff=report_FC_"$seed"_diff_"$smooth"_sample-"$cov"_type-HCP_ML.txt # ML = Macro Label

      # extract corresponding cluster size k
      # k=$(1d_tool.py -infile "$in".CSimA.NN1_bisided.1D -csim_show_clustsize -csim_pthr $p_val -csim_alpha $alpha_val -verb 0)

      # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates mask
      3dClusterize -inset "$in"+tlrc -ithr 3 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map mask.nii.gz > $zval_diff # extract peak coordinate and n-ary mask

      # if any cluster survive thresholding
      if [ -f "mask.nii.gz" ]; then

        # determine location of cluster: run report through atlas
        whereami -omask mask.nii.gz -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_diff # ML = macro label
        whereami -omask mask.nii.gz -atlas CA_ML_18_MNI > $ml_caez_diff # ML = macro label

        # extract thresholded effect size map
        3dClusterize -inset "$in"+tlrc -ithr 3 -idat 2 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_dat effsize.nii.gz -pref_map bin.nii.gz -binary

        # create plot
        func_range=0.5
        fig=lenient_"$cov"_"$seed"
        #fig=lenient_CDMB_aHPC.axi.jpg # debug
        op=8

        n=$(3dinfo -dmax mask.nii.gz)
        n=`expr $n + $n`


        @chauffeur_afni -ulay $template -olay effsize.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx $n -monty 1 -box_focus_slices bin.nii.gz -no_cor -no_sag

        rm mask.nii.gz
        rm effsize.nii.gz
        rm bin.nii.gz

      fi


    done

  done


done
