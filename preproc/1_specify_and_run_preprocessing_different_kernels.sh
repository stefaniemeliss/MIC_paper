#!/bin/tcsh
source ~/.cshrc

################################################################################
# pre-processing for functional data using anfi_proc.py 
################################################################################

#module load afni19.3.03

# --------------------------------------------------------------------
# Script: s.2016_ChenEtal_02_ap.tcsh
#
# From:
# Chen GC, Taylor PA, Shin Y-W, Reynolds RC, Cox RW (2016). Untangling
# the Relatedness among Correlations, Part II: Inter-Subject
# Correlation Group Analysis through Linear Mixed-Effects
# Modeling. Neuroimage (in press).
#
# Originally run using: AFNI_16.1.16
# --------------------------------------------------------------------

# further modified based on Example 11. (see afni_proc.py -help)

# FMRI processing script, ISC movie data.
# Assumes previously run FS and SUMA commands, respectively: 
# $ recon-all -all -subject $subj -i $anat
# $ @SUMA_Make_Spec_FS -sid $subj -NIFTI

# Set top level directory structure
set topdir = /storage/shared/research/cinn/2018/MAGMOT #study folder
echo $topdir

set derivroot = $topdir/derivatives
set outroot = $derivroot/afni_proc

# define subject listecho $
set BIDSdir = $topdir/rawdata

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-experimental04*`) # this creates an array containing all subjects in the BIDS directory
#echo $subjects
echo $#subjects

#set subjects	= (sub-experimental016) 
#set subjects	= sub-control002
#set subjects	= sub-control017

# set different smoothing kernels
set kernels = (0 4 6 8)
#set kernels = (0)

# define tasks to loop through
set tasks = (magictrickwatching rest)
#set tasks = (rest)
#set tasks = (magictrickwatching)



# for each task
foreach task ($tasks)

    # for each subject in the subjects array
    foreach subject ($subjects)
	    
	    ############################# SET UP PRE-PROCESSING ############################# 
	    
	    #set subject	= "sub-experimental005"
	    echo $subject
	    
	    set subj = "$subject"_task-"$task"
	    
	    # define output dir of original afni_proc.py command
	    set output_dir = "$subj".results
	    cd $outroot/$subject/$output_dir
	    
	    # make output folder
        set anal_out = $derivroot/analysis/rest_paper/data/$subject
        mkdir -p $anal_out
	    
	    # define number of runs
	    if ($task == magictrickwatching) then
	        if ($subject == sub-experimental016) then
		        set runs = (`count -digits 2 1 4`)
	        else
		        set runs = (`count -digits 2 1 3`)
	        endif
	    else 
	        set runs = (`count -digits 2 1 2`)
	    endif
	    

	    ############################# RUN RE-PROCESSING ############################# 
		    
	    foreach blur ($kernels)

	        set subjstr = "${subj}"_s"$blur"
	        

	        if ($blur != 0) then

			    # ================================== blur ==================================
			    # blur each volume of each run
			    foreach run ( $runs )
				    3dBlurToFWHM -FWHM $blur -mask mask_epi_anat.$subj+tlrc \
						         -input pb03.$subj.r$run.volreg+tlrc    \
						         -prefix pb00.$subjstr.r$run.blur 
			    end
			        
		    endif
		    
		    
		    # ================================= scale ==================================
		    # scale each voxel time series to have a mean of 100
		    # (be sure no negatives creep in)
		    # (subject to a range of [0,200])
		    foreach run ( $runs )
		    

		        if ($blur == 0) then
		            # create average time course based on either smoothed or unsmoothed data		        
		            3dTstat -prefix rm.mean_r$run pb03.$subj.r$run.volreg+tlrc		      
		            
		            # scale
			        3dcalc -a pb03.$subj.r$run.volreg+tlrc -b rm.mean_r$run+tlrc \
			               -c mask_epi_anat.$subj+tlrc                           \
			               -expr 'c * min(200, a/b*100)*step(a)*step(b)'         \
			               -prefix pb00.$subjstr.r$run.scale
      
		        else
		            # create average time course based on either smoothed or unsmoothed data
		            3dTstat -prefix rm.mean_r$run pb00.$subjstr.r$run.blur+tlrc
		            
		            # scale
			        3dcalc -a pb00.$subjstr.r$run.blur+tlrc -b rm.mean_r$run+tlrc \
			               -c mask_epi_anat.$subj+tlrc                           \
			               -expr 'c * min(200, a/b*100)*step(a)*step(b)'         \
			               -prefix pb00.$subjstr.r$run.scale

		        endif
		    
		    end

		    # ================================ regress =================================
		    # ------------------------------
		    # run the regression analysis
		    
		    if ($task == magictrickwatching) then

	            if ($subject == sub-experimental016) then

			        3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                      \
				        -mask mask_epi_anat.$subj+tlrc                                        \
				        -censor censor_${subj}_combined_2.1D                                  \
				        -ortvec bandpass_rall.1D bandpass                                     \
				        -ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				        -ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				        -ortvec ROIPC.FSvent.r03.1D ROIPC.FSvent.r03                          \
				        -ortvec ROIPC.FSvent.r04.1D ROIPC.FSvent.r04                          \
				        -ortvec mot_demean.r01.1D mot_demean_r01                              \
				        -ortvec mot_demean.r02.1D mot_demean_r02                              \
				        -ortvec mot_demean.r03.1D mot_demean_r03                              \
				        -ortvec mot_demean.r04.1D mot_demean_r04                              \
				        -ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				        -ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				        -ortvec mot_deriv.r03.1D mot_deriv_r03                                \
				        -ortvec mot_deriv.r04.1D mot_deriv_r04                                \
				        -polort 2 -float                                                      \
				        -num_stimts 0                                                         \
				        -fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				        -x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				        -fitts fitts.$subjstr                                                 \
				        -errts errts.${subjstr}                                               \
				        -x1D_stop                                                             \
				        -bucket stats.$subjstr

		        else

			        3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                      \
				        -mask mask_epi_anat.$subj+tlrc                                        \
				        -censor censor_${subj}_combined_2.1D                                  \
				        -ortvec bandpass_rall.1D bandpass                                     \
				        -ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				        -ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				        -ortvec ROIPC.FSvent.r03.1D ROIPC.FSvent.r03                          \
				        -ortvec mot_demean.r01.1D mot_demean_r01                              \
				        -ortvec mot_demean.r02.1D mot_demean_r02                              \
				        -ortvec mot_demean.r03.1D mot_demean_r03                              \
				        -ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				        -ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				        -ortvec mot_deriv.r03.1D mot_deriv_r03                                \
				        -polort 2 -float                                                      \
				        -num_stimts 0                                                         \
				        -fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				        -x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				        -fitts fitts.$subjstr                                                 \
				        -errts errts.${subjstr}                                               \
				        -x1D_stop                                                             \
				        -bucket stats.$subjstr

		        endif
		      
		    else
		    
		        # resting state
		        3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                      \
                    -mask mask_epi_anat.$subj+tlrc                                        \
                    -censor censor_${subj}_combined_2.1D                                  \
                    -ortvec bandpass_rall.1D bandpass                                     \
                    -ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
                    -ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
                    -ortvec mot_demean.r01.1D mot_demean_r01                              \
                    -ortvec mot_demean.r02.1D mot_demean_r02                              \
                    -ortvec mot_deriv.r01.1D mot_deriv_r01                                \
                    -ortvec mot_deriv.r02.1D mot_deriv_r02                                \
                    -polort 2 -float                                                      \
                    -num_stimts 0                                                         \
                    -fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
                    -x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
                    -fitts fitts.$subjstr                                                 \
                    -errts errts.${subjstr}                                               \
                    -x1D_stop                                                             \
                    -bucket stats.$subjstr
                    
            endif
		    
		    # -- use 3dTproject to project out regression matrix --
		    #    (make errts like 3dDeconvolve, but more quickly)
		    3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				       -mask mask_epi_anat.$subj+tlrc                                 \
				       -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				       -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.tproject

		    # if 3dDeconvolve fails, terminate the script
		    if ( $status != 0 ) then
			    echo '---------------------------------------'
			    echo '** 3dDeconvolve error, failing...'
			    echo '   (consider the file 3dDeconvolve.err)'
			    exit
		    endif

		    # ============================ local WM regres =============================		
		    # -- use 3dTproject to project out regression matrix --
		    #    (make errts like 3dDeconvolve, but more quickly)
		    3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				       -mask mask_epi_anat.$subj+tlrc                                 \
				       -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				       -dsort Local_FSWMe_rall+tlrc                                   \
				       -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.fanaticor
				       
		    # ============================ blur estimation =============================
		    
            # compute blur estimates
            touch blur_est.$subjstr.1D   # start with empty file

            # create directory for ACF curve files
            mkdir files_ACF_s$blur

            # -- estimate blur for each run in errts --
            touch blur.errts.$subjstr.1D
            
            # restrict to uncensored TRs, per run
            foreach run ( $runs )
                set trs = `1d_tool.py -infile X.s${blur}.xmat.1D -show_trs_uncensored encoded                   \
                                      -show_trs_run $run`
                if ( $trs == "" ) continue
                3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                                        \
                        -ACF files_ACF_s$blur/out.3dFWHMx.ACF.errts.r$run.1D                                  \
                        errts.$subjstr.fanaticor+tlrc"[$trs]" >> blur.errts.$subjstr.1D

            end
            
            # compute average FWHM blur (from every other row) and append
            set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{0..$(2)}'\'` )
            echo average errts FWHM blurs: $blurs
            echo "$blurs   # errts FWHM blur estimates" >> blur_est.$subjstr.1D

            # compute average ACF blur (from every other row) and append
            set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{1..$(2)}'\'` )
            echo average errts ACF blurs: $blurs
            echo "$blurs   # errts ACF blur estimates" >> blur_est.$subjstr.1D

		    # ================================ convert =================================

		    # define niifile
		    set niifile = $subj"_desc-s"$blur"preproc_bold.nii.gz"

		    # do the AFNI to .nii conversion
	        3dAFNItoNIFTI -prefix $anal_out/$niifile errts.${subjstr}.fanaticor+tlrc

		    # remove temporary files
		    \rm -f rm.*
		    
	    end # end blur

    end # end subject
    
end # end task
