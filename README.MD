# SeeVR on CBRAIN

Using SeeVR (https://github.com/abhogal-lab/seeVR) version 1.5.5.

## Useful commands:

Compilation:

    mcc -m Gas_CVR_wrapper.m -o gas_cvr -a ./code -d standalone

Docker build command:

    docker build -t dvalev/seevr-cbrain:0.1.0 .

Run Docker container interactively:

    docker run -it -v $PWD/data:/opt/seevr/data dvalev/seevr-cbrain:0.1.0

Command to run gas_cvr within container:

    ./run_gas_cvr.sh /opt/MCR-2021b/v911/ \[arguments\]

    ./run_gas_cvr.sh /opt/MCR-2021b/v911/ data time.txt CO2.txt BOLD_applytopup.nii.gz BOLD_mean_brain_seg_0.nii.gz BOLD_mean_brain_seg_1.nii.gz BOLD_mean_brain_seg_2.nii.gz BOLD_mean_brain_mask.nii.gz BOLD_mcf.par 0.3 5,10,15 1,350 97 3 5,25 3 0.6 -2,2 -3,25

Run in MATLAB:

    Gas_CVR_wrapper('data', 'time.txt', 'CO2.txt', 'BOLD_applytopup.nii.gz', 'BOLD_mean_brain_seg_0.nii.gz', 'BOLD_mean_brain_seg_1.nii.gz', 'BOLD_mean_brain_seg_2.nii.gz', 'BOLD_mean_brain_mask.nii.gz', 'BOLD_mcf.par', '0.3', '5,10,15', '1,350', '97', '3', '5,25', '3', '0.6', '-2,2', '-3,25')


