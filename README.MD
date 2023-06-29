# SeeVR on CBRAIN

Using SeeVR (https://github.com/abhogal-lab/seeVR) version 1.5.5.

Current container version: dvalev/seevr-cbrain:0.1.0

## How to set up the container:

### Compile MATLAB code:

The MATLAB code for SeeVR + the Gas_CVR wrapper must be compiled using the MATLAB compiler toolkit add-on. From the root project directory, use this command in MATLAB to compile the code to the 'standalone' folder:

    mcc -m Gas_CVR_wrapper.m -o gas_cvr -a ./code -d standalone

### Build Docker container:

From the project root directory, run a Docker Build command to construct the image. It uses mcin/fsl (containing FSL) and dvalev/seevr-base (containing MATLAB runtime) as base layers.

Docker build command:

    docker build -t dvalev/seevr-cbrain:\[version tag] .

### Testing locally:

To test the container, you can run it locally in interactive mode and run commands inside it. Here is a sample command to start the container in interactive mode with the current working directory mounted to /opt/seevr:

Run Docker container interactively:

    docker run -it -v $PWD:/opt/seevr dvalev/seevr-cbrain:\[version tag]

Command to run preprocess_BIDS.sh in container:
    
    ./preprocess_BIDS.sh data/sub-006 outdir data/CO2.txt data/time.txt 0.3 5,10,15 1,350 97 3 5,25 3 0.6 -2,2 -3,25

Command to run gas_cvr within container:

    ./run_gas_cvr.sh /opt/MCR-2021b/v911/ \[arguments\]

    ./run_gas_cvr.sh /opt/MCR-2021b/v911/ data time.txt CO2.txt BOLD_applytopup.nii.gz BOLD_mean_brain_seg_0.nii.gz BOLD_mean_brain_seg_1.nii.gz BOLD_mean_brain_seg_2.nii.gz BOLD_mean_brain_mask.nii.gz BOLD_mcf.par 0.3 5,10,15 1,350 97 3 5,25 3 0.6 -2,2 -3,25

Run in MATLAB:

    Gas_CVR_wrapper('data', 'time.txt', 'CO2.txt', 'BOLD_applytopup.nii.gz', 'BOLD_mean_brain_seg_0.nii.gz', 'BOLD_mean_brain_seg_1.nii.gz', 'BOLD_mean_brain_seg_2.nii.gz', 'BOLD_mean_brain_mask.nii.gz', 'BOLD_mcf.par', '0.3', '5,10,15', '1,350', '97', '3', '5,25', '3', '0.6', '-2,2', '-3,25')
    
### Boutiques descriptor (instructions):

A boutiques descriptor needs to call the Docker image 'dvalev/seevr-cbrain:\[version]' and run in in Singularity. The entrypoint command for the container is:

    ./preprocess_BIDS.sh \[BIDS input] \[outdir] \[CO2 input] \[time input] \[motioncorr] \[dispersion] \[lim\_DoI] \[LVpercentile] \[spatialdim] \[normWindow] \[interp\_factor] \[corrthresh] \[lagthresh] \[-lag]
    
For example:

    ./preprocess_BIDS.sh data/sub-006 outdir data/CO2.txt data/time.txt 0.3 5,10,15 1,350 97 3 5,25 3 0.6 -2,2 -3,25
    
List type arguments can use the Boutiques "list" type to automatically take multiple values and format them as a comma-separated, spaceless string. These are parsed into arrays by Gas_CVR_wrapper().

There are no optional or non-positional arguments in this tool, so every argument must be either required or have a default value.
    



