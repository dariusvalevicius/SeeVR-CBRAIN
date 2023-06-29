function Gas_CVR_wrapper(dir, time_CO2, data_CO2, mdcBOLD, GMmask, WMmask, CSFmask, mbBOLDmask, mpfilename, motioncorr, dispersion, lim_DoI, LVpercentile, spatialdim, normWindow, interp_factor, corrthresh, lagthresh,lag)
%%
% Inputs:
   % - dir = directory of all MRI data
   % - time_CO2 = time points of end-tidal CO2 data
   % - data_CO2 = end-tidal CO2 data
   % - mdcBOLD = motion and distortion corrected BOLD series (output
   %             of MCFLIRT and TOPUP pacages from FSL)
   % - GMmask = Grey matter mask
   % - WMmask = white matter mask
   % - CSFmask = CSF mask
   % - mbBOLDmask = mean BOLD brain mask
   % - mpfilename = Find MCFLIRT motion parameter file
   % - motioncorr = motion correction correlation threshold i.e. 0.3
   % - dispersion = dispersion for calculating HRF_probe_CO2 using
   %                convHRF function i.e. [5 10 15 25]
   % - lim_DoI = indices of the vector data to chop data of interest
   %             using chopTimeseries function i.e. [1 390]
   % - LVpercentile = percentile to remove Large vessel effect
   % - spatialdim = spatial dimension for 3D gaussian smoothing i.e. 2
   % - normWindow = normalization window for data i.e. [5 25]
   % - interp_factor = temporally interpolation factor i.e. 2
   % - corrthresh = correlation threshold for generating regressor default is 0.7
   % - lagthresh = Thresholds for refinening optimized regressors
   %               [lowlagthresh highlagthresh] i.e. [-2 2]
   %  -lag = Lag thresholds (in units of TR) for lag map creation
   %         [lowerlg upperlag] i.e. [-3 25]

dispersion = str2double(strsplit(dispersion, ','));
lim_DoI = str2double(strsplit(lim_DoI, ','));
normWindow = str2double(strsplit(normWindow, ','));
corrthresh = str2double(strsplit(corrthresh, ','));
lagthresh = str2double(strsplit(lagthresh, ','));
lag = str2double(strsplit(lag, ','));

corrthresh = str2double(corrthresh);
interp_factor = str2double(interp_factor);
LVpercentile = str2double(LVpercentile);
motioncorr = str2double(motioncorr);
spatialdim = str2double(spatialdim);


%%
% addpath(genpath('/NAS/home/al_reza/CBRAIN/1-seeVR/Resource/seeVR-main'))

dir = [pwd '/' dir '/']; % I assume that the current folder is the folder including all files

% initialize the opts structure (!)
global opts
opts.verbose = 0;
% % Set the location for the end-tidal data
% opts.seqpath = gas_dir;

cd(dir)
% Load motion corrected data
[BOLD,INFO] = loadTimeseries(dir, mdcBOLD);

% Load GM segmentation
[GMmask,INFOmask] = loadMask(dir, GMmask);

% Load WM segmentation
[WMmask,~] = loadMask(dir, WMmask);

% Load CSF segmentation
[CSFmask,~] = loadMask(dir, CSFmask);

% Load brain mask
[WBmask,~] = loadMask(dir, mbBOLDmask);

%% 2) Setup directories for saving output
  
% specify the root directory to save results
opts.savedir = dir; 

% specify a sub directory to save parameter maps - this
% directory can be changed for multiple runs etc.
if ispc
opts.resultsdir = [dir,'RESULTS\']; mkdir(opts.resultsdir);
else
opts.resultsdir = [dir,'RESULTS/']; mkdir(opts.resultsdir);    
end

% specify a sub directory to save certain figures - this
% directory can be changed for multiple runs etc.
if ispc
opts.figdir = [dir,'FIGURES\']; mkdir(opts.figdir);
else
opts.figdir = [dir,'FIGURES/']; mkdir(opts.figdir);    
end

%% 3) Load end-tidal gas traces and temporally align with MRI data
% Use the GM time-series for alignment
TS = meanTimeseries(BOLD,GMmask);

% resample input gas data to TR
[resamp_time_CO2, CO2_corr] = resampletoTR_new(opts.TR,time_CO2,data_CO2); 
% [resamp_time_O2, O2_corr] = resampletoTR_new(opts.TR,time_O2,data_O2);

% Automaticaaly align the resampled to TR data to Timepoints of GM BOLD data
[probe1a,probe2a] = autoAlign(CO2_corr, TS); % for this example the alignment offset was 7
% If both CO2 & O2 are supplied, the GUI returns two probe vectors to the
% Matlab workspace. If only CO2 or O2 are provided, then 1 probe is
% returned.
CO2trace = probe1a;
% O2trace = probe2a;
clear TS probe1a probe2a

% This part is useful if we want to give some options to user whioxch is
% not feasible in CBRAIN
% if isequal(gen,'Respiract 4th Generation')
%     % Load breathing data based on RespirAct system
%     [CO2_corr, O2_corr] = loadRAMRgen4(opts); 
% 
%     % Correlation with the time-series is done using the first input.
%     [probe1a,probe2a] = autoAlign(CO2_corr, O2_corr, TS); % for this example the alignment offset was 7
%     % If both CO2 & O2 are supplied, the GUI returns two probe vectors to the
%     % Matlab workspace. If only CO2 or O2 are provided, then 1 probe is
%     % returned.
%     CO2trace = probe1a;
%     O2trace = probe2a;
%     clear TS probe1a probe2a
%     
% elseif isequal(gen, 'Respiract 3rd Generation')
%     [corrvec_CO2, corrvec_O2] = loadRAMRgen3(opts);
% 
%     % Correlation with the time-series is done using the first input.
%     [probe1a,probe2a] = autoAlign(corrvec_CO2, corrvec_O2, TS); % for this example the alignment offset was 7
%     % If both CO2 & O2 are supplied, the GUI returns two probe vectors to the
%     % Matlab workspace. If only CO2 or O2 are provided, then 1 probe is
%     % returned.
%     CO2trace = probe1a;
%     O2trace = probe2a;
%     clear TS probe1a probe2a
% elseif isequal(gen, '') % arbitrary end tidal gas data 
%     [resamp_time_CO2, CO2_corr] = resampletoTR(opts.TR,time_CO2,data_CO2);
%     [resamp_time_O2, O2_corr] = resampletoTR(opts.TR,time_O2,data_O2);
%     
%     % Correlation with the time-series is done using the first input.
%     [probe1a,probe2a] = autoAlign(CO2_corr, O2_corr, TS); % for this example the alignment offset was 7
%     % If both CO2 & O2 are supplied, the GUI returns two probe vectors to the
%     % Matlab workspace. If only CO2 or O2 are provided, then 1 probe is
%     % returned.
%     CO2trace = probe1a;
%     O2trace = probe2a;
%     clear TS probe1a probe2a
% end


%% 4) Regress out motion parameters
cd(dir)  % Go to our data directory
nuisance = load(mpfilename); % Load nuisance regressors (translation, rotation etc.)

% Calculate motion derivatives
dtnuisance =  gradient(nuisance); % temporal derivative
sqnuisance = nuisance.*nuisance; % square of motion
motion = [nuisance dtnuisance sqnuisance];

% Demean/rescale regressors
for ii=1:(size(motion,2)); motion(:,ii) = demeanData(rescale(motion(:,ii))); end

% Filter regressors
reference = meanTimeseries(BOLD,GMmask);
opts.motioncorr = motioncorr; % correlation (r) threshold
[keep, leave] = filtRegressor(motion, reference, opts);

% Attempt to clean all voxels contained in the WB mask using nuisance
% regressors. In this example, there are not many motion related artifacts.
opts.disp = dispersion; %dispersion
opts.plot = 0; %shows hrf
[~,~,HRF_probe_CO2] = convHRF(CO2trace, opts);

%regression
% add a linear term to acount for signal drift


% Data scrubbing
[cleanData] = scrubData(BOLD,WBmask, keep, HRF_probe_CO2, opts);

% You can also consider to include the nuisance probes with high
% correlation to the reference signal as data probes or a linear term
%L = rescale(LegendreN(1,xdata));
%[cleanData] = scrubData(BOLD,WBmask, [L' keep], [HRF_probe_CO2' leave], opts);


%% 5) Isolate data of interest
% Similar to the normTimeseries function, chopTimeseries can be supplied
% with a 2 element vector [idx1 idx2]. If no indicies are provided, manual
% input will be requested (try it to isolate different data)
[idx, rBOLD] = chopTimeseries(cleanData, GMmask, lim_DoI); 

%establish new xdata vector for plotting
xdata = opts.TR:opts.TR:opts.TR*size(rBOLD,4);

% IMPORTANT use the index values to isolate the corresponding CO2/O2 values
PetCO2 = CO2trace(idx(1):idx(2));
% PetO2 = O2trace(1,idx(1):idx(2));

%% 6) Removing contributions using a large vessel mask
% Supply necessary options
% Define the cutoff percentile (higher values are removed; default = 98);
opts.LVpercentile = LVpercentile;

% If the stats toolbox is not present, then a manual threshold must 
% be applied. This can vary depending on the data (i.e. trial and error).
[mWBmask] = remLV(rBOLD,WBmask,opts);


%% 7) Temporally de-noise data
%de-noise
denBOLD = denoiseData(rBOLD, WBmask, opts);

clear cleanData
%% 8) Smooth and normalize data
% 3d Gaussian smoothing
opts.spatialdim = spatialdim; %default = 2

sBOLD = smthData(denBOLD, mWBmask, opts);

% Normalize data to first 15 baseline images
nBOLD = normTimeseries(sBOLD,mWBmask,normWindow);
clear sBOLD denBOLD
%% 9) Calculate CVR and hemodynamic lag
% lagCVR has the option to add a nuissance regressor (for example a PetCO2
% trace when focusing on using PetO2 as the main explanatory regressor).
% This regressor is only used when opts.glm_model = 1 (default = 0). For
% now we won't include anything.
L = [];

% Setup some extra options - most are default but we will initialized them
% for the sake of this tutorial (see manual)

% We would like to generate CVR maps
opts.cvr_maps = 1; %default is 1

% Factor by which to temporally interpolate data. Better for picking up
% lags between TR. Higher value means longer processing time and more RAM
% used (so be careful)
opts.interp_factor = interp_factor; %default is 4

% The correlation threshold is an important parameter for generating the
% optimized regressor. For noisy data, a too low value will throw an error.
% Ideally this should be set as high as possible, but may need some trial
% and error.
opts.corrthresh = corrthresh; %default is 0.7

% Thresholds for refinening optimized regressors. If you make this range too large 
% it smears out your regressor and lag resolution is lost. When using a CO2
% probe, the initial 'bulk' alignment becomes important here as well. A bad
% alignment will mean this range is also maybe not appropriate or should be
% widened for best results. (ASSUMING TR ~1s!, careful). 
opts.lowerlagthresh = lagthresh(1); %default is -3
opts.upperlagthresh = lagthresh(2); %default is 3

% Lag thresholds (in units of TR) for lag map creation. Since we are looking at a healthy
% brain, this can be limited to between 20-60TRs. For impariment you can consider
% to raise the upper threshold to between 60-90TRs (ASSUMING TR ~1s!, careful). 
opts.lowlag = lag(1); %setup lower lag limit; negative for misalignment and noisy correlation
opts.highlag = lag(2); %setups upper lag limit; allow for long lags associated with pathology

% For comparison we can also run the lagged GLM analysis (beware this can
% take quite some time)
opts.glm_model = 1; %default is 0
opts.corr_model = 1; %default is 1
% If you already have a probe from a previous run (with the correct
% length), you can opts.load_probe = 1 and the existing probe will be
% loaded - probe optimization will then be skipped
opts.load_probe = 0; %default is 0
% Produces the optimized regressor (default = 1) using the RAPIDTIDE approach.
% When this is set to 0; a straight correlation with the input probe is
% done. 

% For creating the optimized BOLD regressor, we will only use GM voxels for
% maximum CO2 sensitivity. We will also remove large vessel contributions.
mGMmask = int16(GMmask).*int16(mWBmask);
%remove CSF from analysis
mCSFmask = CSFmask -1; mCSFmask = abs(mCSFmask);
mask = int16(mCSFmask).*int16(mWBmask);

%Perform hemodyamic analysis
% The lagCVR function saves all maps and also returns them in a struct for
% further analysis. It also returns the optimized probe when applicable.
[newprobe, maps] = lagCVR(mGMmask, mask, nBOLD, PetCO2, L, opts);

% The advantage of using the global opts struct is that the variables used
% for a particular processing run (including all defaults set within
% functions themselves) can be saved to compare between runs.

save([opts.resultsdir, 'processing_options.mat'], 'opts');
close all;

disp("Finished Gas_CVR analysis.");

end
