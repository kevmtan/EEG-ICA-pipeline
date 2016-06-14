function procICAdipfit_wrapper(s)
% s = subject number
    % This function is a wrapper that calls procICAdipfit.m
    % Path, directories, and parameters for procICAdipfit.m are contained here
    % procICAdipfit.m performs ICA and DIPFIT localization on procPREP.m outputs
    % Kevin Tan - Jun 8, 2016

p = struct;
%% Path & Directories

% Custom matlab path
cd('/home/kevinalt/MATLAB');
path(pathdef);

% Root folder where the proc functions and subject files are
p.rootFolder = '/data2/tarrlab/kevinalt/Elissa_NEIL';
addpath([p.rootFolder]);

% CD to folder where PREP file is
cd([rootFolder '/' int2str(s) '/PREP']); % Change to subject directory

% Get PREP filename
prepname = 'oneBack_PREP'; % Searchable Filename of PREP file
p.prepFile = strtrim(ls(['*' prepname '.bdf']));

% Folder where preproc files will be saved (an ICA subfolder will be made here)
folder = '20160608';
mkdir([rootFolder '/' int2str(s) '/' folder]);
cd([rootFolder '/' int2str(s) '/' folder]);

%% Parameters

% Channel params
p.scalpCh = 1:128; % all scalp channels
p.exCh = 129:136; % all external channels
p.EOGch = 131:135; % EOG channels

% Epoch params
p.epMin = -.5; % Epoch start for ICA dataset (long as possible w/o overlap)
p.epMax = 1.25; % Epoch end for final & ICA dataset
p.epEvents = num2cell(1:200); % Event types (triggers) to epoch to

%% Initial ICA params (for finding epochs with one-of-a-kind artifacts)

% AMICA params
p.ICAi.on = 1; % (0 = off, 1 = on)
p.ICAi.kVal = 30; % K-value in determining number of ICA components (default = 30)
p.ICAi.threads = 36; % Threads to use for ICA (Psych-O processors have 6-8 cores)
p.ICAi.iter = 2000; % Max iterations
p.ICAi.rej = 0; % Reject data that doesn't fit AMICA's model
p.ICAi.rejSD = 99; % Z-score threshold for data rejection (outliers from AMICA's model)
p.ICAi.rejNum = 99; % Number of times data is rejected
p.ICAi.rejStart = 99; % First iteration at which to reject data
p.ICAi.rejInt = 99; % Intervals from 1st reject iteration at which to reject data

% Used to identify stereotyped artifactual ICs to not use for epoch rejection
    % .measure = enable (1) or disable (0)
    % .z = Z-score cutoff
p.CRi.on = 1; % enable (1) or disable (0)
p.CRi.o.measure(1) = 0; % Median Gradient (disabled, finds non-sterotyped ICs)
p.CRi.o.z(1) = 3;
p.CRi.o.measure(2) = 0; % Slope around LPF band (disabled, finds non-sterotyped ICs)
p.CRi.o.z(2) = 3;
p.CRi.o.measure(3) = 0; % Spatial Kurtosis (disabled, finds non-sterotyped ICs)
p.CRi.o.z(3) = 3;
p.CRi.o.measure(4) = 0; % Hurst Exponent (disabled, finds non-sterotyped ICs)
p.CRi.o.z(4) = 3;
p.CRi.o.measure(5) = 1; % EOG/ECG correlation (identifies stereotyped EOG/ECG ICs)
p.CRi.o.z(5) = 3;

%% Epoch rejection params (on initial ICs for clean final ICA decomp)
p.epRej.on = 1; % Reject epochs for final ICA decomp, yes(1) or no(0)
p.epRej.uvMin = -1000; % Low microvolt threshold
p.epRej.uvMax = 1000; % High microvolt threshold
p.epRej.probLoc = 8.75; % Per-component probability threshold
p.epRej.probGlb = 999; % All-component probability threshold
p.epRej.kurtLoc = 7.25; % Per-component kurtosis z-threshold
p.epRej.kurtGlb = 999; % All-component kurtosis z-threshold
% p.epRej.fLim = [20 40]; % Muscle frequency range
% p.epRej.fThresh = [-100 25]; % Muscle noise deviation threshold

%% Final ICA params (on clean epochs)

% AMICA params
p.ICAf.on = 1; % (0 = off, 1 = on)
p.ICAf.kVal = 30; % K-value in determining number of ICA components (default = 30)
p.ICAf.threads = 36; % Threads to use for ICA (Psych-O processors have 6-8 cores)
p.ICAf.iter = 2000; % Max iterations
p.ICAf.rej = 1; % Reject data that doesn't fit AMICA's model
p.ICAf.rejSD = 5; % Z-score threshold for data rejection (outliers from AMICA's model)
p.ICAf.rejNum = 5; % Number of times data is rejected
p.ICAf.rejStart = 4; % First iteration at which to reject data
p.ICAf.rejInt = 3; % Intervals from 1st reject iteration at which to reject data

% MARA params (artifactual IC classifier)
p.MARA = 1; % enable (1) or disable (0)
p.MARAprob = 0.25; % Posterior probability threshold (Kevin = 0.25) 

% ICA component rejection options (if used after MARA, be conservative!)
    % .measure = enable (1) or disable (0)
    % .z = Z-score cutoff
p.CRf.on = 1; % enable (1) or disable (0)
p.CRf.o.measure(1) = 1; % Median Gradient
p.CRf.o.z(1) = 6;
p.CRf.o.measure(2) = 1; % Slope around LPF band (not relevant since LPF wasn't performed)
p.CRf.o.z(2) = 6;
p.CRf.o.measure(3) = 1; % Spatial Kurtosis
p.CRf.o.z(3) = 3;
p.CRf.o.measure(4) = 1; % Hurst Exponent
p.CRf.o.z(4) = 6;
p.CRf.o.measure(5) = 1; % EOG/ECG correlation
p.CRf.o.z(5) = 6;

%% DIPFIT params
p.dipfit = 1; % (0 = off, 1 = on)

% Location of models/position info. Default files are located in the DIPFIT plugin dir in your EEGLAB dir
p.hdmfile = '/home/kevinalt/MATLAB/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/standard_BESA.mat'; % Location of desired head model
p.mrifile = '/home/kevinalt/MATLAB/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/avg152t1.mat'; % Location of anatomical MRI data (standardized MNI is used by default)
p.chanfile = '/home/kevinalt/MATLAB/eeglab13_5_4b/functions/resources/Standard-10-5-Cap385_noEx.elp'; % Location of channel positions (standardized 10/20 locations by default)

% Manually check for dual-dipole ICs after autofit? (1 = yes, 0 = no)
p.dipfitManual = 1; % (use importDIPFIT after manually checking)

%% Prepare for & call main ICA/dipfit preprocessing function

% Call ICA/dipfit preprocessing function
procICAdipfit(s, p);
end
