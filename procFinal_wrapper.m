function procFinal_wrapper(s)
% s = subject number
    % This function is a wrapper that calls procFinal.m
    % Path, directories, and parameters for procFinal.m are contained here
    % procFinal.m requires outputs from procPREP.m, procICAdipfit.m
    % Kevin Tan - Jun 8, 2016

%% Parameters

% Load custom pathdef (Make a custom pathdef! Psych-o default paths messed up)
cd('/home/kevinalt/MATLAB');
path(pathdef);

% Root folder where the proc functions and subject files are
p.rootFolder = '/data2/tarrlab/kevinalt/Elissa_NEIL';
addpath([p.rootFolder]);

% Change to subject dir
cd([p.rootFolder '/' int2str(s)]);

% Name of dir where ICA subdir is in
folder = '20160608';

% Task name
p.task = 'oneBack';

% Filenames for structs to import
p.prepFile = [p.rootFolder '/' int2str(s) '/PREP/' int2str(s) '_' p.task '_PREP.set']; % PREP output fule
p.AMICA = [p.rootFolder '/' int2str(s) '/' folder '/ICA/' int2str(s) '_ICAf.mat']; % final ICA matrix
p.dipfit = [p.rootFolder '/' int2str(s) '/' folder '/ICA/' int2str(s) '_dipfit.mat']; % DIPFIT model
p.p_decomp = [p.rootFolder '/' int2str(s) '/' folder '/ICA/' int2str(s) '_' p.task '_PREP_1hz_ep_ICAi_epRej_ICAf_ICrej_DIPFIT_preproc.mat']; % Preproc struct output/parameters from ICA dataset
% p.PTB = strtrim(ls([p.rootFolder '/' int2str(s) '/S' int2str(s) '_' p.task 'NeilEEG_*_struct.mat'])); % Struct output from psychtoolbox paradigm

% Channel params
p.scalpCh = 1:128; % all scalp channels
p.exCh = 129:136; % all external channels
p.EOGch = 131:135; % EOG channels

% Epoch params
p.epMin = -.5; % Epoch start (seconds before event/trigger)
p.epMax = 1.25; % Epoch end (seconds after event/trigger)
p.epEvents = num2cell(1:200); % Event types (triggers) to epoch to
p.bMin = -200; % baseline correction start time (ms)

%% Filtering params

p.HP.freq = 0.1; % Hi-pass frequency
p.HP.name = '0-1hz'; % String to add to setname/filname

%% Epoch rejection params

p.epRej.on = 1; % Reject epochs for ICA, yes(1) or no(0)
p.epRej.uvMin = -1000; % Low microvolt threshold
p.epRej.uvMax = 1000; % High microvolt threshold
p.epRej.probLoc = 6.75; % Per-component probability threshold
p.epRej.probGlb = 3.5; % All-component probability threshold
p.epRej.kurtLoc = 6.75; % Per-component kurtosis z-threshold
p.epRej.kurtGlb = 3.5; % All-component kurtosis z-threshold
p.epRej.fLim = [20 40]; % Muscle frequency range
p.epRej.fThresh = [-100 25]; % Muscle noise deviation threshold


%% Additional IC rejection based on DIPFIT residual variance (RV)
p.ICrej.hiRV.on = 1; % 1=on, 0=off

p.ICrej.hiRV.rv = 0.20; % RV > threshold for IC rejection
p.ICrej.hiRV.rvRej = 0.13; % RV < threshold for ICs to use for epoch rejection

%% Prepare for & call main ICA/dipfit preprocessing function

% Call ICA/dipfit preprocessing function
procFinal(s, p);
end

