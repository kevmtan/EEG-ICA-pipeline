function procPREP(s)
% s = subject number
    % Runs the PREP pipeline for a BDF file (Bigdely-Shamlo et al., 2015)
    % Edit parameters below so this function knows how to find BDF file
    % Kevin Tan - Jun 8, 2016

%% Path, Directories & Parameters

% Load your custom matlab path
cd('/home/kevinalt/MATLAB');
path(pathdef);

% Change to folder where BDF file is
rootFolder = '/data2/tarrlab/kevinalt/Elissa_NEIL';
cd([rootFolder '/' int2str(s)]); % Change to subject directory

% Naming
taskbdfname = 'oneBack'; % Searchable BDF filename
task = 'oneBack'; % Base output filename

% Get BDF filename from folder
bdf = strtrim(ls(['*' taskbdfname '.bdf']));

% Make PREP directory in subject directory
mkdir('PREP');
cd('PREP');

%% PREP parameters (Kevin's preferences for CMU Psy Biosemi)
% Parameters explained in PrepPipeline/utilities/getPipelineDefaults.m
params = struct();
params.name = [num2str(s) '_' task '_PREP']; % Output filename

% Ignore boundary events
params.ignoreBoundaryEvents = true; % false if there are breaks in your recording

% Channel parameters for deterending and line noise removal
params.detrendChannels = 1:136; % scalp+external channels to detrend for processing
params.lineNoiseChannels = 1:136; % scalp+external channels for CleanLine

% Referencing parameters
params.referenceType = 'robust'; % Robust referencing (cleaned scalp electrode average)
params.meanEstimateType = 'huber'; % Use huber to find means
params.interpolationOrder = 'post-reference'; % Re-interpolate bad channels after referencing
params.referenceChannels = 1:128; % scalp channels to use for initial average reference
params.evaluationChannels = 1:128; % scalp channels that will contribute to the robust reference
params.rereferencedChannels = 1:136; % scalp+external channels that will be referenced to the robust reference

% Bad channel detection parameters (made more conservative since too many ch removed)
params.badTimeThreshold = 0.1;
params.highFrequencyNoiseThreshold = 8;
params.robustDeviationThreshold = 13;
params.ransacOff = true; % turned off, way too aggressive


%% Start of PREP

% Load EEGLab
eeglab

% EEGLab Options
pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 0, 'option_saveversion6', 0,...
    'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
    'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
    'option_checkversion', 1, 'option_chat', 0);

% Load BDF into EEGlab
disp(['Loading ' bdf]);
EEG = pop_biosig(bdf, 'channels', 1:136);
EEG.setname = params.name; % rename

% Save OG chan labels
EEG.etc.urchanlocs = EEG.chanlocs;

% Load channel locations
EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'lookup', '/home/kevinalt/MATLAB/eeglab13_5_4b/functions/resources/Standard-10-5-Cap385_noEx.elp');

% Trim start of dataset (anything prior to 1sec before first event)
firstEventPt = EEG.event(1, 1).latency;
try
    EEG = pop_select(EEG, 'nopoint', [1, firstEventPt - 513]);
    disp('trimmed beggining of dataset');
catch
end

% Trim end of dataset (anything past 1 sec after last event)
lastEventPt = EEG.event(1, size(EEG.event, 2)).latency;
try
    EEG = pop_select(EEG, 'nopoint', [lastEventPt + 1, EEG.pnts]);
    disp('trimmed end of dataset');
catch
end

% PREP pipeline
disp('Running PREP:');
[EEG, ~] = prepPipeline(EEG, params);

% Save PREP-interpolated channel names
[~, interpChNames] = eeg_decodechan(EEG.chanlocs,...
    EEG.etc.noiseDetection.reference.interpolatedChannels.all);
interpCh = EEG.etc.noiseDetection.reference.interpolatedChannels;
EEG.etc.PREPinterpChNames = interpChNames;
EEG.etc.PREPinterpCh = interpCh;
disp('interpolated channels:');
disp(interpChNames);
disp(interpCh);
save([num2str(s) '_' task '_interpChNames.mat'], 'interpChNames', '-mat', '-v7.3');
save([num2str(s) '_' task '_interpCh.mat'], 'interpCh', '-mat', '-v7.3');
disp('saved interp channel structs:');

% Save channel locations file
chanlocs = EEG.chanlocs;
save([num2str(s) '_' task '_chanlocs.mat'], 'chanlocs', '-mat', '-v7.3');
disp('saved channel location struct');

% Save PREP'd dataset
fname = [num2str(s) '_' task '_PREP.set'];
save(fname, 'EEG', '-mat', '-v7.3');
disp(['saved:' fname]);
end
