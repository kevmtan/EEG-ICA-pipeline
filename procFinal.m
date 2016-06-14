function procFinal(s, p)
% s = subject number; p = parameters struct
    % This function is called by procFinal_wrapper.m, which contains its parameters etc.
    % This performs final preprocessing using outputs from procPREP.m, procICAdipfit.m
    % Kevin Tan - Jun 8, 2016

% Load logFile
logFile = fopen([int2str(s) '_' p.task '_log.txt'], 'w+');
fprintf(logFile, ['START: ' datestr(now) '\n \n']);

% % Load Psychtoolbox paradigm data
% load(p.PTB);
% if strcmp(p.task, 'jig')
%     PTB = jigResults; %#ok<*NODEF>
% elseif strcmp(p.task, 'mem')
%     PTB = memResults;
% end

% Load structs from ICA stream
amica = load(p.AMICA);
pd = load(p.p_decomp);
pd = pd.preprocstruct;
EEG.subject = s;
EEG.etc.preproc_decomp = pd;
EEG.etc.preproc = p;
try
    dipfit = load(p.dipfit);
    dipfit = dipfit.dipfit;
catch
end


%% Start of EEGlab

% Load EEGlab
tic;
eeglab;
pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1,...
    'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
    'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
    'option_checkversion', 1, 'option_chat', 0);

% Load PREP'd datset
EEG = pop_loadset('filename', p.prepFile);
disp(['Loading ' p.prepFile]);
fprintf(logFile, ['%.2f - Loading ' p.prepFile '\n'], toc);

% Get interpolated channels from PREP
interpCh = EEG.etc.noiseDetection.reference.interpolatedChannels.all;
[~, interpChNames] = eeg_decodechan(EEG.chanlocs, interpCh);
EEG.etc.PREPinterpCh = interpCh; EEG.etc.PREPinterpChNames = interpChNames;
disp('PREP interpolated channels:');
disp(interpChNames);
disp(interpCh);
fprintf(logFile, '%.2f - PREP interpolated channels: ', toc);
fprintf(logFile, '%s\t', interpChNames{:});
fprintf(logFile, '\n');

% Remove large 'etc.noiseDetection' field from PREP
EEG.etc = rmfield(EEG.etc, 'noiseDetection');
disp('Removed large EEG.etc.noiseDetection field from PREP pipeline');
fprintf(logFile, '%.2f - removed large EEG.etc.noiseDetection field from PREP pipeline \n', toc);

% Remove mastoids to reduce rank and align with ICA processing
EEG = pop_select(EEG, 'nochannel', {'M1', 'M2'});
chanlocs = EEG.chanlocs;
disp('Removed mastoids to reduce rank');
fprintf(logFile, '%.2f - Removed mastoids to reduce rank \n', toc);

% Correct channel indices
[~, exChNames] = eeg_decodechan(EEG.chanlocs, p.exCh(1):size(EEG.data, 1));
[ECGch, ~] = eeg_decodechan(EEG.chanlocs, 'ECG'); %#ok<ASGLU>
EEG.etc.preproc.exCh_ChRm = exChNames;

% High pass filter FIR1
EEG = pop_eegfiltnew(EEG, p.HP.freq, 0);
EEG.setname = [EEG.setname '_' p.HP.name];
fprintf(logFile, ['%.2f - Hi-passed ' p.HP.name ' FIR1 zero-order \n'], toc);

% Epoch no baseline correction
EEG = pop_epoch(EEG, p.epEvents, [p.epMin p.epMax], 'newname', [EEG.setname '_ep']);
disp(['Epoched from ' num2str(p.epMin) ' to ' num2str(p.epMax) ', no baseline correction']);
fprintf(logFile, ['%.2f - Epoched from ' num2str(p.epMin) ' to ' num2str(p.epMax)...
    ', no baseline correction \n'], toc);

% % Insert psychtoolbox trial data in epoch struct
% if s == 2 && strcmp(p.task, 'jig')
%     epOffset = 1;
% else
%     epOffset = 0;
% end
% 
% epCheck = zeros(2, size(EEG.epoch, 2));
% for e = 1:size(EEG.epoch, 2)
%     for fn = fieldnames(PTB)' % don't get rid of transpose!
%         EEG.epoch(1, e).(fn{1}) = PTB(1, e+epOffset).(fn{1});
%     end
%     % Save bdf trigger and ptb stimtype for consistency check
%     if iscell(EEG.epoch(1, e).eventlatency) % bdf event type
%         eLock = find(cell2mat(EEG.epoch(1, e).eventlatency) == 0);
%         epCheck(1, e) = EEG.epoch(1, e).eventtype{1, eLock}; %#ok<*FNDSB>
%     else
%         epCheck(1, e) = EEG.epoch(1, e).eventtype;
%     end
%     epCheck(2, e) = EEG.epoch(1, e).stimID;
%     if isfield(PTB, 'jiggle') % stim ID
%         if EEG.epoch(1, e).jiggle == 1;
%             epCheck(2, e) = EEG.epoch(1, e).stimID + 100; % add 100 for jiggles
%         end
%     end
%     clear eLock
% end
% epConsist = find(diff(epCheck, 1, 1)); % Consistency check
% EEG.etc.preproc.epCheck = epCheck; % Save to struct vv
% EEG.etc.preproc.epConsist = epConsist;
% 
% if isempty(epConsist)
%     disp('Imported Psychtoolbox trial data, BDF/PTB data consistent');
%     fprintf(logFile, '%.2f - Imported Psychtoolbox trial data, BDF/PTB data consistent \n', toc);
% else
%     disp(['Imported Psychtoolbox trial data, found ' num2str(length(epConsist)) ' inconsistent trials:']);
%     disp(num2str(epConsist));
%     fprintf(logFile, '%.2f - Imported Psychtoolbox trial data, Found %d inconsistent trials: \n', toc, length(epConsist));
%     fprintf(logFile, ' %d', epConsist);
%     fprintf(logFile, '\n');
% end
% pop_saveset(EEG, 'filename', [EEG.setname]); %save
% fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);

% Remove PREP-interpolated chans, to be interpolated from ICA'd and cleaned data
EEG = pop_select(EEG, 'nochannel', interpChNames);
disp('Removed PREP-interpolated channels');
fprintf(logFile, '%.2f - Removed PREP-interpolated channels \n', toc);

% Temporarily remove 1 scalp channel if no PCA
if isfield(pd, 'scalpRm')
    scalpRm = pd.scalpRm;
    EEG = pop_select(EEG, 'nochannel', scalpRm); % Remove chan
    EEG.etc.preproc.scalpRm = scalpRm; % Save to EEG struct
    disp('Removed 1 scalp channel for ICA');
    disp(scalpRm);
    fprintf(logFile, '%.2f - Removed 1 scalp channel for ICA: ', toc);
    fprintf(logFile, '%s\n', scalpRm{:});
end

% Determine current # of chans
currentCh = size(EEG.data, 1);
disp(['Current # of channels: ' num2str(currentCh)]);
fprintf(logFile, ['%.2f - Current # of channels: ' num2str(currentCh) '\n'], toc);

% Generate ICA activations from ica matrices of ICA dataset
EEG.icaweights = amica.weights;
EEG.icasphere = amica.sphere;
EEG.etc.amica = amica.mods;
EEG = eeg_checkset(EEG, 'ica');
EEG.setname = [EEG.setname '_ICA']; %rename
disp(['Generated ' num2str(size(EEG.icaweights, 1)) ' IC activations using ICA matricies from ICA dataset']);
fprintf(logFile, ['%.2f - Generated ' num2str(size(EEG.icaweights, 1)) ' IC activations using ICA matricies from ICA dataset \n'], toc);
pop_saveset(EEG, 'filename', [EEG.setname]); %save
fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);

% Reject bad ICs
if isfield(pd, 'badICs') || (p.ICrej.hiRV.on == 1)
    % Reject MARA/FASTER bad ICs
    if isfield(pd, 'badICs')
        if ~isempty(pd.badICs)
            EEG = pop_subcomp(EEG, pd.badICs, 0);
            disp('Rejected bad ICs');
            fprintf(logFile, '%.2f - Rejected bad ICs \n', toc);
        end
    else
        disp('Rejected no bad ICs');
        fprintf(logFile, '%.2f - Rejected no bad ICs \n', toc);
    end
    
    % Reject high DIPFIT residual variance ICs
    if p.ICrej.hiRV.on == 1
        hiRV = find([dipfit.model.rv] >= p.ICrej.hiRV.rv);
        if ~isempty(hiRV)
            % Prepare for rejection
            EEG.etc.preproc.ICrej.hiRV.ICs = hiRV;
            remICs = setdiff(1:size(EEG.icaweights, 1), hiRV);
            
            % Rejection
            EEG = pop_subcomp(EEG, hiRV, 0);
            disp('Rejected high RV ICs:');
            disp(num2str(hiRV));
            fprintf(logFile, '%.2f - Rejected %d high-RV ICs: ', toc, length(hiRV));
            fprintf(logFile, ' %d', hiRV);
            fprintf(logFile, '\n');
        end
    else
        disp('Rejected no high-RV ICs');
        fprintf(logFile, '%.2f - Rejected no high-RV ICs \n', toc);
    end
    EEG.setname = [EEG.setname '_ICrej'];
end

% ICs to use for epoch rejection
goodICs = find([dipfit.model(remICs).rv] <= p.ICrej.hiRV.rvRej);
EEG.etc.preproc.goodICs = goodICs;

% Interpolate removed channels
EEG = pop_interp(EEG, chanlocs, 'spherical');
disp('Re-interpolated removed channels');
fprintf(logFile, '%.2f - Re-interpolated removed channels \n', toc);

% Baseline correction
EEG = pop_rmbase(EEG, [p.bMin 0]);
EEG.setname = [EEG.setname '_BL'];  %rename
disp(['Baseline corrected from: ' num2str(p.bMin) ' to 0']);
fprintf(logFile, ['%.2f - baseline corrected from: ' num2str(p.bMin) ' to 0 \n'], toc);
pop_saveset(EEG, 'filename', [EEG.setname]); %save
fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);

% Epoch rejection based on good ICA components
if p.epRej.on == 1
    
    % Find epochs w/ extreme voltages
    EEG = pop_eegthresh(EEG, 0, goodICs, p.epRej.uvMin, p.epRej.uvMax,...
        EEG.xmin, EEG.xmax, 0, 0);
    thresh = find(EEG.reject.icarejthresh);
    EEG.etc.preproc.epRej.thresh = thresh;
    fprintf(logFile, '%.2f - Found %d epochs with extreme voltages:',...
        toc, length(thresh));
    fprintf(logFile, ' %d', thresh);
    fprintf(logFile, '\n');
    
    % Find improbable epochs
    [EEG, ~, ~, ~, ~] = pop_jointprob(EEG, 0, goodICs,...
        p.epRej.probLoc, p.epRej.probGlb, 0, 0, 0);
    prob = find(EEG.reject.icarejjp);
    EEG.etc.preproc.epRej.prob = prob;
    fprintf(logFile, '%.2f - Found %d improbable epochs:',...
        toc, length(prob));
    fprintf(logFile, ' %d', prob);
    fprintf(logFile, '\n');
    
    % Find epochs with high temporal kurtosis
    EEG = pop_rejkurt(EEG, 0, goodICs, p.epRej.kurtLoc, p.epRej.kurtGlb, 0, 0, 0);
    kurt = find(EEG.reject.icarejkurt);
    EEG.etc.preproc.epRej.kurt = kurt;
    fprintf(logFile, '%.2f - Found %d epochs with high kurtosis:',...
        toc, length(kurt));
    fprintf(logFile, ' %d', kurt);
    fprintf(logFile, '\n');
    
    % Find trials via muscle frequency range deviation
    [EEG, ~] = pop_rejspec(EEG, 0, 'elecrange', goodICs,...
        'freqlimits', p.epRej.fLim, 'threshold', p.epRej.fThresh,...
        'method', 'fft', 'eegplotreject', 0);
    freqrej = find(EEG.reject.icarejfreq);
    EEG.etc.preproc.epRej.freq = freqrej;
    fprintf(logFile, '%.2f - Found %d epochs with high muscle noise:',...
        toc, length(freqrej));
    fprintf(logFile, ' %d', freqrej);
    fprintf(logFile, '\n');
    
    % Concactenate bad epochs
    EEG.etc.preproc.epRej.all = horzcat(EEG.etc.preproc.epRej.thresh,...
        EEG.etc.preproc.epRej.prob, EEG.etc.preproc.epRej.kurt,...
        EEG.etc.preproc.epRej.freq);
    EEG.etc.preproc.epRej.all = unique(EEG.etc.preproc.epRej.all, 'sorted'); % Remove duplicates
    
    % Reject bad epochs
    EEG = pop_select(EEG, 'notrial', EEG.etc.preproc.epRej.all);
    EEG.setname = [EEG.setname '_epRej'];  %rename
    fprintf(logFile, '%.2f - Removed %d bad epochs \n', toc, length(EEG.etc.preproc.epRej.all));
    pop_saveset(EEG, 'filename', [EEG.setname]); %save
    fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);
end

% Import IC dipole locations from ICA/DIPFIT dataset
if pd.dipfit == 1
    EEG.dipfit = dipfit;
    EEG.dipfit.model = dipfit.model(remICs);
    EEG = eeg_checkset(EEG, 'ica');
    disp('Imported IC dipole locations from ICA/DIPFIT dataset');
    fprintf(logFile, '%.2f - Imported IC dipole locations from ICA/DIPFIT dataset \n', toc);
    EEG.setname = [EEG.setname '_DIPFIT']; %rename
    pop_saveset(EEG, 'filename', [EEG.setname]); %save
    fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);
end

%% Separate into conditions
% Conditions need to be separated into their own .set files for EEGLAB's STUDY
% structure for group/statistical analysis. To do this you need to import 
% Psychtoolbox data into the EEG.epochs structure


% % Make conds dir
% mkdir('conds');
% cd('conds');
% disp('Making condition datasets:');
% fprintf(logFile, '%.2f - Making condition datasets: \n', toc);
%
% load('/data2/tarrlab/kevinalt/Elissa_NEIL/indoorStim.mat');
% outdoorStim = find(~ismember(1:100, indoorStim));
%
% % Separate the conditions for eeglab STUDY structure
% for c = 1:2
%
%     if c == 1
%         cName = 'indoor';
%         cStim = indoorStim;
%     elseif c == 2
%         cName = 'outdoor';
%         cStim = outdoorStim;
%     end
%
%     condEps = zeros(size(EEG.epoch, 2));
%     for e = 1:size(EEG.epoch, 2)
%         if any(EEG.epoch(1, e).stimID == cStim)
%             condEps(e) = 1;
%         end
%     end
%     tempEEG = pop_select(EEG, 'trial', find(condEps));
%     tempEEG.condition = cName;
%     tempEEG.setname = [EEG.setname '_' cName];
%     pop_saveset(tempEEG, 'filename', [tempEEG.setname]);
%     fprintf(logFile, ['Saved dataset: ' tempEEG.setname '\n']);
%     clear tempEEG e eLock condEps
% end


fprintf(logFile, '\n \n');
fprintf(logFile, ['END: ' datestr(now)]);

end

