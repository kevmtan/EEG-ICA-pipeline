function procICAdipfit(s, p)
% s = subject number; p = parameters struct
    % This function is called by procICAdipfit_wrapper.m, which contains its parameters etc.
    % This performs ICA and DIPFIT localization on PREP output files
    % Kevin Tan - Jun 8, 2016

%%

% Load EEGlab
tic;
eeglab;
pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1,...
    'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
    'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
    'option_checkversion', 1, 'option_chat', 0);

% Make new folder
mkdir('ICA');
cd('ICA');

% Start log
logFile = fopen([int2str(s) '_decomp_log.txt'], 'w+');
fprintf(logFile, ['START: ' datestr(now) '\n \n \n']);

% Load PREP'd datset
disp(['Loading ' p.prepFile]);
EEG = pop_loadset('filename', p.prepFile);
fprintf(logFile, ['%.2f - Loaded ' p.prepFile '\n'], toc);

% Get interpolated channels from PREP
interpCh = EEG.etc.noiseDetection.reference.interpolatedChannels.all;
[~, interpChNames] = eeg_decodechan(EEG.chanlocs, interpCh);
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

% Save things to etc struct
EEG.subject = s;
EEG.etc.preproc = p;
EEG.etc.PREPinterpChNames = interpChNames;
EEG.etc.PREPinterpCh = interpCh;

% High pass filter 1hz FIR1
EEG = pop_eegfiltnew(EEG, 1, 0);
EEG.setname = [EEG.setname '_1hz'];
fprintf(logFile, '%.2f - Hi-passed 1hz FIR1 zero-order \n', toc);

% Epoch no baseline correction
EEG = pop_epoch(EEG, p.epEvents, [p.epMin p.epMax], 'newname', [EEG.setname '_ep']);
disp(['Epoched from ' num2str(p.epMin) ' to ' num2str(p.epMax)...
    ', no baseline correction']);
fprintf(logFile, ['%.2f - Epoched from ' num2str(p.epMin) ' to ' num2str(p.epMax)...
    ', no baseline correction \n \n \n'], toc);

% Remove mastoids to reduce rank
EEG = pop_select(EEG, 'nochannel', {'M1', 'M2'});
chanlocs = EEG.chanlocs;
fprintf(logFile, '%.2f - Removed mastoids to reduce rank \n', toc);
disp('Removed mastoids to reduce rank');

% Correct channel indices
[~, exChNames] = eeg_decodechan(EEG.chanlocs, p.exCh(1):size(EEG.data, 1));
EEG.etc.preproc.exCh_ChRm = exChNames;

% Remove PREP-interpolated chans to reduce non-linearity for ICA decomposition
EEG = pop_select(EEG, 'nochannel', interpChNames);

% Remove 1 scalp channel if maxIC > currentCh to try & avoid PCA
disp('Attempting to avoid PCA, temporarily removing 1 scalp channel to correct robust reference rank deficiency');
fprintf(logFile, '%.2f - Attempting to avoid PCA, temporarily removing 1 scalp channel to correct robust reference rank deficiency \n', toc);

% Find scalp channel to remove (this could be done better?)
try
    [~, scalpRm] = eeg_decodechan(EEG.chanlocs, 'Cz');
catch %#ok<*CTCH>
    try
        [~, scalpRm] = eeg_decodechan(EEG.chanlocs, 'FCz');
    catch
        try
            [~, scalpRm] = eeg_decodechan(EEG.chanlocs, 'CPz');
        catch
            try
                [~, scalpRm] = eeg_decodechan(EEG.chanlocs, 'C1');
            catch
                try
                    [~, scalpRm] = eeg_decodechan(EEG.chanlocs, 'C2');
                catch
                    warning('Unable to remove a scalp channel, PCA necessary'); %#ok<*WNTAG>
                    fprintf(logFile, '%.2f - WARNING: Unable to remove a scalp channel, PCA necessary \n', toc);
                end
            end
        end
    end
end

% Temporarily remove 1 scalp channel to correct robust reference rank deficiency
if exist('scalpRm', 'var')
    EEG = pop_select(EEG, 'nochannel', scalpRm); % Remove chan
    EEG.etc.preproc.scalpRm = scalpRm; % Save to EEG struct
    disp('Removed 1 scalp channel for ICA');
    disp(scalpRm);
    fprintf(logFile, '%.2f - Removed 1 scalp channel for ICA: ', toc);
    fprintf(logFile, '%s\n', scalpRm{:});
    
    % Determine current # of chans
    currentCh = size(EEG.data, 1);
    disp(['Current # of channels: ' num2str(currentCh)]);
    fprintf(logFile, ['%.2f - Current # of channels: ' num2str(currentCh) '\n'], toc);
end

% Correct external channel indices
[exCh, exChNames] = eeg_decodechan(EEG.chanlocs, exChNames);
[ECGch, ~] = eeg_decodechan(EEG.chanlocs, 'ECG');


%% Initial AMICA decomposition for artifactual epoch rejection
disp('Initial AMICA decomposition for artifactual epoch rejection:');

% Determine max # of ICs (depends on k-value, dataset length, # interpolated chans)
maxIC = floor(sqrt(size(EEG.data(:, :), 2) / p.ICAi.kVal));
currentCh = size(EEG.data, 1);
disp(['Current # of channels: ' num2str(currentCh)]);
fprintf(logFile, ['%.2f - Current # of channels: ' num2str(currentCh) '\n'], toc);
disp(['Data length supports up to ' num2str(maxIC) ' IC decompositions']);
fprintf(logFile, ['%.2f - Data length supports up to ' num2str(maxIC) ' IC decompositions \n'], toc);

% Set AMICA folder
amicaout = [pwd '/amica_init'];

% Find if PCA necessary
if maxIC < currentCh % PCA necessary
    % Further dimension reduction necessary due to data length
    nmIC = maxIC;
    disp(['PCA necessary, reducing data to ' num2str(nmIC) ' components then running ICA:']);
    fprintf(logFile, ['%.2f - PCA necessary, reducing data to ' num2str(nmIC) ' components then running ICA: \n'], toc);
else % Normal AMICA
    nmIC = currentCh;
    disp(['No PCA necessary, running ICA for ' num2str(nmIC) ' components']);
    fprintf(logFile, ['%.2f - No PCA necessary, running ICA for ' num2str(nmIC) ' components \n'], toc);
end

% Run AMICA
[weights, sphere, mods] = runamica12(EEG.data(:, :), 'outdir', amicaout,...
    'use_queue', 0, 'max_threads', p.ICAi.threads, 'numprocs', p.ICAi.threads,...
    'max_iter', p.ICAi.iter, 'pcakeep', nmIC, 'do_reject', p.ICAi.rej,...
    'rejsig', p.ICAi.rejSD, 'numrej', p.ICAi.rejNum, 'rejstart', p.ICAi.rejStart,...
    'rejint', p.ICAi.rejInt);

% Add AMICA matrices to EEG data & compute ICA activations
EEG.icaweights = weights;
EEG.icasphere = sphere;
EEG.etc.ICAi.weights = weights;
EEG.etc.ICAi.sphere = sphere;
EEG.etc.ICAi.amica = mods;
EEG = eeg_checkset(EEG, 'ica'); % Compute ICA activations
EEG.setname = [EEG.setname '_ICAi'];

% Save dataset/AMICA matrices for later use
save([num2str(s) '_ICAi.mat'], 'weights', 'sphere', 'mods');
disp(['Ran ICA for ' num2str(nmIC) ' components, matrices saved to ICAi.mat']);
fprintf(logFile, ['%.2f - Ran ICA for ' num2str(nmIC) ' components, matrices saved to ICAi.mat \n'], toc);

% Find EOG & ECG ICs to not include for epoch rejection
if p.CRi.on == 1
    list_properties = component_properties(EEG, exCh);
    [exICs] = min_z(list_properties, p.CRi.o);
    exICs = find(exICs);
    erICs = 1:size(EEG.icasphere, 1);
    erICs = setdiff(erICs, exICs);
    EEG.etc.preproc.epRej.exICs = exICs; % Save to EEG struct
    EEG.etc.preproc.epRej.erICs = erICs; % Save to EEG struct
    disp(['Found' num2str(length(exICs)) 'external ICs:']);
    disp(exICs);
    fprintf(logFile, '%.2f - Found %d EOG/ECG ICs:', toc, length(exICs));
    fprintf(logFile, ' %d', exICs);
    fprintf(logFile, '\n');
else
    erICs = 1:size(EEG.icasphere, 1);
    EEG.etc.preproc.epRej.erICs = erICs; % Save to EEG struct
end

if p.epRej.on == 1
    % Find epochs w/ extreme voltages
    EEG = pop_eegthresh(EEG, 0, erICs, p.epRej.uvMin, p.epRej.uvMax,...
        EEG.xmin, EEG.xmax, 0, 0);
    EEG.etc.preproc.epRej.thresh = find(EEG.reject.icarejthresh);
    fprintf(logFile, '%.2f - Found %d epochs with extreme voltages:',...
        toc, length(EEG.etc.preproc.epRej.thresh));
    fprintf(logFile, ' %d', EEG.etc.preproc.epRej.thresh);
    fprintf(logFile, '\n');
    
    % Find improbable epochs
    [EEG, ~, ~, ~, ~] = pop_jointprob(EEG, 0, erICs,...
        p.epRej.probLoc, p.epRej.probGlb, 0, 0, 0);
    EEG.etc.preproc.epRej.prob = find(EEG.reject.icarejjp);
    fprintf(logFile, '%.2f - Found %d improbable epochs:',...
        toc, length(EEG.etc.preproc.epRej.prob));
    fprintf(logFile, ' %d', EEG.etc.preproc.epRej.prob);
    fprintf(logFile, '\n');
    
    % Find epochs with high temporal kurtosis
    EEG = pop_rejkurt(EEG, 0, erICs, p.epRej.kurtLoc, p.epRej.kurtGlb, 0, 0, 0);
    EEG.etc.preproc.epRej.kurt = find(EEG.reject.icarejkurt);
    fprintf(logFile, '%.2f - Found %d epochs with high kurtosis:',...
        toc, length(EEG.etc.preproc.epRej.kurt));
    fprintf(logFile, ' %d', EEG.etc.preproc.epRej.kurt);
    fprintf(logFile, '\n');
    
    % Concactenate bad epochs
    EEG.etc.preproc.epRej.all = horzcat(EEG.etc.preproc.epRej.thresh,...
        EEG.etc.preproc.epRej.prob, EEG.etc.preproc.epRej.kurt);
    EEG.etc.preproc.epRej.all = unique(EEG.etc.preproc.epRej.all, 'sorted'); % Remove duplicates
    
    % Save dataset
    pop_saveset(EEG, 'filename', [EEG.setname]); % Save
    fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);
    
    % Reject bad epochs
    EEG = pop_rejepoch(EEG, EEG.etc.preproc.epRej.all, 0);
    EEG.setname = [EEG.setname '_epRej'];  %rename
    fprintf(logFile, '%.2f - Removed %d bad epochs \n', toc, length(EEG.etc.preproc.epRej.all));
else
    % Save dataset
    pop_saveset(EEG, 'filename', [EEG.setname]); % Save
    fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);
end

%% Final ICA decomposition & DIPFIT on clean epochs
disp('Final AMICA decomposition:');
fprintf(logFile, '\n \n');
fprintf(logFile, '%.2f - Final AMICA decomposition: \n \n', toc);

% Remove old ICA stuff
clear weights sphere mods
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icawinv = [];
EEG.icachansind = [];
EEG.chaninfo.icachansind = [];
EEG.icaact = [];

% Determine max # of ICs (depends on k-value, dataset length, # interpolated chans)
maxIC = floor(sqrt(size(EEG.data(:, :), 2) / p.ICAf.kVal));
currentCh = size(EEG.data, 1);
disp(['Current # of channels: ' num2str(currentCh)]);
fprintf(logFile, ['%.2f - Current # of channels: ' num2str(currentCh) '\n'], toc);
disp(['Data length supports up to ' num2str(maxIC) ' IC decompositions']);
fprintf(logFile, ['%.2f - Data length supports up to ' num2str(maxIC) ' IC decompositions \n'], toc);

% Set AMICA folder
amicaout = [pwd '/amica'];

% Find if PCA necessary
if maxIC < currentCh % PCA necessary
    % Further dimension reduction necessary due to data length
    nmIC = maxIC;
    disp(['PCA necessary, reducing data to ' num2str(nmIC) ' components then running ICA:']);
    fprintf(logFile, ['%.2f - PCA necessary, reducing data to ' num2str(nmIC) ' components then running ICA: \n'], toc);
else % Normal AMICA
    nmIC = currentCh;
    disp(['No PCA necessary, running ICA for ' num2str(nmIC) ' components']);
    fprintf(logFile, ['%.2f - No PCA necessary, running ICA for ' num2str(nmIC) ' components \n'], toc);
end

% Run AMICA
[weights, sphere, mods] = runamica12(EEG.data(:, :), 'outdir', amicaout,...
    'use_queue', 0, 'max_threads', p.ICAf.threads, 'numprocs', p.ICAf.threads,...
    'max_iter', p.ICAf.iter, 'pcakeep', nmIC, 'do_reject', p.ICAf.rej,...
    'rejsig', p.ICAf.rejSD, 'numrej', p.ICAf.rejNum, 'rejstart', p.ICAf.rejStart,...
    'rejint', p.ICAf.rejInt);

% Add AMICA matrices to EEG data & compute ICA activations
EEG.icaweights = weights;
EEG.icasphere = sphere;
EEG.etc.amica = mods;
EEG = eeg_checkset(EEG, 'ica'); % Compute ICA activations
EEG.setname = [EEG.setname '_ICAf'];

% Save dataset/AMICA matrices for later use
save([num2str(s) '_ICAf.mat'], 'weights', 'sphere', 'mods');
disp(['Ran ICA for ' num2str(nmIC) ' components, matrices saved to ICAf.mat']);
fprintf(logFile, ['%.2f - Ran ICA for ' num2str(nmIC) ' components, matrices saved to ICAf.mat \n'], toc);

% Classify artifactual ICs using MARA
if p.MARA == 1
    [EEG.etc.preproc.badICs_MARAclassifier, MARAinfo] = MARA(EEG);
    badICs_MARA = find(MARAinfo.posterior_artefactprob > p.MARAprob);
    EEG.etc.preproc.badICs_MARA = badICs_MARA; % Save to EEG struct
    EEG.etc.preproc.MARAinfo = MARAinfo; % Save to EEG struct
    disp('MARA-classified artifactual ICs:');
    disp(num2str(badICs_MARA));
    fprintf(logFile, '%.2f - MARA classified %d artifactual ICs:', toc, length(badICs_MARA));
    fprintf(logFile, ' %d', badICs_MARA);
    fprintf(logFile, '\n');
end

% Find remaining bad ICs using FASTER z-thresholds
if p.CRf.on == 1
    remICs = setdiff(1:size(EEG.icaweights, 1), badICs_MARA); % find remaining ICs
    ICstats = findICstats(EEG, ECGch, remICs); % Calculate IC stats
    [badICs_FASTER] = min_z(ICstats, p.CRf.o); % Z-threshold
    badICs_FASTER = find(badICs_FASTER)';
    EEG.etc.preproc.badICs_FASTER = badICs_FASTER; % Save to EEG struct
    disp('FASTER bad ICs:');
    disp(num2str(badICs_FASTER));
    fprintf(logFile, '%.2f - FASTER found %d bad ICs:', toc, length(badICs_FASTER));
    fprintf(logFile, ' %d', badICs_FASTER);
    fprintf(logFile, '\n');
end

% Concactenate bad ICs
badICs = [];
if p.MARA == 1 && p.CRf.on == 1
    badICs = horzcat(badICs_MARA, badICs_FASTER);
    badICs = unique(badICs, 'sorted'); % Remove duplicates (shouldn't happen)
elseif p.MARA == 1
    badICs = badICs_MARA;
elseif p.FASTER.ica_options.component_rejection == 1
    badICs = badICs_FASTER;
end

% Save dataset
pop_saveset(EEG, 'filename', [EEG.setname]); % Save
fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);

% Reject bad ICs
if ~isempty(badICs)
    EEG = pop_subcomp(EEG, badICs, 0);
    EEG.etc.preproc.badICs = badICs;
    disp('Rejected bad ICs:');
    disp(num2str(badICs));
    fprintf(logFile, '%.2f - Rejected %d ICs \n', toc, length(badICs));
else
    EEG.etc.preproc.badICs = [];
    disp('Rejected no ICs');
    fprintf(logFile, '%.2f - Rejected no ICs \n', toc);
end
EEG.setname = [EEG.setname '_ICrej'];

% Re-interpolate removed channels
EEG = pop_interp(EEG, chanlocs, 'spherical');
disp('Re-interpolated removed channels');
fprintf(logFile, '%.2f - Re-interpolated removed channels \n', toc);

% Save dataset
pop_saveset(EEG, 'filename', [EEG.setname]); % Save
fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);

% DIPFIT
if p.dipfit == 1
    disp('Finding component dipole locations with DIPFIT:');
    
    % Create temporary dataset for DIPFIT without external channels
    EEG = pop_select(EEG, 'nochannel', exChNames);
    disp('Removed external channels \n');
    fprintf(logFile, '%.2f - Removed external channels \n', toc);
    
    % Find dipole locations
    EEG = pop_dipfit_settings(EEG, 'hdmfile', p.hdmfile, 'coordformat', 'Spherical',...
        'mrifile', p.mrifile, 'chanfile', p.chanfile);
    EEG = pop_multifit(EEG, 1:size(EEG.icaweights, 1), 'dipoles', 1, 'rmout', 'on');
    
    % Transfer dipfit files to EEG struct
    dipfit = EEG.dipfit;
    
    % Save dataset
    EEG.setname = [EEG.setname '_DIPFIT'];
    pop_saveset(EEG, 'filename', [EEG.setname]); % Save
    fprintf(logFile, ['Saved dataset: ' EEG.setname '\n \n']);
    
    % Save dipfit as separate mat file
    if p.dipfitManual == 0
        save([num2str(s) '_dipfit.mat']', 'dipfit');
        fprintf(logFile, '%.2f - Found component dipole locations with DIPFIT, saved to dipfit.mat \n', toc);
    else
        fprintf(logFile, '%.2f - Found component dipole locations with DIPFIT, awaiting manual check for dual-dipole ICs \n', toc);
    end
end

% Save preproc stuff to EEG struct
preprocstruct = EEG.etc.preproc;
save([EEG.setname '_preproc.mat'], 'preprocstruct');
end

%% Find independent component stats
function ICstats = findICstats(EEG, exCh, remICs)

ICstats = zeros(size(EEG.icaact,1),5); %This 5 corresponds to number of measurements made.

for u = remICs
    measure = 1;
    % TEMPORAL PROPERTIES
    
    % 1 Median gradient value, for high frequency stuff
    ICstats(u,measure) = median(diff(EEG.icaact(u,:)));
    measure = measure + 1;
    
    % 2 Mean slope around the LPF band (spectral) [DEPRECIATED]
    ignore_lpf=1;
    ICstats(u,measure) = 0;
    measure = measure + 1;
    
    % SPATIAL PROPERTIES
    
    % 3 Kurtosis of spatial map (if v peaky, i.e. one or two points high
    % and everywhere else low, then it's probably noise on a single
    % channel)
    ICstats(u,measure) = kurt(EEG.icawinv(:,u));
    measure = measure + 1;
    
    % OTHER PROPERTIES
    
    % 4 Hurst exponent
    ICstats(u,measure) = hurst_exponent(EEG.icaact(u,:));
    measure = measure + 1;
    
    % 10 Eyeblink correlations
    if (exist('exCh','var') && ~isempty(exCh))
        for v = 1:length(exCh)
            if ~(max(EEG.data(exCh(v),:))==0 && min(EEG.data(exCh(v),:))==0);
                f = corrcoef(EEG.icaact(u,:),EEG.data(exCh(v),:));
                x(v) = abs(f(1,2)); %#ok<*AGROW>
            else
                x(v) = v;
            end
        end
        ICstats(u,measure) = max(x);
        measure = measure + 1;
    end
end

for u = 1:size(ICstats,2)
    ICstats(isnan(ICstats(:,u)),u)=nanmean(ICstats(:,u));
    ICstats(:,u) = ICstats(:,u) - median(ICstats(:,u));
end
end
