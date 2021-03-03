
% Robin Haak,  25 February 2021

%% METHODS:

% The extracellular voltage traces are first preprocessed using common average referencing (optional):
% substracting each channel's median to remove baseline offsets, then substracting the median across
% all channels at each timepoint to remove artifacts. The data are spike sorted using Kilosort2
% with standard parameters, except for a maximal number of clusters that is equal to 4 x number of channels.
% Kilosort2 is able to track spikes of a neuron when its location relative to the probe
% changes (i.e., during drift). Additionaly, Kilosort2 performs automated splits/merges similar to what
% a human would do, based on spike waveform similarity, the bimodality of the distribution of
% waveform features, and auto- and cross correlograms. Kilosort2's output is then automatically curated
% based on a unit's presence troughout the recording, refractory period violations, and amplitude cutoff.

%% - CONFIGURATION(1) - Most of the parameters in this section you'll have to change ONCE

clear
close all

HeleroUserPath = 'D:\Users\user6\Documents';
addpath(genpath('D:\Users\user6\Documents\V1GRIA3_SingleUnits\SpikeSorting')) % path to 'SpikeSorting' folder
addpath('D:\Users\user6\Documents\npy-matlab-master\npy-matlab') % for converting to Phy
rootH = 'E:\'; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = 'D:\Users\user6\Documents\V1GRIA3_SingleUnits\SpikeSorting\Kilosort-2.0\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)

%% - CONFIGURATION(2) - Do not change anything here

run(fullfile(pathToYourConfigFile, 'StandardConfig_MOVEME.m'))
    ops.fproc = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    ops.trange = [0 Inf];
    
rhdDir = uigetdir(HeleroUserPath, 'Select directory that contains .rhd files to sort (GOOD FILES ONLY!)');

ops.chanMap = uigetfile(pathToYourConfigFile, 'Select channelmap');
    chanMapName = [pathToYourConfigFile filesep ops.chanMap];
    load(chanMapName, 'chanMap');
    ops.NchanTOT = numel(chanMap);
    ops.Nfilt = ops.NchanTOT*4; % maximum number of clusters

prompt = {'Do you want to apply common average referencing? 0 = no, 1 = yes'};
dlgtitle = 'Do you want to apply common average referencing?';
dims = [1 50];
definput = {'1'};
applyCAR = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput))); %#ok<ST2NM>

%% - RHDtoBinary

fprintf('Combining .rhd files into a binary file...\n');
[ops.fs, datFileName, nameDate] = RHDtoBinary(rhdDir);% Combine all .rhd recordings files into one binary file. Additionally, extract digital board inputs.
units.name_date = nameDate;
rootZ = [rhdDir, '\ksOutput'];
fprintf('Done\n');   

%% - COMMON AVERAGE REFERENCING

if applyCAR == 1
    fprintf('Common average referencing...\n'); 
    applyCARtoDat(datFileName, ops.NchanTOT, rootZ);
    fprintf('Done\n');
else
    fprintf('No common average referencing applied\n')
end


%% - RUN KILOSORT2

fprintf('Running Kilosort2...\n');
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
% Modified this section so that it searches for the _CAR.dat file when
% applyCAR == 1

if applyCAR == 1
    fs          = [dir(fullfile(rootZ, '*_CAR.bin')) dir(fullfile(rootZ, '*_CAR.dat'))]; % We do not use the .bin format
    ops.fbinary = fullfile(rootZ, fs(1).name);
elseif applyCAR == 0 
    fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
    ops.fbinary = fullfile(rootZ, fs(1).name);
end

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort2/issues/29
%rez = remove_ks2_duplicate_spikes(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');

%% - AUTOMATED CURATION OF KILOSORT'S OUTPUT + SAVING SPIKESTIMES / WAVEFORMS FOR THESE CLUSTERS

% Curate clustering output based on (1) presence troughout the recording, 
% (2) refractory period violations, and (3) amplitude distribution cutoff.
% N.B., We do not use measures based on PC components, as they are
% sensitive to drift and more difficult to compare to between recordings 

qc.rpv = 0.5; % <
qc.amplitude = 0.2; % < 
qc.presence = 0.9; % > 

[units] = postKS(rootZ, qc, units); % Creates the 'units' struct
fprintf('Press space after inspection of the results to manually exclude clusters (optional) and/or save results\n');
pause
            
% manually exclude a cluster based on, e.g., waveform shape (optional)
prompt = {'Manually exclude cluster(s)'};
dlgtitle = '(cluster ID(s):';
dims = [1 50];
definput = {''};
badClu = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput))); %#ok<ST2NM
cluKeep = ~ismember(units.cids, badClu);
fprintf('Removed %d cluster(s) \n', numel(units.cids) - sum(cluKeep));

% save only the clusters you want to keep
units.cids = units.cids(cluKeep);
units.spike_times = units.spike_times(cluKeep);
units.spike_samples = units.spike_samples(cluKeep);
units.waveforms = units.waveforms(cluKeep, :);
units.trough_to_peak = units.trough_to_peak(cluKeep);
units.rp_violations = units.rp_violations(cluKeep);
units.presence_ratio = units.presence_ratio(cluKeep);
units.amplitudes = units.amplitudes(cluKeep);
units.amplitude_cutoff = units.amplitude_cutoff(cluKeep);
units.amplitude_pdf = units.amplitude_pdf(cluKeep,:);
units.amplitude_dis = units.amplitude_dis(cluKeep,:);

units.presence = units.presence(cluKeep, :);

fname = fullfile(rhdDir, 'units.mat');
save(fname, 'units', '-v7.3');
close all

fprintf('----FINISHED!----\n');
