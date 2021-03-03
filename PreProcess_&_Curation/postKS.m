
% Robin Haak, 1 March 2021

function [units] = postKS(ksDir, qc, units)

% Quality metrics:
% - REFRACTORY PERIOD VIOLATIONS (0-1) This metric is based on: Hill DN, Mehta SB, Kleinfeld D. Quality metrics 
%   to accompany spike sorting of extracellular signals. Journal of Neuroscience. 2011 Jun 15;31(24):8699-705.
% - PRESENCE RATIO (0-0.99) This metric comes from the Allen Institute Neuropixels Visual Coding
%   dataset. It measures the fraction of time during a session in which a unit is spiking, and ranges from 0 to 0.99
% - AMPLITUDE CUTOFF (0-0.5) This metric is inspired by Hill et al. 2011 and the Allen Institute Neuropixels Visual
%   Coding dataset.

%% First, extract spiketimes/amplitudes for each cluster

fprintf('Extracting spiketimes and waveform features for each cluster...\n')

clu = readNPY(fullfile(ksDir, 'spike_templates.npy')) + 1; % Because the data is not manually sorted; correct zero-indexing
cids = unique(clu);

spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));
spkSamples = readNPY(fullfile(ksDir, 'spike_times.npy'));
spkTimes = double(spkSamples)/spikeStruct.sample_rate;
amplitudes = readNPY(fullfile(ksDir, 'amplitudes.npy')); % Not in physical units
spkTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy'));

units.cids = cids;
units.spike_times = cell(numel(cids), 1);
units.spike_samples = cell(numel(cids), 1);
units.amplitudes = cell(numel(cids), 1);

for i = 1:numel(cids)
    spks = clu==cids(i);
    units.spike_samples{i} = spkSamples(ismember(spks, 1));
    units.spike_times{i} = spkTimes(ismember(spks, 1));
    units.amplitudes{i} = amplitudes(ismember(spks, 1));
end

temps = readNPY(fullfile(ksDir, 'templates.npy'));
winv = readNPY(fullfile(ksDir,  'whitening_mat_inv.npy'));

tempsUnW = zeros(size(temps)); % Unwhiten all the templates
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

[~,max_site] = max(max(abs(temps),[],2),[],3); % Get the channel with the maximum amplitude and extract waveform from this channel
templates_max = nan(size(temps,1),size(temps,2));
for curr_template = 1:size(temps,1)
    templates_max(curr_template,:) = ...
        temps(curr_template,:,max_site(curr_template));
end

[~,waveform_trough] = min(templates_max,[],2);
[~, templateDuration] = arrayfun(@(x) ...
    max(templates_max(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(templates_max,1)));

units.waveforms = [];
units.trough_to_peak = [];

for i = 1:numel(cids)
    units.waveforms(i, :) = templates_max(i, :);
    units.trough_to_peak(i, 1) = (double(templateDuration(i))/spikeStruct.sample_rate)*1e6; % trough-to-peak latency in microseconds9
end

fprintf('Done\n')

%% Cluster quality control

fprintf('Filtering out clusters that miss a lot of spikes and/or are highly contaminated... \n')

% Refractory period violations
units.rp_violations = [];
refDur = 0.0015; % Refractory period
minISI = 0.0005; % Censored/'shadow' period

for i = 1:numel(cids)
    spks = units.spike_times{i};
    totalRate = length(spks)/spks(end);
    nViolations = sum(diff(spks) <=refDur);
    violationTime = 2*length(spks)*(refDur-minISI);
    
    violationRate = nViolations/violationTime;
    units.rp_violations(i, 1) = violationRate/totalRate;
    
    if units.rp_violations(i, 1) > 1
       units.rp_violations(i, 1) = 1;
    end
end

% Presence ratio
units.presence = [];
units.presence_ratio = [];

for i = 1:numel(cids)
    spks = units.spike_times{i};
    min_time = min(spks);
    max_time = max(spks);
    a = histcounts(spks, (linspace(min_time, max_time)));
    units.presence(i,:) = a;
    units.presence_ratio(i, 1) = sum(a > 0)/100;
end

% Amplitude cutoff
units.amplitude_pdf = [];
units.amplitude_cutoff = [];

for i = 1:numel(cids)
    amplClu = units.amplitudes{i};
    
    [b, edges] = histcounts(amplClu, 500, 'Normalization', 'pdf');
    pdf = smoothdata(b, 'Gaussian');
            units.amplitude_pdf (i,:) = pdf;
            units.amplitude_dis(i,:) = b;
    binsize = mean(diff(edges));
    [~, peak_idx] = max(pdf);
    
    cutOff = min(abs(pdf(peak_idx:end) - pdf(1)));
    g = cutOff + peak_idx;
    G = fix(g + peak_idx);
    fractionMiss = sum(pdf(G:end))*binsize;

    if fractionMiss <= 0.5
       units.amplitude_cutoff(i, 1) = fractionMiss;
    else
       units.amplitude_cutoff(i, 1) = 0.5;
    end
end

% Determine which clusters pass
QCpass = [];

for i = 1:numel(cids)
    if units.rp_violations(i) < qc.rpv && units.presence_ratio(i) > qc.presence && units.amplitude_cutoff(i) < qc.amplitude
        QCpass(i, 1) = 1;
    else
        QCpass(i, 1) = 0;
    end
end

fprintf('Done\n');
fprintf('%d clusters passed quality control \n', sum(QCpass));

%% Plot the 'good' clusters

if sum(QCpass) > 0  

    cluKeep = find(QCpass);
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

    for i = 1:numel(units.cids)
        figure;
        cluTitle = ['Cluster ' num2str(units.cids(i))];
        sgtitle(cluTitle);
               
            subplot(2,2,1)
            yline(0, '--', 'Color', '#696969', 'Linewidth', 1); hold on
            waveform = units.waveforms(i,:);
            z = ~waveform;
            waveform(z) = [];
            plot(waveform, 'k', 'Linewidth', 1.5);
            title('Average waveform');
            ylabel('Amplitude (ab.units)');
            xlabel('Samples'); hold off
            xlim([1, Inf]);


            subplot(2,2,2)
            presence = units.presence(i,:);
            plot(presence, 'k', 'Linewidth', 1.5);
            title('Unit presence');
            ylabel('Spike count');
            xlabel('Total recording time (binned)');
                        xlim([1, Inf]);


            subplot(2,2,3)
            st = units.spike_times{i};
            [K, ~, ~, ~, ~] = ccg(st, st, 50, 0.001);
            K(51) = 0; % remove central bin for plotting
            bar(K, 'FaceColor', '#696969','EdgeColor','none', 'BarWidth', 1);
            title('Autocorrelogram');
            ylabel('# of events');
            xlabel('Latency (ms)');
            xlim([1, 100]);
            xticks([1 51 100]);
            xticklabels({'-50', '0', '50'}); 

            subplot(2,2,4);
            amps = units.amplitude_dis(i,:);
            z2 = ~amps;
            amps(z2) = [];
            bar(amps, 'FaceColor', '#696969', 'EdgeColor','none', 'BarWidth', 1); hold on    
            pdf2 = units.amplitude_pdf(i, :);
            pdf2(z2) = [];
            plot(pdf2, 'k', 'Linewidth', 1.5 )
            title('Amplitude');
            ylabel('Probability density');
            xlabel('Amplitude (binned)');
            xlim([1, Inf]); hold off
            
    end
end

fprintf('Done\n')
      
