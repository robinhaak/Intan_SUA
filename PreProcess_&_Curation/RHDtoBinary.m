
function [sample_rate, datFileName, nameDate] = RHDtoBinary(rhdDir)

% Input = rhdDir, the directory where the .rhd files (of one(!) recording
% session are stored). N.B., only include the 'good' recordings, as this
% function combines all the .rhd files into one binary file for sorting.

% By Robin Haak, 24 February 2021

rhdFiles = dir([rhdDir '\*.rhd']); % list all rhd files
[~, idx] = sort({rhdFiles.date});
rhdFiles = rhdFiles(idx);

for i = 1:length(rhdFiles) % print files to check order
    fprintf('%s\n',rhdFiles(i).name)
end

mkdir(rhdDir, 'ksOutput');
nameDate = rhdFiles(1).name(1:end-7);
datFileName = [rhdDir '\ksOutput\' nameDate '.dat'];

fid = fopen(datFileName, 'w'); % open .dat file for writing

for i = 1:length(rhdFiles)
    fprintf('Loading file %i of %i, %s\n',i, length(rhdFiles),fullfile(rhdDir,rhdFiles(i).name));
    [amplifier_data, board_dig_in_data, ~, sample_rate] = read_Intan_RHD2000_file2(fullfile(rhdDir, rhdFiles(i).name));
   % [amplifier_data, board_dig_in_data, board_adc_data, sample_rate] = read_Intan_RHD2000_file2(fullfile(rhdDir, rhdFiles(i).name));
    
    diFileName = [rhdDir filesep rhdFiles(i).name(1:end-4),'_DigitalInputs.mat'];
    save(diFileName,'board_dig_in_data') % save digital inputs
    
%     if exist('board_adc_data')
%         aiFileName = [rhdDir filesep rhdFiles(i).name(1:end-4),'_AnalogInputs.mat'];
%         save(aiFileName,'board_adc_data') % save analog inputs
%     end
    
    fwrite(fid, amplifier_data(:),'int16'); % append to .dat file
end

fclose('all');
fprintf('Done\n');