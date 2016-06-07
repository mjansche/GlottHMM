clc

% Analysis-synthesis of a speech file

% Define directories and files
bin_dir = '../bin/';
cfg_dir = '../cfg/';
wav_dir = '../wav/';
syn_dir = '../syn/';
wav_file = 'mv_0001.wav';
cfg_file = 'config_default';

% Analysis
str_analysis = [bin_dir 'Analysis ' wav_dir wav_file ' ' cfg_dir cfg_file];
disp(str_analysis);
system(str_analysis);

% Synthesis
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_file(1:end-4) ' ' cfg_dir cfg_file];
disp(str_synthesis);
system(str_synthesis);

% Move synthesis file to folder "syn"
if ~exist(syn_dir,'dir')
    system(['mkdir ' syn_dir]);
end
system(['mv ' wav_dir wav_file(1:end-3) 'syn.wav ' syn_dir]);