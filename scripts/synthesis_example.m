clc

% Analysis-synthesis of speech using various voice source modeling methods.
% Check the synthesized speech samples in folder "syn".

% Define directories and data
bin_dir = '../bin/';
wav_dir = '../wav/';
syn_dir = '../syn/';
cfg_dir = '../cfg/';
wav_name = 'mv_0001.wav';
config_filename = 'config_default';
config_user_plib = 'config_user_plib';
config_user_pca = 'config_user_pca';
config_user_dnn = 'config_user_dnn';
config_user_2pp = 'config_user_2pp';
config_user_2pp_new = 'config_user_2pp_new';
config_user_noiaif = 'config_user_noiaif';

% Create folder "syn"
if ~exist(syn_dir,'dir')
    system(['mkdir ' syn_dir]);
end

% Analysis
str_analysis = [bin_dir 'Analysis ' wav_dir wav_name ' ' cfg_dir config_filename];
disp(str_analysis);
[~,~] = system(str_analysis);

% Synthesis using single pulse
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' cfg_dir config_filename];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) 'singlepulse.wav']);

% Synthesis using pulse library
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_plib];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) 'pulselib.wav']);

% Synthesis using PCA pulse
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_pca];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) 'pcapulse.wav']);

% Synthesis using DNN-based voice source modeling
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_dnn];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) 'dnnpulse.wav']);

% Synthesis using single pre-computed 2-pitch period pulse (2pp)
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_2pp];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) '2pp.wav']);

% Synthesis using single 2-pitch period pulse separated from the pulse
% library (2pp_new)
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_2pp_new];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) '2pp_new.wav']);

% Synthesis using simple LPC
str_synthesis = [bin_dir 'Synthesis ' wav_dir wav_name(1:end-4) ' ' ...
    cfg_dir config_filename ' ' cfg_dir config_user_noiaif];
disp(str_synthesis);
[~,~] = system(str_synthesis);
system(['mv ' wav_dir wav_name(1:end-3) 'syn.wav ' syn_dir...
    wav_name(1:end-3) 'noiaif.wav']);
disp('Done')


