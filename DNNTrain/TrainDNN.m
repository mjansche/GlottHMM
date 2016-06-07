%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TrainDNN
%
% 19.1.2015
% (c) Tuomo Raitio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you use GlottHMM within an article, paper, report or any other work
% that you wish to publish, please cite it as:
% 
% T. Raitio, A. Suni, J. Yamagishi, H. Pulakka, J. Nurminen, M. Vainio,
% and P. Alku, "HMM-based speech synthesis utilizing glottal inverse
% filtering," IEEE Trans. on Audio, Speech and Lang. Proc., vol. 19, no. 1,
% pp. 153-165, Jan. 2011.
% 
% T. Raitio, A. Suni, H. Pulakka, M. Vainio, and P. Alku, "Utilizing
% glottal source pulse library for generating improved excitation signal
% for HMM-based speech synthesis," in Proc. ICASSP, 2011, pp. 4564-4567.
% 
% T. Raitio, A. Suni, L. Juvela, M. Vainio, and P. Alku, "Deep neural
% network based trainable voice source model for synthesis of speech with
% varying vocal effort," in Proc. of Interspeech, 2014, pp. 1969-1973.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% DNN settings
numlayer = 3;       % Number of layers (excl. input, 3 yields a deep net)
numhid = [200 200]; % Units in hidden layer (length: numlayer-1)

% Training settings
maxepoch = 200;     % Maximum number of training epochs
maxepoch_rbm = 50;  % Maximum number of training epochs in RBM pretraining
batchsize = 50;     % Number of data points in each batch
Ncbatch = 10;       % Number of combined batches in minimization
max_iter = 3;       % Number of line searches
learn_rate = 5;     % Learning rate

% Data
matfile = 'pulselib1.mat'; % For loading new data matrix
p_data = 1.0;              % Proportion of data used (default 1.0)
p_train = 0.9;             % Prop. of train data, rest is used for eval.
scale_input_data = 1;      % Scale data to range [0.1, 0.9]
discard_outlier_data = 0;  % Discard pulses whose RMSE from average > limit
rmse_limit = 1.0;          % RMSE limit for discarding pulses

% Train phases
prepare_data = 1;                  % Set 1 if loading new data matrix
rbm_pretrain = 0;                  % Set 1 for RBM pre-training
continue_previous_training = 0;    % Set 1 if continuing existing training

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display message
clc
disp('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
disp('<>             Training of Deep Neural Network                  <>')
disp('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
disp(' ')

% Load data
if prepare_data == 1 && continue_previous_training == 0
    prepare_data_sel_scale;
end

% Make batches
makebatches;

% Perform RBM pretraining
if rbm_pretrain == 1 && continue_previous_training == 0
    rbm_pretraining;
end

% Perform backpropagation
backprop;

