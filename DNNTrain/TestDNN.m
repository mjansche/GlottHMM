%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TestDNN
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load train and test data
load TestData.mat
load TestDataTar.mat
load('DNN_Weights.mat');
load('DNN_Errors.mat');

% Plot errors
start_index = 2;
figure(1)
plot(train_err(start_index:end))
hold on
plot(test_err(start_index:end),'r')
legend('Train error','Test error');
xlabel('Training epochs')
ylabel('Error')
disp('Press any key to continue...');
pause
close all


%%%%%%%%%%%%%%%%%%%% COMPUTE PARMETER(S) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
numlayer = length(W);
wprobs = cell(numlayer-1);

% Set bias term to data
data = D;
N = size(data,1);
data = [data ones(N,1)];

% Layer 1
wprobs{1} = 1./(1 + exp(-data*W{1}));
wprobs{1} = [wprobs{1} ones(N,1)];

% Layer 2...numlayer-1
for i = 2:numlayer-1
    wprobs{i} = 1./(1 + exp(-wprobs{i-1}*W{i}));
    wprobs{i} = [wprobs{i} ones(N,1)];
end

% Top layer
Ddnn = wprobs{numlayer-1}*W{numlayer};


%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot results
figure(1)
for i = 1:N
    clf
    plot(DTar(i,:))
    hold on
    plot(Ddnn(i,:),'k')
    xlabel('Samples')
    ylabel('Amplitude')
    legend('Original','DNN','Location','SouthEast')
    pause
end
disp('Done')




