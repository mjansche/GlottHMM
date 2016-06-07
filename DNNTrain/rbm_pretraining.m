% RBM pre-training of each layer
% Tuomo Raitio
% Nov. 6 2013

% Display message
disp('Performing RBM pretraining...');

% Sanity check
if length(numhid) ~= numlayer-1
    disp(['Error: Vector ''numhid'' must be of length ' num2str(numlayer-1)]);
    break;
end

% Train layer 1
disp(['   Pretraining Layer 1 with RBM: ' num2str(numdims) ' x ' num2str(numhid(1))]);
restart = 1;
grbm; % Gaussian RBM for real-valued input
hidrecbiases = hidbiases; 
save rbm_layer1 vishid hidrecbiases visbiases;

% Train hidden layers
for i = 1:numlayer-2
    disp(['   Pretraining Layer ' int2str(i+1) ' with RBM: ' num2str(numhid(i)) ' x ' num2str(numhid(i+1))]);
    batchdata_tmp = batchposhidprobs;
    restart = 1;
    rbm;
    hidpen = vishid;
    penrecbiases = hidbiases;
    hidgenbiases = visbiases;
    save(['rbm_layer' num2str(i+1) '.mat'],'hidpen','penrecbiases','hidgenbiases');
end

% Train top layer 
disp(['   Pretraining Layer ' int2str(numlayer) ' with RBM: ' num2str(numhid(end)) ' x ' num2str(numdims_out)]);
batchdata_tmp = batchposhidprobs;
numhid_out = numdims_out;
restart = 1;
rbmhidlinear;
hidtop = vishid;
toprecbiases = hidbiases;
topgenbiases = visbiases;
save(['rbm_layer' num2str(numlayer) '.mat'],'hidtop','toprecbiases','topgenbiases');
disp(' ')

