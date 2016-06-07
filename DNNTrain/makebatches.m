% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.

% Modified by Tuomo Raitio for HMM-based speech synthesis applications
% Nov. 6 2013


% Display message
disp('Making batches...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make batches for train data
load('TrainData.mat');
totnum = size(D,1);
numdims = size(D,2);
numbatches = floor(totnum/batchsize);
if numbatches == 0
    disp('Error: Not enough train data or too large batch size');
    break;
end
disp(['   ' num2str(numbatches) ' minibatches of ' num2str(batchsize) ' cases each']);
disp(['   In backpropagation, combine ' num2str(Ncbatch) ' minibatches into 1 larger batch'])
disp('   Minibatches:')
disp(['      Size of training dataset: ' num2str(totnum) ' x ' num2str(numdims)]);
rand('state',0); % So we know the permutation of the training data
randomorder = randperm(totnum);
batchdata = zeros(batchsize, numdims, numbatches);
for b = 1:numbatches
    batchdata(:,:,b) = D(randomorder(1+(b-1)*batchsize:b*batchsize), :);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make batches for train target data
load('TrainDataTar.mat');
totnum = size(DTar,1);
numdims = size(DTar,2);
numbatches = floor(totnum/batchsize);
disp(['      Size of target train dataset: ' num2str(totnum) ' x ' num2str(numdims)]);
rand('state',0); % So we know the permutation of the training data
randomorder = randperm(totnum);
batchdataTar = zeros(batchsize, numdims, numbatches);
for b = 1:numbatches
    batchdataTar(:,:,b) = DTar(randomorder(1+(b-1)*batchsize:b*batchsize), :);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make batches for test data
load('TestData.mat');
totnum = size(D,1);
numdims = size(D,2);
numbatches = floor(totnum/batchsize);
if numbatches == 0
    disp('Error: Not enough test data or too large batch size');
    break;
end
disp(['      Size of testing dataset: ' num2str(totnum) ' x ' num2str(numdims)]);
rand('state',0); % So we know the permutation of the training data
randomorder = randperm(totnum);
testbatchdata = zeros(batchsize, numdims, numbatches);
for b = 1:numbatches
    testbatchdata(:,:,b) = D(randomorder(1+(b-1)*batchsize:b*batchsize), :);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make batches for test target data
load('TestDataTar.mat');
totnum = size(DTar,1);
numdims  =  size(DTar,2);
numbatches = floor(totnum/batchsize);
disp(['      Size of target test dataset: ' num2str(totnum) ' x ' num2str(numdims)]);
rand('state',0); % So we know the permutation of the training data
randomorder = randperm(totnum);
testbatchdataTar = zeros(batchsize, numdims, numbatches);
for b = 1:numbatches
    testbatchdataTar(:,:,b) = DTar(randomorder(1+(b-1)*batchsize:b*batchsize), :);
end


% Dimensions
[numcases numdims numbatches] = size(batchdata);
numdims_out = size(batchdataTar,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset random seeds 
rand('state',sum(100*clock)); 
randn('state',sum(100*clock));
disp(' ')


