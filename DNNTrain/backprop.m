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

% This program fine-tunes an autoencoder with backpropagation.
% Weights of the autoencoder are going to be saved in mnist_weights.mat
% and trainig and test reconstruction errors in mnist_error.mat
% You can also set maxepoch, default value is 200 as in our paper.  

% Modified by Tuomo Raitio for HMM-based speech synthesis applications
% Nov. 6 2013


%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity check
if length(numhid) ~= numlayer-1
    disp(['Error: Vector ''numhid'' must be of length ' num2str(numlayer-1)]);
    break;
end

% Do not initialize if previous training data from workspace is used
if continue_previous_training == 0

    % DNN architecture
    [numcases numdims numbatches] = size(batchdata);
    numdims_out = size(batchdataTar,2);
    N = numcases;
    
    % Create DNN weight structure
    W = cell(numlayer,1);
    
    % Initialize weights either with RBM pretrain values
    if rbm_pretrain == 1
        for i = 1:numlayer
            load(['rbm_layer' int2str(i)]);
            if i == 1
                W{i} = [vishid; hidrecbiases];
            elseif i == numlayer
                W{i} = [hidtop; toprecbiases];
            else
                W{i} = [hidpen; penrecbiases];
            end
        end
    else
        % Initialize weights with random values
        for i = 1:numlayer
            if i == 1
                W{i} = 0.1*randn(size(D,2)+1,numhid(i));
            elseif i == numlayer
                W{i} = 0.1*randn(numhid(numlayer-1)+1,size(DTar,2));
            else
                W{i} = 0.1*randn(numhid(i-1)+1,numhid(i));
            end
        end
    end
    
    % Initialize
    start_epoch = 1;
    test_err = [];
    train_err = [];
else
    % Load previous training data
    load TrainData
    load TestData
    load DNN_Errors
    load DNN_Weights
    start_epoch = length(test_err) + 1;
end

% Dimensions
L = cell(numlayer+1,1);
for i = 1:numlayer
    L{i} = size(W{i},1)-1;
end
L{numlayer+1} = size(DTar,2);
Dim = zeros(numlayer+1,1);
for i = 1:numlayer+1
	Dim(i) = L{i};
end

% Print DNN architecture details
disp('DNN architecture:')
disp('   Layers:')
for i = 1:numlayer+1
    if i == 1
        disp(['      Layer ' int2str(i) ': ' num2str(numdims) ' (input)'])
    elseif i == numlayer+1
        disp(['      Layer ' int2str(i) ': ' num2str(numdims_out) ' (output)'])
    else
        disp(['      Layer ' int2str(i) ': ' num2str(numhid(i-1)) ' (hidden)'])
    end
end
disp('   Weights:')
for i = 1:numlayer
    disp(['      w' int2str(i) ': ' num2str(size(W{i},1)) ' x ' num2str(size(W{i},2))])
end
disp(' ')



%%%%%%%%%%%%%%%%%%%%%%% END OF INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% Display message
disp('Performing backpropagation...');
disp('   Fine-tuning DNN by minimizing MSE');

% Start training
wprobs = cell(numlayer-1,1);
for epoch = start_epoch:maxepoch

    %%%%%%%%%%%%%%%%%%%% COMPUTE TRAINING RECONSTRUCTION ERROR %%%%%%%%%%%%
    err = 0; 
    [numcases numdims numbatches] = size(batchdata);
    N = numcases;
    for batch = 1:numbatches

        % Set data
        data = [batchdata(:,:,batch) ones(N,1)];
        dataY = [batchdataTar(:,:,batch) ones(N,1)];

        % Layer 1
        wprobs{1} = 1./(1 + exp(-data*W{1}));
        wprobs{1} = [wprobs{1} ones(N,1)];
        
        % Layer 2...numlayer-1
        for i = 2:numlayer-1
            wprobs{i} = 1./(1 + exp(-wprobs{i-1}*W{i}));
            wprobs{i} = [wprobs{i} ones(N,1)];
        end
        
        % Top layer
        dataout = wprobs{numlayer-1}*W{numlayer};

        % Evaluate error
        err = err + 1/N*sum(sum((dataY(:,1:end-1)-dataout).^2)); 
    end
    train_err(epoch) = err/numbatches;

    %%%%%%%%%%%%%% END OF COMPUTING TRAINING RECONSTRUCTION ERROR %%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%%%% COMPUTE TEST RECONSTRUCTION ERROR %%%%%%%%%%%%%%%%
   
    [testnumcases testnumdims testnumbatches] = size(testbatchdata);
    N = testnumcases;
    err = 0;
    wprobs = cell(numlayer-1,1);
    for batch = 1:testnumbatches
        
        % Set data
        data = [testbatchdata(:,:,batch) ones(N,1)];
        dataY = [testbatchdataTar(:,:,batch) ones(N,1)];

        % Layer 1
        wprobs{1} = 1./(1 + exp(-data*W{1}));
        wprobs{1} = [wprobs{1} ones(N,1)];
        
        % Layer 2...numlayer-1
        for i = 2:numlayer-1
            wprobs{i} = 1./(1 + exp(-wprobs{i-1}*W{i}));
            wprobs{i} = [wprobs{i} ones(N,1)];
        end
        
        % Top layer
        dataout = wprobs{numlayer-1}*W{numlayer};

        % Evaluate error
        err = err +  1/N*sum(sum((dataY(:,1:end-1)-dataout).^2)); 
    end
    test_err(epoch) = err/testnumbatches;
    
    % Display errors
    disp(' ')
    disp(['   Before epoch ' num2str(epoch) ' Train squared error: ' num2str(train_err(epoch)) ...
        ' Test squared error: ' num2str(test_err(epoch))]);
    disp(' ')

    
    %%%%%%%%%%%%%% END OF COMPUTING TEST RECONSTRUCTION ERROR %%%%%%%%%%%%%
    
    % This defines the number of batches that are combined into a larger
    % batch. NCB = 1 means that no combining takes place
    % (10 was used earlier, which may be slower)
    NCB = Ncbatch;
    tt = 0;
    for batch = 1:numbatches/NCB 
        disp(['      epoch ' num2str(epoch) ' batch ' num2str(batch) '/' num2str(floor(numbatches/NCB))]);

        %%%%%%%%%%% COMBINE NCB MINIBATCHES INTO 1 LARGER BATCH %%%%%%%%%%%
        tt = tt + 1; 
        data = zeros(NCB*batchsize,numdims);
        dataY = zeros(NCB*batchsize,numdims_out);
        for kk = 1:NCB
            data((kk-1)*batchsize+1:kk*batchsize,:) = batchdata(:,:,(tt-1)*NCB+kk); 
            dataY((kk-1)*batchsize+1:kk*batchsize,:) = batchdataTar(:,:,(tt-1)*NCB+kk); 
        end

        %%%%%%%%%%%%%%% PERFORM CONJUGATE GRADIENT WITH M LINESEARCHES %%%%

        % Reshape for minimization
        VV = [];
        for i = 1:numlayer
            VV = [VV W{i}(:)'];
        end
        VV = VV';

        % Minimize
        [X, fX] = minimize(VV,'CG_MNIST',max_iter,Dim,data,dataY,learn_rate);

        % Reshape to form weight matrices
        W{1} = reshape(X(1:(L{1}+1)*L{2}),L{1}+1,L{2});
        xxx = (L{1}+1)*L{2};
        for i = 2:numlayer
            W{i} = reshape(X(xxx+1:xxx+(L{i}+1)*L{i+1}),L{i}+1,L{i+1});
            xxx = xxx + (L{i}+1)*L{i+1};
        end

        %%%%%%%%%%%%%%% END OF CONJUGATE GRADIENT WITH 3 LINESEARCHES %%%%%
    end

    % Save results
    save DNN_Weights.mat W data_min data_max
    save DNN_Errors.mat test_err train_err

end



