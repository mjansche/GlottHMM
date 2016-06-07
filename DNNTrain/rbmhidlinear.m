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

% This program trains Restricted Boltzmann Machine in which
% visible, binary, stochastic pixels are connected to
% hidden, tochastic real-valued feature detectors drawn from a unit
% variance Gaussian whose mean is determined by the input from 
% the logistic visible units. Learning is done with 1-step Contrastive Divergence.
% The program assumes that the following variables are set externally:
% maxepoch_rbm  -- maximum number of epochs
% numhid_out    -- number of hidden units
% batchdata_tmp -- the data that is divided into batches (numcases numdims numbatches)
% restart       -- set to 1 if learning starts from beginning

% Modified by Tuomo Raitio for HMM-based speech synthesis applications
% Nov. 6 2013

epsilonw      = 0.001; % Learning rate for weights 
epsilonvb     = 0.001; % Learning rate for biases of visible units
epsilonhb     = 0.001; % Learning rate for biases of hidden units 
weightcost  = 0.0002;  
initialmomentum  = 0.5;
finalmomentum    = 0.9;

[numcases numdims numbatches] = size(batchdata_tmp);

if restart == 1
    % Initializing symmetric weights and biases
    restart = 0;
    epoch = 1;
    vishid     = 0.1*randn(numdims, numhid_out);
    hidbiases  = zeros(1,numhid_out);
    visbiases  = zeros(1,numdims);
    poshidprobs = zeros(numcases,numhid_out);
    neghidprobs = zeros(numcases,numhid_out);
    posprods    = zeros(numdims,numhid_out);
    negprods    = zeros(numdims,numhid_out);
    vishidinc  = zeros(numdims,numhid_out);
    hidbiasinc = zeros(1,numhid_out);
    visbiasinc = zeros(1,numdims);
    sigmainc = zeros(1,numhid_out);
    batchposhidprobs = zeros(numcases,numhid_out,numbatches);
end

for epoch = epoch:maxepoch_rbm
    errsum = 0;
    for batch = 1:numbatches
        disp(['      Layer ' int2str(i+2) ': epoch ' int2str(epoch) ' batch ' int2str(batch)]);

        %%%%%%%%% START POSITIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data = batchdata_tmp(:,:,batch);
        poshidprobs =  (data*vishid) + repmat(hidbiases,numcases,1);
        batchposhidprobs(:,:,batch)=poshidprobs;
        posprods    = data' * poshidprobs;
        poshidact   = sum(poshidprobs);
        posvisact = sum(data);

        %%%%%%%%% END OF POSITIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        poshidstates = poshidprobs+randn(numcases,numhid_out);

        %%%%%%%%% START NEGATIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        negdata = 1./(1 + exp(-poshidstates*vishid' - repmat(visbiases,numcases,1)));
        neghidprobs = (negdata*vishid) + repmat(hidbiases,numcases,1);
        negprods  = negdata'*neghidprobs;
        neghidact = sum(neghidprobs);
        negvisact = sum(negdata); 

        %%%%%%%%% END OF NEGATIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        err = sum(sum((data-negdata).^2)); 
        errsum = err + errsum;
        if epoch > 5
        	momentum = finalmomentum;
        else
        	momentum = initialmomentum;
        end

        %%%%%%%%% UPDATE WEIGHTS AND BIASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vishidinc = momentum*vishidinc + ...
                    epsilonw*( (posprods-negprods)/numcases - weightcost*vishid);
        visbiasinc = momentum*visbiasinc + (epsilonvb/numcases)*(posvisact-negvisact);
        hidbiasinc = momentum*hidbiasinc + (epsilonhb/numcases)*(poshidact-neghidact);
        vishid = vishid + vishidinc;
        visbiases = visbiases + visbiasinc;
        hidbiases = hidbiases + hidbiasinc;

        %%%%%%%%%%%%%%%% END OF UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    disp(' ')
    disp(['   Before epoch ' num2str(epoch) ' error: ' num2str(errsum)])
    disp(' ')
end
