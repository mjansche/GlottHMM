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

function [f, df] = CG_MNIST(VV,Dim,XX,YY,learn_rate)

% Dimensions
N = size(XX,1);
numlayer = length(Dim)-1;
for i = 1:length(Dim)
    L{i} = Dim(i);
end

% Reshape to form weight matrices
W{1} = reshape(VV(1:(L{1}+1)*L{2}),L{1}+1,L{2});
xxx = (L{1}+1)*L{2};
for i = 2:numlayer
    W{i} = reshape(VV(xxx+1:xxx+(L{i}+1)*L{i+1}),L{i}+1,L{i+1});
    xxx = xxx + (L{i}+1)*L{i+1};
end

% Initialize
XX = [XX ones(N,1)];
wprobs = cell(numlayer);

% Layer 1
wprobs{1} = 1./(1 + exp(-XX*W{1}));
wprobs{1} = [wprobs{1} ones(N,1)];

% Layer 2...numlayer-1
for i = 2:numlayer-1
    wprobs{i} = 1./(1 + exp(-wprobs{i-1}*W{i}));
    wprobs{i} = [wprobs{i} ones(N,1)];
end

% Top layer
XXout = wprobs{numlayer-1}*W{numlayer};
YY = [YY ones(N,1)];


% BACK-PROPAGATION

% MAX ENTROPY
%f = -1/N*sum(sum( YY(:,1:end-1).*log(XXout) + (1-YY(:,1:end-1)).*log(1-XXout)));

% MSE
Dw = cell(numlayer,1);
Ix = cell(numlayer,1);

% Error
f = 1/N*sum((sum( (YY(:,1:end-1)-(XXout)).^2 )));
IO = learn_rate/N*(XXout-YY(:,1:end-1));

% Output layer
Ix{numlayer} = IO;
Dw{numlayer} = wprobs{numlayer-1}'*Ix{numlayer};

% Lower layers
for i = numlayer-1:-1:2
    Ix{i} = (Ix{i+1}*W{i+1}').*wprobs{i}.*(1-wprobs{i}); 
    Ix{i} = Ix{i}(:,1:end-1);
    Dw{i} = wprobs{i-1}'*Ix{i};
end

% Input layer
Ix{1} = (Ix{2}*W{2}').*wprobs{1}.*(1-wprobs{1}); 
Ix{1} = Ix{1}(:,1:end-1);
Dw{1} =  XX'*Ix{1};

% Reshape
df = [];
for i = 1:numlayer
    df = [df Dw{i}(:)'];
    
end
df = df';




