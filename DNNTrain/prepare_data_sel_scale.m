% Tuomo Raitio
% Nov. 6 2013

% Display message
disp('Preparing data...');

% Load data (format: data_in (N x dim_in), data_out (N x dim_out))
disp('   Loading raw data...')
load(matfile)
N = size(data_in,1);
dim_in = size(data_in,2);
dim_out = size(data_out,2);

% Discard data that is too far from the mean
if discard_outlier_data == 1
    rmse_limit = 1.0;
    sel = zeros(N,1);
    m = mean(data_out,1);            
    rmse = zeros(N,1);
    for i = 1:N
        rmse(i) = sum((data_out(i,:)-m).^2);
    end
    sel(rmse <= rmse_limit) = 1;
    data_in = data_in(sel == 1,:);
    data_out = data_out(sel == 1,:);
    N = size(data_in,1);
end

% Use only part of the data if p_train indicates so (random selection)
if p_data < 1.0
    N = round(p_data*N);
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
    sel = randperm(size(data_in,1));
    sel = sel(1:N);
    data_in = data_in(sel,:);
    data_out = data_out(sel,:);
end

% Scale input data to range [0.1 0.9]
if scale_input_data == 1
    data_min = min(data_in);
    data_max = max(data_in);
    data_in = 0.1 + 0.8*(data_in - repmat(data_min, size(data_in,1),1))./repmat(data_max - data_min, size(data_in,1),1);
end

% Divide data to train and test data (random selection)
n_train = round(p_train*N);
n_test = N - n_train;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
allind = randperm(N);
train_ind = allind(1:n_train);
test_ind = allind(n_train+1:end);

% Save train data
disp('   Saving train data...')
D = data_in(train_ind,:);
DTar = data_out(train_ind,:);
save TrainData.mat D data_min data_max -v7.3
save TrainDataTar.mat DTar -v7.3

% Save test data
disp('   Saving test data...')
D = data_in(test_ind,:);
DTar = data_out(test_ind,:);
save TestData.mat D data_min data_max -v7.3
save TestDataTar.mat DTar -v7.3

% Clear data
clear data_in data_out
disp(' ')



