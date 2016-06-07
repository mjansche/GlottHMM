clc
clear
close all

% Settings
dnn_name = 'dnnw';
description = 'DNN trained from mv_0001 file';
inputpar = 'F0, Gain, HNR, LSFsource, LSF';
dim_in = [1 1 5 10 30];
outputpar = 'Glottal flow pulse derivative';
dim_out = 400;
author = 'Author Name';

% Load DNN trained data
disp('Saving DNN weights...')
load DNN_Weights
n = size(W,1);

% Load DNN errors
load DNN_Errors
epochs = length(test_err);

% Make directory
[~,curdir] = system('pwd');
dnn_dir = [curdir(1:end-1) '/' dnn_name];
if ~exist([curdir(1:end-1) '/' dnn_name],'dir')
    system(['mkdir ' dnn_dir]);
end

% Save weights to file
for i = 1:n
    w = W{i};
    w = reshape(w',size(w,1)*size(w,2),1);
    fid = fopen([dnn_dir '/dnn.w' num2str(i)],'wt');
    fprintf(fid,'%1.9f\n',w);
    fclose(fid);
end

% Save inout data min and max
load TestData.mat
fid = fopen([dnn_dir '/input.minmax'],'wt');
fprintf(fid,'%1.9f\n',[data_min data_max]);
fclose(fid);

% Save info
fid = fopen([dnn_dir '/readme'],'wt');
fprintf(fid,'%s\n',description);
fprintf(fid,'%s\n',author);
fprintf(fid,'%s\n',[date ' ' datestr(now, 'HH:MM:SS')]);
fprintf(fid,'\n');
fprintf(fid,'Input parameters : %s\n',inputpar);
fprintf(fid,'Input dimensions : %s\n',int2str(dim_in));
fprintf(fid,'Output parameters : %s\n',outputpar);
fprintf(fid,'Output dimensions : %s\n',int2str(dim_out));
fprintf(fid,'\n');
fprintf(fid,'Training epochs: %d\n', epochs);
fprintf(fid,'DNN size:\n');
for i = 1:n
    s1 = size(W{i},1);
    s2 = size(W{i},2);
    fprintf(fid,'  Layer %s: ',int2str(i));
    fprintf(fid,'%d x %d\n',s1,s2);
end
fclose(fid);
disp('Done')





  

  


