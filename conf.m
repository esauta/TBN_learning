%%% CONFIGURATION SCRIPT %%%%

threshold= 10^(-3); % end criterion: threshold for BIC difference of consecutive iterations
%this threshold should be changed depending on the global score of the initial network structure  

learning= 1; %to create the starting workspace 

max_osc= 10; % #of trials to escape from a local minimum situation

max_iter= 0; % you can set this variable to 2 for calculating the threshold for the end criterion

n_run = 1:100; %#of algorithm iterations
% where to save outputs
run_date= datestr(datetime('today')) ;
fold_name= strcat('Run_',run_date);
% set the directory in which all codes and input file are stored
path_def= pwd;
