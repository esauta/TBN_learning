clear all
format long
tic
path_def= pwd;
% where to import input data
path_input=strcat(path_def,'\','Input\'); %change slashes depending on the environment used (linux or windows)

%% IMPORT ALL THE INPUT FILES
file_GenesSbinding = fopen(strcat(path_input,'TF_genes_binding.txt'),'r'); % weighted edges TF:Gene 
dataArray = textscan(file_GenesSbinding, '%s%s%f%[^\n\r]', 'Delimiter', '\t');
dataArray(3) = cellfun(@(x) num2cell(x), dataArray(3), 'UniformOutput', false);
bindingP_Genes = [dataArray{1:end-1}];

Gedges= bindingP_Genes(:,1:2);

Gnames = unique(Gedges(:,2));
TFnames = unique(Gedges(:,1));

file_whitelist = fopen(strcat(path_input,'whitelist.txt'),'r');
fwhiteL = textscan(file_whitelist, '%s%s%[^\n\r]', 'Delimiter', '\t');
whitelist = [fwhiteL{1:end-1}];


file_TFbinding = fopen(strcat(path_input,'edges_TFTF_initialDAG_binding.txt'),'r');
dataArray_tf = textscan(file_TFbinding, '%s%s%f%[^\n\r]', 'Delimiter', '\t');
dataArray_tf(3) = cellfun(@(x) num2cell(x), dataArray_tf(3), 'UniformOutput', false);
bindingP_tf = [dataArray_tf{1:end-1}];

TFedges = bindingP_tf(:,1:2);

file_TFedges_all = fopen(strcat(path_input,'All_TFTF_edges.txt'),'r'); 
fedgesTF_all = textscan(file_TFedges_all, '%s%s%f%[^\n\r]', 'Delimiter', '\t');
fedgesTF_all(3) = cellfun(@(x) num2cell(x), fedgesTF_all(3), 'UniformOutput', false);
all_arcs_list = [fedgesTF_all{1:end-1}];

file_ExpressG = importdata(strcat(path_input,'Expression_matrix.txt'));
expr_values= file_ExpressG.data;
Gene_Expres_list = file_ExpressG.textdata;
expr_values_original = expr_values;

clearvars file_GenesSbinding dataArray file_whitelist fwhiteL file_TFbinding dataArray_tf file_TFedges_all fedgesTF_all file_ExpressG;

disp 'import and input preparation: done'
fprintf(['number of source TFs:  %i \n'],length(TFnames))
fprintf(['number of target Genes:  %i \n'],length(Gnames))

TFedges_original= TFedges;
%% create the initial total matrixes and the relative network structure object
[TF_starting_matr, total_matrix, total_matrix_binding, TFnames_topSorted,indexes_top_Order] = create_TMatrix(TFedges,TFnames, bindingP_tf,bindingP_Genes,Gedges,Gnames);
disp 'matrix creation: done'

Numero_nodi= length(Gnames) + length(TFnames);
numb_cont_values =  ones(1,Numero_nodi);
nodes_names = ({TFnames_topSorted{:},Gnames{:}}).'; 
nodes_index= 1:Numero_nodi;

%create the structure object in which all model parameter will be stored
bnet = mk_bnet(total_matrix, numb_cont_values,'names',nodes_names,'observed', nodes_index,'discrete', []);
disp 'bnet: done'

%% Assignment of CPDs and learning parameters on the initial Model

dataE_complete = zeros(length(total_matrix),size(expr_values,2));            
for nodo =1:length(nodes_names)
    nameN = nodes_names(nodo);
    indice_gExpr = find(strcmp(nameN, Gene_Expres_list)==1);
    evaluesN= (expr_values(indice_gExpr,:));
    dataE_complete(nodo,:) = evaluesN;
end
 
bnet_wo_CPD= bnet; %save the pre CPD's assignment bnet
for n=1:length(nodes_names)
    bnet.CPD{n}= gaussian_CPD(bnet,n);
end 
disp('Assignment_CPDs: done')

bnet_updt = learn_iparams(bnet, dataE_complete); 
disp 'Learning on the initial network: done'

% calculate the global_score
[global_Score,gs_old] = global_score_calc(nodes_names, bnet_updt,Gene_Expres_list,expr_values);
disp('global_score calculation: done')
%% calculate correlations among TF-TF
% each gene does not change the number of parents; they are fixed from the
% ChIP seq experiments; for this reason the corr is calculated only
% between TF-TF and not among TF-target_genes

for n=1:length(TFnames)
    ind_tfs(n)= find(strcmp(TFnames(n), Gene_Expres_list));
end

evals_tfs= expr_values(ind_tfs,:);
C= corrcoef(evals_tfs');

correlation_whitelist = cell(size(whitelist));
for edge=1:length(whitelist)    
    ind_source= find(strcmp(whitelist{edge,1}, Gene_Expres_list));
    ind_sC= find(ind_source == ind_tfs);
    ind_target = find(strcmp(whitelist{edge,2}, Gene_Expres_list));
    ind_tC= find(ind_target == ind_tfs);
    corr_val = abs(C(ind_sC,ind_tC));
    correlation_whitelist{edge,1} = whitelist{edge,1};
    correlation_whitelist{edge,2} = whitelist{edge,2};
    correlation_whitelist{edge,3} = corr_val;  
end
correlation_whitelist_dagArcs = cell(size(bindingP_tf));
for edge=1:length(bindingP_tf)    
    ind_source= find(strcmp(bindingP_tf{edge,1}, Gene_Expres_list));
    ind_sC= find(ind_source == ind_tfs);
    ind_target = find(strcmp(bindingP_tf{edge,2}, Gene_Expres_list));
    ind_tC= find(ind_target == ind_tfs);
    corr_val = abs(C(ind_sC,ind_tC));
    correlation_whitelist_dagArcs{edge,1} = bindingP_tf{edge,1};
    correlation_whitelist_dagArcs{edge,2} = bindingP_tf{edge,2};
    correlation_whitelist_dagArcs{edge,3} = corr_val;  
end

disp('Correlation: done')

% create copies of original objs
bnet_updt1= bnet_updt;
correlation_whitelist_original =  vertcat(correlation_whitelist_dagArcs, correlation_whitelist);
total_matrix_original = total_matrix;
total_matrix_binding_original = total_matrix_binding;
TFnames_topSorted_original= TFnames_topSorted;

sorted_nodeNames= vertcat(TFnames_topSorted,Gnames);

[new_bnet, CPD_struct] = mk_newBnet(bnet_updt, sorted_nodeNames);
bnet_updt= new_bnet;
bnet_old = bnet_updt;

save Starting_wksp.mat
disp 'workspace saved!'
toc

