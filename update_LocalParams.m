function [bnet, CPD_structure] = update_LocalParams(target_node,all_nodes_names_sort_updt,bnet,CPD_structure,expr_values,Gene_Expres_list)
ind_target= find(strcmp(target_node,all_nodes_names_sort_updt)==1);
% target_node must be mapped in this vect with the new topological order!
tarN_parents = parents(bnet.dag,ind_target);

ind_ev_parents=[];
for p=1:length(tarN_parents)
    p_name = all_nodes_names_sort_updt(tarN_parents(p));
    ind_ev_parents(p) = find(strcmp(p_name,Gene_Expres_list)==1);
end
ind_ev_TgNode= find(strcmp(target_node,Gene_Expres_list)==1); 
Y_matrix= zeros(size(expr_values,2),1); %Y column_vect
X_matrix= zeros(size(expr_values,2),size(ind_ev_parents,2)); %regression matrix
for obs=1:size(expr_values,2)
    Y_matrix(obs,1) = expr_values(ind_ev_TgNode,obs);         
    for pev=1:length(ind_ev_parents)            
        X_matrix(obs,pev)= expr_values(ind_ev_parents(pev),obs);                
    end
end

% add a column of ones to the regression_matrix
intercept = ones(1,size(expr_values,2))';
X_matrixi= [intercept X_matrix];

%b = inv(X'X)* X'y
A= (X_matrixi' * X_matrixi);
B= (X_matrixi' * Y_matrix);
Bpred_matrix = A\B; %to find parameters w/o inversion

% cov calculation
% error_Var=(sum(y-(a+bx))^2)/n-k
N= size(expr_values,2);
k_parents = length(tarN_parents);
cov= sum((Y_matrix-(X_matrixi*Bpred_matrix)).^2)/ (N-k_parents);


% update locally the bnet with these new values 
% violate bnet_struct and the target_Node's CPD
ind2map= find(strcmp([CPD_structure.nodes],target_node));
s2updt= CPD_structure(ind2map).CPD;
if ~isempty(tarN_parents)
    s2updt.cps= tarN_parents;
    s2updt.weights = Bpred_matrix(2:size(Bpred_matrix,1),:)';
    s2updt.mean= Bpred_matrix(1,1);
    
else
    s2updt.cps= [];
    s2updt.weights = [];    
    s2updt.mean= Bpred_matrix;    
end
s2updt.cov=cov;

% re assign the updated CPD
CPD_structure(ind2map).CPD= s2updt; 

end
