function [tfnodes_sorted_updt, all_nodes_names_sort_updt,updated_edges_list,new_totalmatrix, bnet_updt,CPD_structure] = make_New_Model(target_node,bnet_old, CPD_structure, TFedges, TFnames,Gnames,Gedges,expr_values,Gene_Expres_list)

bnet_updt= bnet_old;
% create a topological_ordered vect of TFnames
[tf_tf_sorted, tfnodes_sorted_updt, tf_sort_updt_indexes] = mk_adj_mat_topOrder(TFedges, TFnames, 1);

% create the new_total matrix 
all_nodes_names_sort_updt = ({tfnodes_sorted_updt{:},Gnames{:}}).';  
updated_edges_list = vertcat(TFedges,Gedges);    
dg_total_matr = digraph(updated_edges_list(:,1), updated_edges_list(:,2));

% reorder the matrix depending on the new_sorted_nodeNames_vect (all_nodes_names_sort_updt)
dg_total_new_sorted= reordernodes(dg_total_matr,all_nodes_names_sort_updt);

%% update the structure:
new_totalmatrix = adjacency(dg_total_new_sorted);  

bnet_updt.dag = new_totalmatrix; %update the dag
bnet_updt.names = all_nodes_names_sort_updt;

n = length(new_totalmatrix); %update parents bnet's field
bnet_updt.parents = cell(1,n);
for i=1:n
  bnet_updt.parents{i} = parents(new_totalmatrix, i);
end

% isolate the local_model in which the arc has been added
% target_node: extract name and index of the local target
[bnet_updt, CPD_structure] = update_LocalParams(target_node,all_nodes_names_sort_updt,bnet_updt,CPD_structure,expr_values,Gene_Expres_list);
end
