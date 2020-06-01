function [score_eval_model,arc] = arcs_evaluation(arc,n,bnet_old_copy,CPD_struct_copy, TFedges_2update, TFedges_binding_2update,sampling_edges_copy, Gnames, TFnames,Gedges,Gene_Expres_list,expr_values)

% calculate the old_model_structure
target_node_ev = sampling_edges_copy{arc,2};
[~, tfnodes_sorted_old_ev, ~] = mk_adj_mat_topOrder(TFedges_2update, TFnames, 1);

% check if the reverse exists in the sampling_list of this added arc
ind_reverse = find(strcmp(sampling_edges_copy{arc,2},TFedges_2update(:,1))& strcmp(sampling_edges_copy{arc,1},TFedges_2update(:,2)));
ind_reverse_binding = find(strcmp(sampling_edges_copy{arc,2},TFedges_binding_2update(:,1))& strcmp(sampling_edges_copy{arc,1},TFedges_binding_2update(:,2)));    
if ~isempty(ind_reverse) && ~isempty(ind_reverse_binding)
    % calculate the old_RSS
    targetO_node_ev = sampling_edges_copy{arc,1};
    [RSS_old_ev,tarN_parents_old_ev] = calculate_RSS(targetO_node_ev,tfnodes_sorted_old_ev,bnet_old_copy,CPD_struct_copy,expr_values,Gene_Expres_list);
    
    % delete the reverse of the arc
    TFedges_2update(ind_reverse,:) = [];    
    TFedges_binding_2update(ind_reverse_binding,:) = [];    
    % create the eliminated_arc_Model
    D_tf = digraph(TFedges_2update(:,1),TFedges_2update(:,2));
    if ~isdag(D_tf)    
%         disp('not a dag, stop the arc evaluation!') 
        score_eval_model = NaN;
        return
    end       
else
    % calculate the RSS_old_ev on the involved local model
    [RSS_old_ev,tarN_parents_old_ev] = calculate_RSS(target_node_ev,tfnodes_sorted_old_ev,bnet_old_copy,CPD_struct_copy,expr_values,Gene_Expres_list);
end

% fix the arc tested at this iteration into the model
TFedges_2update(end+1,:) = {sampling_edges_copy{arc,1},sampling_edges_copy{arc,2}};

% check if it is DAG
D = digraph(TFedges_2update(:,1),TFedges_2update(:,2));
if ~isdag(D)    
%     disp('not a dag, stop the arc evaluation!') 
    score_eval_model = NaN;
    return
end
[~, all_nodes_names_sort_updt_ev,updated_edges_list,total_matrix_copy, bnet_updt_copy,CPD_struct_copy] = make_New_Model(target_node_ev,bnet_old_copy,CPD_struct_copy,TFedges_2update, TFnames,Gnames,Gedges,expr_values,Gene_Expres_list);
dg_total_matr = digraph(updated_edges_list(:,1), updated_edges_list(:,2));
if ~isdag(dg_total_matr)  
    score_eval_model = NaN;
    return
end            
% calculate RSS on the new model        
[RSS_new_ev,tarN_parents_new_ev] = calculate_RSS(target_node_ev,all_nodes_names_sort_updt_ev,bnet_updt_copy,CPD_struct_copy,expr_values,Gene_Expres_list);
if ~isempty(ind_reverse) && ~isempty(ind_reverse_binding)
    if RSS_new_ev < RSS_old_ev   
        score_eval_model = n*(log(RSS_old_ev - RSS_new_ev));      
    else
%         disp('adding a reverse but RSS_new_ev > RSS_old_ev') 
        score_eval_model = NaN;
        return
    end
else
    if RSS_new_ev < RSS_old_ev         
        score_eval_model = n*(log(RSS_old_ev-RSS_new_ev)) - (log(n));        
    else
%         disp('adding this arc does not improve the local model: RSS_new_ev > RSS_old_ev')
        score_eval_model = NaN;
        return
    end        
end
end
     
                