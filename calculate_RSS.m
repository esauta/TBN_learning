function [RSS_new, tarN_parents] = calculate_RSS(target_node,all_nodes_names_sort_updt,bnet_updt, CPD_structure,expr_values,Gene_Expres_list)
% isolate the local_model in which the arc has been added
ind_target_CPD= find(strcmp([CPD_structure.nodes], target_node)==1);
ind_target_bnet= find(strcmp([bnet_updt.names], target_node)==1);
ind_ev_TarN= find(strcmp(target_node,Gene_Expres_list)==1);
% save target_Node's CPD
s_tarNode= CPD_structure(ind_target_CPD).CPD;
tarN_parents = parents(bnet_updt.dag,ind_target_bnet);
interc_tarN= s_tarNode.mean;
if ~isempty(tarN_parents)
    parents_weights = s_tarNode.weights;
    ind_ev_parents=[];
    for p=1:length(tarN_parents)
        p_name = all_nodes_names_sort_updt(tarN_parents(p));
        ind_ev_parents(p) = find(strcmp(p_name,Gene_Expres_list)==1);
    end
    % calculate RSS = sum(obs-y_pred)^2;
    vect_sq_resid=[];
    for obs=1:size(expr_values,2)
        tarN_evidence = expr_values(ind_ev_TarN,obs);
        vect_ps_evidence=[];
        for pev=1:length(ind_ev_parents)
            ev = expr_values(ind_ev_parents(pev),obs);
            vect_ps_evidence= [vect_ps_evidence ev];                
        end
        % calculate y_pred
        y_pred = (interc_tarN + sum(vect_ps_evidence.*parents_weights));
        vect_sq_resid(obs)= (tarN_evidence-y_pred)^2;
    end
else
    vect_sq_resid=[];
    for obs=1:size(expr_values,2)
        tarN_evidence = expr_values(ind_ev_TarN,obs);
        % y_pred = interc_tarN
        vect_sq_resid(obs)= (tarN_evidence-interc_tarN)^2;
    end
end
   
% local RSS
RSS_new= sum(vect_sq_resid);
end