function [global_Score,gs] = global_score_calc(all_node_names, bnet_updt1,Gene_Expres_list,expr_values)

global_Score.node= 1:length(all_node_names);
global_Score.parents =[];
global_Score.local_RSS = [];
global_Score.local_BIC = [];

n = size(expr_values,2);
BIC_vect=length(all_node_names);
for node=1:length(all_node_names)  
    node_source= all_node_names(node);
    global_Score(node).node = node_source;
    struct_Node= struct(bnet_updt1.CPD{node});
    node_parents = parents(bnet_updt1.dag,node);
    ind_ev_TarNode= find(strcmp(node_source,Gene_Expres_list)==1);
    interc_tar_node= struct_Node.mean;
    vect_sq_resid=[];
    if ~isempty(node_parents)        
        parents_weights= struct_Node.weights;
        ind_ev_parents=[];
        for p=1:length(node_parents)
            global_Score(node).parents = node_parents(p);
            p_name = all_node_names(node_parents(p));
            ind_ev_parents(p) = find(strcmp(p_name,Gene_Expres_list)==1);
        end
        % calculate RSS         
        for obs=1:size(expr_values,2)
            tarN_evidence = expr_values(ind_ev_TarNode,obs);
            vect_ps_evidence=[];
            for pev=1:length(ind_ev_parents)
                ev = expr_values(ind_ev_parents(pev),obs);
                vect_ps_evidence= [vect_ps_evidence ev];                
            end            
            y_pred = (interc_tar_node + sum(vect_ps_evidence.*parents_weights));
            vect_sq_resid(obs)= (tarN_evidence-y_pred)^2;
        end      
    else
        global_Score(node).parents=[];
        for obs=1:size(expr_values,2)
            tarN_evidence = expr_values(ind_ev_TarNode,obs);            
         
            vect_sq_resid(obs)= (tarN_evidence-interc_tar_node)^2;
        end
    end
    % local RSS
    RSS_local= sum(vect_sq_resid);
    global_Score(node).local_RSS = RSS_local;
    % local BIC
    local_BIC = n*log(RSS_local/n)+length(node_parents)*log(n);
    BIC_vect(node)= local_BIC;
    global_Score(node).local_BIC = local_BIC;                                 
end
gs=sum(BIC_vect);%sum of all BIC_scores
     

    
    