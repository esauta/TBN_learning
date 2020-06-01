
rng('shuffle') %seeds the random number generator based on the current time
format long
w=warning('off','all');
conf

if learning
    Initial_learning
    disp('initial_workspace created!')
end

conf

for ind=1:length(n_run)
    counter_run = n_run(ind);
    disp(['the ', num2str(counter_run),' run starts'])
    path2save = strcat(path_def, fold_name,'_',num2str(counter_run),'\');
    mkdir(path2save);
    globalwk_name = strcat(num2str(counter_run),'run_');
    reportfile= strcat('Report_',run_date,'_',num2str(counter_run),'.txt');
    logfile=strcat('log_',run_date,'_',num2str(counter_run),'.txt');  
    % initialize vectors
    Added_arcs_vect = {};
    reversed_added_arcs_vect = {};
    vect_dag_scores = [];    
    RSS_new=0;
    RSS_old=0;
    delta_gBIC=0;
    gs_prov= -5.63579412686210e+04; %has to be initialized depending on threshold set in conf
    dag_score_old = 2;
    dag_score_new= 3;   
    ngiro=1;
    resampling_cnt=0;
    osc_param= 1;

    % Initialize report file and log file
    fid= fopen(reportfile, 'w');
    fprintf(fid, 'Iteration\tOld_gBIC\tNew_gBIC\tDelta_gBIC\tlocal_deltaBIC_old\tlocal_deltaBIC_new\tDelta_delta\tOld_RSS\tNew_RSS\tDelta_RSS\n');
    fclose(fid);
    logf= fopen(logfile, 'w');
    
    load('Starting_wksp.mat')      
    fprintf(logf, 'BIC of the initial_model: ');
    fprintf(logf, '%.15e\n\n', gs_old);   
    n = size(expr_values,2);

    while 1
        if (ngiro == max_iter)
            disp('Exiting because max_iter reached')
            break
        end
        if gs_prov > gs_old
            if resampling_cnt==0 || resampling_cnt <= 5
                disp('Try to resampling..')
                resampling_cnt = resampling_cnt+1;
                disp(['Resampling: ',num2str(resampling_cnt)])
                disp(gs_prov)            
                numb_to_extract= 100;
                %numb_to_extract= 10; %set for example data
                fprintf(logf, 'Try to resampling.. \n');

           elseif (5 < resampling_cnt) && (resampling_cnt < 10)
               numb_to_extract= 150;  
               %numb_to_extract= 15; %set for example data
               disp('Try to increase sampling_size and test more arcs..')
               fprintf(logf,['Try to increase sampling_size and ' ...
                              'test more arcs..\n']);
               resampling_cnt = resampling_cnt+1;
               disp(['Resampling with more arcs: ', ...
                      num2str(resampling_cnt)])
               disp(gs_prov)
               f_res =strcat('globalwkspace_',num2str(ngiro), ...
                              '_Resampling_',num2str(resampling_cnt),'.mat');
               save(fullfile(path2save,f_res))

           elseif resampling_cnt == 10
               disp('After 10 consecutive trials,')
               disp('Exiting: old model is better')
               disp(gs_prov)
               f_last_wksp = strcat('globalwkspace_',num2str(ngiro),'_exiting.mat');
               save(fullfile(path2save,f_last_wksp))
               break
           end
        else           
            numb_to_extract= 100;
%             numb_to_extract= 5; %set for example data
            gs_new = gs_prov;
            resampling_cnt=0;
        end
        delta = (dag_score_new - dag_score_old)/dag_score_old;
        delta_bic = (gs_new - gs_old)/gs_old;
        delta_gBIC = [delta_bic; delta_gBIC];
    %     disp(['new_gBIC: ', num2str(gs_new)])
    %     disp(['old_gBIC: ', num2str(gs_old)])
    %     disp(['delta_gBIC: ', num2str(delta_bic)])
    %     disp(['delta_localBIC: ', num2str(delta)])  

        if (abs(delta_bic))~=0
            if abs(delta_bic) < threshold
                if osc_param == max_osc
                    disp('Exiting because threshold reached')
                    f_last_wksp_ext = strcat('globalwkspace_',num2str(globalwk_name),num2str(ngiro),'_exiting.mat');
                     save(fullfile(path2save,f_last_wksp_ext))
                    break
                else
                    osc_param= osc_param+1;
                end
            else
                osc_param= 1;
            end
        end

        fid= fopen(reportfile, 'a');
        fprintf(fid, '%i\t%.15e\t%.15e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n', ngiro, ...
                gs_old,gs_new,delta_bic,dag_score_old, dag_score_new, delta, RSS_old, RSS_new);
        fclose(fid);
        fprintf(logf, 'giro: %i\n\n',ngiro);
        fprintf(logf, 'archi da testare: %i\n', numb_to_extract);
        fprintf(logf, 'model_gScore_new: %.15e\n',gs_new);
        fprintf(logf, 'model_gScore_old: %.15e\n',gs_old);
        fprintf(logf, 'gScore_difference: %.15e\n\n', delta_bic);

        dag_score_old = dag_score_new;
        disp(['arcs to test: ', num2str(numb_to_extract)])            
        shuffling = randperm(length(correlation_whitelist)); 
        new_correl_whitelist= correlation_whitelist(shuffling,:);
        % sampling w/o replacement
        sampling_edges= datasample(new_correl_whitelist,min(numb_to_extract,length(new_correl_whitelist)),'Replace',false, 'Weights',[new_correl_whitelist{:,3}]);

        scores_vect =[]; 
        ind_arcs= [];
        % create copy of vects used in the parfor
        TFedges_2update= TFedges;
        TFedges_binding_2update = bindingP_tf;
        CPD_struct_copy = CPD_struct;
        sampling_edges_copy = sampling_edges;
        bnet_old_copy = bnet_old;
        parfor arc=1:length(sampling_edges)        
            [score_evalModel]= arcs_evaluation(arc,n,bnet_old_copy,CPD_struct_copy,TFedges_2update,TFedges_binding_2update,sampling_edges_copy,Gnames,TFnames,Gedges,Gene_Expres_list,expr_values);
            scores_vect(arc) = score_evalModel;
            ind_arcs(arc)= arc;
        end        
        disp 'end arcs evaluation!'       
    %     if ~any(scores_vect) %check if all scores are zeros
    %         disp 'all model_scores are zeros!'
    %         disp(scores_vect)
    %         continue
    %     end
        [best_score_model, ind_best]= max(scores_vect);
        if ind_arcs(ind_best) ~= ind_best
            break
        end
        disp(best_score_model)
        vect_dag_scores = [vect_dag_scores,best_score_model];

        % save results    
        fprintf(logf, 'best_model_score(delta_BIC): ');
        fprintf(logf, '%.5e\n', best_score_model);

        Added_arcs_vect = [Added_arcs_vect; {sampling_edges{ind_best,:}}];
    %     disp('this arc has been added: '); disp(sampling_edges(ind_best,:))
    %     disp(['numb. of added_arcs: ',num2str(size(Added_arcs_vect,1))])    
        fprintf(logf,'best_arc: ');
        fprintf(logf, '%s\t%s\t%.5f\n\n', sampling_edges{ind_best,:});   

        disp('updating the model..')
        target_node = sampling_edges{ind_best,2};
        % add the best_arc into the model_arcs  
        TFedges(end+1,:) = {sampling_edges{ind_best,1},sampling_edges{ind_best,2}};
        ind_best_all_arcList = find(strcmp(sampling_edges{ind_best,1},all_arcs_list(:,1))& strcmp(sampling_edges{ind_best,2},all_arcs_list(:,2)));
        best_bs_2add= all_arcs_list{ind_best_all_arcList,3};
        bindingP_tf(end+1,:) = {sampling_edges{ind_best,1},sampling_edges{ind_best,2},best_bs_2add};

        % check if a reverse_arc was added  
        ind_reverse_best = find(strcmp(sampling_edges{ind_best,2},TFedges(:,1))& strcmp(sampling_edges{ind_best,1},TFedges(:,2)));
        ind_reverse_best_binding = find(strcmp(sampling_edges{ind_best,2},bindingP_tf(:,1))& strcmp(sampling_edges{ind_best,1},bindingP_tf(:,2)));    
        if ~isempty(ind_reverse_best) && ~isempty(ind_reverse_best_binding)% se esiste il reverse del best_arc 
            targetO_node = sampling_edges{ind_best,1};
            % delete the reverse arc that will be added into the model
            TFedges(ind_reverse_best,:) = [];
            bs_add = bindingP_tf{ind_reverse_best_binding,3};
            bindingP_tf(ind_reverse_best_binding,:) = [];    
            % create the eliminated_arc_Model
            D_tf = digraph(TFedges(:,1),TFedges(:,2));
            if ~isdag(D_tf)    
                disp('not a dag, check it!')    
                break
            end   
            disp 'a reverse in the model has been added!'
            % mapping this arc with the correct correlation_value,add it to the corr_whitelist
            ind_corr = find(strcmp(sampling_edges{ind_best,2}, correlation_whitelist_original(:,1)) & strcmp(sampling_edges{ind_best,1},correlation_whitelist_original(:,2)));
            corr_value = correlation_whitelist_original{ind_corr,3};
            correlation_whitelist(end+1,:)= {sampling_edges{ind_best,2},sampling_edges{ind_best,1},corr_value};            
            reversed_added_arcs_vect(end+1,:) = {sampling_edges{ind_best,2},sampling_edges{ind_best,1},bs_add};   

            fprintf(logf, 'a reverse in the model has been added!\n\n');
            fprintf(logf, 'a reverse of this arc: %s\t%s\t', sampling_edges{ind_best,1:end-1});
            fprintf(logf, 'with this binding_prob: %.5f\n', sampling_edges{ind_best,3});
            % update the global_Score struct with the new info
            disp 'updating the global_Score struct..'
           % update the CPD_struct and the global_score    
           [tfnodes_sorted_updt, all_nodes_names_sort_updt,updated_edges_list,total_matrix, bnet_updt,CPD_struct] = make_New_Model(target_node,bnet_old,CPD_struct,TFedges, TFnames,Gnames,Gedges,expr_values,Gene_Expres_list);
           [bnet_updt,CPD_struct] = update_LocalParams(targetO_node,all_nodes_names_sort_updt,bnet_updt,CPD_struct,expr_values,Gene_Expres_list);            
           [RSS_oN, tarN_Oparents] = calculate_RSS(targetO_node,all_nodes_names_sort_updt,bnet_updt,CPD_struct,expr_values,Gene_Expres_list);

            ind_targONode= find(strcmp([global_Score.node],targetO_node)==1);
            global_Score(ind_targONode).parents = tarN_Oparents;
            global_Score(ind_targONode).local_RSS = RSS_oN;

            local_OBIC = n*log(RSS_oN/n)+ log(n)*length(tarN_Oparents);       
            global_Score(ind_targONode).local_BIC= local_OBIC;  
        else
            [tfnodes_sorted_updt, all_nodes_names_sort_updt,updated_edges_list,total_matrix, bnet_updt,CPD_struct] = make_New_Model(target_node,bnet_old,CPD_struct,TFedges, TFnames,Gnames,Gedges,expr_values,Gene_Expres_list);
        end     
        [RSS_new,tarN_parents_new] = calculate_RSS(target_node,all_nodes_names_sort_updt,bnet_updt,CPD_struct,expr_values,Gene_Expres_list);

        % delete this arc from the correlation_whit
        ind_corrWhitelist= find(strcmp(sampling_edges{ind_best,1},correlation_whitelist(:,1))& strcmp(sampling_edges{ind_best,2},correlation_whitelist(:,2)));
        correlation_whitelist(ind_corrWhitelist,:) = [];    

        % update the global_Score struct with the new local model
        ind_targNode= find(strcmp([global_Score.node],target_node)==1);
        global_Score(ind_targNode).parents = tarN_parents_new;    
        global_Score(ind_targNode).local_RSS = RSS_new;

        local_nBIC = n*log(RSS_new/n)+log(n)*length(tarN_parents_new);  
        global_Score(ind_targNode).local_BIC= local_nBIC;  
        % update the gs calculation
        gs_old= gs_new;
        gs_prov = sum([global_Score.local_BIC]);      

        dag_score_new = best_score_model;
        % save a global wkspace
        f_global_wksp = strcat('globalwkspace_',num2str(ngiro),'.mat');
        save(fullfile(path2save,f_global_wksp))
        ngiro = ngiro+1;
        disp(['the ', num2str(ngiro),' iteration starts'])    

        bnet_old = bnet_updt;        
    end
    fclose(logf);

end
fclose('all');

