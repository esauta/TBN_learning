% create a new bnet_struct
function [new_bnet, CPD_struct] = mk_newBnet(bnet, sorted_nodeNames)
    n = length(bnet.dag);
    new_bnet.dag = bnet.dag;
    new_bnet.parents = bnet.parents;
    new_bnet.names = sorted_nodeNames;    
    %new_bnet.order = topological_sort(bnet.dag); % update the nodes' order
    % create another structure to save/update CPDs
    CPD_struct.nodes= 1:length(sorted_nodeNames);
    CPD_struct.CPD = [];
    for n=1:n
        strNode_new = struct;
        strNode_old = struct(bnet.CPD{n});
         
        CPD_struct(n).nodes= sorted_nodeNames(n);
        strNode_new.self = sorted_nodeNames{n};
        strNode_new.cps = bnet.parents{n}; % indicate parents' indexes of each node
        strNode_new.mean = strNode_old.mean;
        strNode_new.cov = strNode_old.cov;
        strNode_new.weights= strNode_old.weights;
        
        CPD_struct(n).CPD= strNode_new;
    end
 
   
end