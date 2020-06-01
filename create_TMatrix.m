function [TF_starting_matrix, total_matrix, total_matrix_binding, TFnames_topSorted,topological_Order_indexes] = create_TMatrix(TFedges,TFnames, bindingP_tf, bindingP_Genes, Gedges, Gnames)
[TF_starting_matrix, TFnames_topSorted,topological_Order_indexes] = mk_adj_mat_topOrder(TFedges, TFnames, 1); 
% 1 means that the nodes are sorted in a topological order (parents before children)
 
%create the binding matrix, observing the topological order of the
    %starting_matrix
TFbinding_Init_matrix =  TOPsorting_TFmatrix(TFnames_topSorted,bindingP_tf);

recycle_matrix=sparse(zeros(size(TFbinding_Init_matrix,1),length(Gnames)));

for r=1:length(Gedges)
  fr = find(strcmp(Gedges{r,1}, TFnames_topSorted));  
  to = find(strcmp(Gedges{r,2}, Gnames));
  bp = bindingP_Genes{r,3};
  recycle_matrix(fr,to) = bp;
end
m_dimens = length(TFnames) + length(Gnames);
total_matrix_binding = sparse(m_dimens,m_dimens);
total_matrix_binding(1:length(TFnames),:) = [TFbinding_Init_matrix,recycle_matrix];

total_matrix= total_matrix_binding;
total_matrix (total_matrix> 0) = 1;

end

    
