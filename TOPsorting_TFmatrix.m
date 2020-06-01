function [TFbinding_Init_matrix] = TOPsorting_TFmatrix(TFnames_topSorted, bindingP_tf)
n=length(TFnames_topSorted);
TFbinding_Init_matrix= zeros(n);
[nr,~] = size(bindingP_tf);
for r=1:nr
  from = strcmp(bindingP_tf{r,1}, TFnames_topSorted);
  assert(~isempty(from));
  to = strcmp(bindingP_tf{r,2}, TFnames_topSorted);
  assert(~isempty(to));
  bs = bindingP_tf{r,3};
  assert(~isempty(bs));
  TFbinding_Init_matrix(from,to) = bs;
end
end