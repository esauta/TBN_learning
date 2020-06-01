function p = find_equiv_posns(vsmall, vlarge)
% 
% THE VECTORS ARE ASSUMED TO BE SORTED.
%
% This function was adapted from Bayes Net Toolbox written by Kevin Murphy
if isempty(vsmall) | isempty(vlarge)
  p = [];
  return;
end
  
bitvec = sparse(1, max(vlarge)); 
bitvec(vsmall) = 1;
p = find(bitvec(vlarge));

