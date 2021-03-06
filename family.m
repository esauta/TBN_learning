function f = family(A,i,t)
% Return the indices of parents and self in sorted order
%
% t is an optional argument: if present, dag is assumed to be a 2-slice DBN
% This function was adapted from Bayes Net Toolbox written by Kevin Murphy

if nargin < 3 
  f = [parents(A,i) i];
else
  if t == 1
    f = [parents(A,i) i];
  else
    ss = length(A)/2;
    j = i+ss;
    f = [parents(A,j) j] + (t-2)*ss;
  end
end
