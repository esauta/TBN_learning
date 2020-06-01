function assert(pred, str)
% Raise an error if the predicate is not true.

if nargin<2, str = ''; end

if ~pred
  s = sprintf('assertion violated: %s', str);
  error(s);
end
