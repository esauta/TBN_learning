function [A, names,order] = mk_adj_mat_topOrder(connections, names, topological)
% Make a directed adjacency matrix from a list of connections between named nodes.
%
% 
% The last argument  indicates that we should topologically sort the nodes (parents before children).
% the graph is only possible if it has no directed cycles.

if nargin < 3, topological = 0; end
  
n=length(names);
A=zeros(n);
[nr nc] = size(connections);
for r=1:nr       
  from = strmatch(connections{r,1}, names, 'exact');
  assert(~isempty(from));
  to = strmatch(connections{r,2}, names, 'exact');
  assert(~isempty(to));
%   fprintf(1, 'from %s %d to %s %d\n', connections{r,1}, from, connections{r,2}, to);
  A(from,to) = 1;
end

if topological
  order = topological_sort(A); 
  A = A(order, order); 
  names = names(order); 
end


