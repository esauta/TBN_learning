function bnet = learn_params(bnet, data)
% Find the maximum likelihood params for a fully observed model
% %
% data(i,m) is the value of node i in case m (can be a cell array)
%
% We set bnet.CPD{i} to its ML/MAP estimate.
%
% Currently we assume no param tying

[n ncases] = size(data);
for j=1:n 
  e = bnet.equiv_class(j);
  assert(e==j);
  if adjustable_CPD(bnet.CPD{e})
    fam = family(bnet.dag,j);
    bnet.CPD{j} = learn_params(bnet.CPD{j}, fam, data, bnet.node_sizes, bnet.cnodes);
  end
end

