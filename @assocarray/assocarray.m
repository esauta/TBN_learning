function A = assocarray(keys, vals)
% Make an associative array
% %
% keys{i} is the i'th string, vals{i} is the i'th value.
% After construction, A('foo') will return the value associated with foo.

A.keys = keys;
A.vals = vals;
A = class(A, 'assocarray');
