function sub = block(blocks, block_sizes)
% Return a vector of subscripts corresponding to the specified blocks.
% 
% e.g., block([2 5], [2 1 2 1 2]) = [3 7 8].
% % This function was adapted from Bayes Net Toolbox written by Kevin Murphy

blocks = blocks(:)';
block_sizes = block_sizes(:)';
skip = [0 cumsum(block_sizes)];
start = skip(blocks)+1;
fin = start + block_sizes(blocks) - 1;
sub = [];
for j=1:length(blocks)
  sub = [sub start(j):fin(j)];
end
