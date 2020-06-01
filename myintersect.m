function C = myintersect(A,B)
% Intersection of two sets of positive integers 
% 
% This function was adapted from Bayes Net Toolbox written by Kevin Murphy
A = A(:)'; B = B(:)';

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 | mb==0
  C = [];
else
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  C = B(logical(bits(B)));  
end


