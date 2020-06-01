function CPD = generic_CPD(clamped)
% Virtual constructor for generic CPD
% 
% This function was adapted from Bayes Net Toolbox written by Kevin Murphy

if nargin < 1, clamped = 0; end

CPD.clamped = clamped;
CPD = class(CPD, 'generic_CPD');
