function p = adjustable_CPD(CPD)
% Check if this CPD have any adjustable params (generic)
% % This function was adapted from Bayes Net Toolbox written by Kevin Murphy

   
p = ~CPD.clamped;
