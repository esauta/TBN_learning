function p = adjustable_CPD(CPD)
% Check if this CPD have any adjustable params (gaussian)

% % This function was adapted from Bayes Net Toolbox written by Kevin Murphy

p = ~CPD.clamped_mean | ~CPD.clamped_cov | ~CPD.clamped_weights;
