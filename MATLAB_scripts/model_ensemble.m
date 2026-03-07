% perturb model parameters for a mechanism, given means, std
% deviation & correlation coeff.
% this can also handle composite rheologies. 

% M : no. of model ensembles
% kn : no. of unknown model parameters(except scaling coeff. & X)
% An: no. of scaling coeff.
% evpos: vector containing positions of E & V in the means vector
% it need not be just E and V; it is the position of all parameters that
% are not permitted to assume negative values, even if the std. dev. is
% large enough to include negative values in the 68% conf. int. 
% mdry: vector of means
% sdry: vector of standard deviations
% Rdry: corr. matrix
% mBidry: mean Bi vector

% format of mdry or sdry : 
% columns: sigp0, Ep(kJ/mol), Vp(cm3/mol), log10(Ap)
% or for composite rheology, for example,
% columns: p1 E1(kJ/mol) V1(cm3/mol) n3 E3 V3 log10(A1) log10(A3)

function [q_para] = model_ensemble(M,kn,An,evpos,means,sdev,cormat,mBi)

means = means([1:kn end-An+1:end]);
sdev = sdev([1:kn end-An+1:end]);
cormat = cormat([1:kn end-An+1:end],[1:kn end-An+1:end]);

[eigvec,eigval] = eig(cormat); 
q_para = perturb_model(M,kn+An,An,evpos,eigval,eigvec,means,sdev,mBi);
% q_para contains [sigP0 E V Bi] or [sigP0 E Bi] if pressure effects are
% ignored
% for composite rheology, q_para can be, for example, [p1 E1 n3 E3 A1 A3]
return;
