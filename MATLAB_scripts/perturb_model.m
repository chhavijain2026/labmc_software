% function to construct matrix containing ensemble of perturbed model 
% parameters

% input parameters: M (# of model ensembles), Km (# of variable flow law
% parameters), Am(# of scaling coeff), evpos(position of Ei & Vi-vector of positions),
% eigval (diagonal matrix of eigenvalues of R),
% eigvec(eigenvector matrix of R) where R is the correlation matrix, meanm 
% (array of means of 1 non-scaling parameter & 1 scaling parameter), 
% stdm (array containing error bounds on
% the non-scaling parameter & on log10(scaling parameter)], 
% mBi (mean of Bi to replace Ai). 

% meanm : p1 / E1 (kJ/mol)/ V1(cm3/mol)/ n3/ E3(kJ/mol)/ V3(cm3/mol) 
% in meanm and stdm, all scaling coefficients should be put at the end of
% the array. Also, use std. dev. in Bi or log10(Ai) rather than in Ai.

% q_pert contains Km columns - (Km-An) perturbed non-scaling parameter and 
% Am perturbed Bi (not Ai)

function q_pert = perturb_model(M,Km,Am,evpos,eigval,eigvec,meanm,stdm,mBi,varargin)

% R = eigvec*eigval*eigvec'; 
% R is semi+ve definite => inv(eigvec)=eigvec' & all eigenvalues>0

pd = makedist('Normal');
P = zeros(M,Km);
for i=1:M
    for j=1:Km
        P(i,j) = random(pd);
%         while P(i,j)>1 % added on 10 October, 2018
%             P(i,j) = random(pd);
%         end;
%         if j==1 
%             P(i,j)
%         end
    end
end
% IP = (1/(M-1))*(P'*P); % just to check that IP ~ Identity when M is v large

Prand = P*(eigval.^0.5)*eigvec; % correlated random no.s
% IP = (1/(M-1))*(Prand'*Prand); % check that IP ~ R when M is v large

q_pert = zeros(M,Km); % matrix for perturbed model parameters 
for k=1:Km % total no. of unknown model parameters
    if k<=Km-Am %
        q_pert(:,k) = meanm(k)+stdm(k)*Prand(:,k);
    else % specifically scaling constants
        q_pert(:,k) = mBi(k-Km+Am)+stdm(1,k)*Prand(:,k);
    end
end

if ~isempty(varargin) && varargin{1}==0
    return;
end

if ~isempty(evpos)
% remove models for which Vi <0 and add that many models again
    q_pert = checkEV3(q_pert,M,Km,Am,evpos,eigval,eigvec,meanm,stdm,mBi);
end

% % convert Bi to Ai
% if strcmp(condition,'dry')
%         A(:,1) = AifromBi(q_pert(:,Km-1),q_pert(:,2),q_pert(:,3));
%         A(:,2) = AifromBi(q_pert(:,Km),q_pert(:,5),q_pert(:,6));
% end;
% if Km==10
%        A(:,1) = AifromBi(q_pert(:,Km-1),q_pert(:,3),q_pert(:,4));
%        A(:,2) = AifromBi(q_pert(:,Km),q_pert(:,7),q_pert(:,8));
% end;
% q_pert(:,Km-1:Km) = [A(:,1) A(:,2)];

return;