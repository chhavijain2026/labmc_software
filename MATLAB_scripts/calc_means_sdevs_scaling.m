% calculate mean & std. dev. of all scaling coefficients
% scaling coeff. are treated differently from 
% non-scaling parameters [based on Eq. (20) in 
% Korenaga & Karato (2008)].
% 
% mechanisms: assumed flow laws- cellstr like {'df','gbs'}
% meansA,sdevsA = empty 1xAn vectors where 
%               An = number of unknown scaling coeff.
% Bcols = (M x An) matrix, where M=no. of MCMC outputs, 
%               Bcols contains columns of the estimated 
%               values of the scaling coeff. for each
%               flow law, if not fixed, normalized to a 
%               const. T & P [see Eq. (20), KK08]
function [meansA, sdevsA, meansB, Bcols, mechanisms] = ...
    calc_means_sdevs_scaling(mechanisms,fixedA,meansA,sdevsA,Bcols)

j = 0;
ii = 0;
fieldsm = fieldnames(mechanisms);

meansB = meansA;
% Bcols = Acols;

for i = 1:numel(fieldsm)
    mech = fieldsm{i};
    if ~mechanisms.(mech).active, continue; end 
    ii = ii+1;
    % scaling parameter
    E = mechanisms.(mech).E;
    V = mechanisms.(mech).V;
    A = mechanisms.(mech).A;
    [mA,stdA,Bi,mBi] = calc_mAi(A,E,V,mean(E),mean(V));
    mechanisms.(mech).B = Bi;
    mechanisms.(mech).means(end-1) = mA;
    mechanisms.(mech).means(end) = mBi;
    mechanisms.(mech).sdevs(end-1) = stdA;
    mechanisms.(mech).sdevs(end) = stdA;
    if fixedA(ii)==0
        j=j+1;
        meansA(j) = mA;
        sdevsA(j) = stdA;
        meansB(j) = mBi;
        %Acols(:,j) = A;
        Bcols(:,j) = Bi;
    end   
    
end
