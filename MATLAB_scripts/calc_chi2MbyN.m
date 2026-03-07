% calculate chi^2_M/N from the MCMC output
% 
% obsdata : input experimental data matrix of dn x 13 size in the format:
%       [T dT P dP edot de sig dsig d dd Cw dCw runs]
%       where dn : number of data points. runs is a 
%       vector of run identifiers from 1 to Nr, where  
%       Nr is the number of unique run numbers.
% outm : MCMC output
% mechanisms: assumed flow laws- cellstr like {'df','gbs'}

function chi2MbyN = calc_chi2MbyN(mechanisms,obsdata,outm,kn)

dn = length(obsdata(:,1));
N = length(outm(:,1));

fieldsm = fieldnames(mechanisms);
nm = 0; % number of active mechanisms
for j=1:numel(fieldsm)
    if mechanisms.(fieldsm{j}).active
        nm=nm+1;
    end
end
runs = obsdata(:,end);
Xn = length(unique(runs));
if Xn==1
    Xn = 0;
end
    
% recalculate chi^2_M/N 
chi2MbyN = zeros(N,1);
Xvec = zeros(dn,1);
relvar = zeros(dn,nm);
ep_pred = zeros(dn,nm);
for i=1:N
    if Xn==0
        Xvec = zeros(dn,1); 
    else
        for ii=1:dn
            Xvec(ii,1) = outm(i,3+kn+runs(ii));
        end
    end
    jj = 0;
    for j=1:numel(fieldsm)
        mech = fieldsm{j};
        if ~mechanisms.(mech).active, continue; end
        jj = jj+1;
        if strcmp(mech,'dg')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).p(i) mechanisms.(mech).n(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        elseif strcmp(mech,'df')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).p(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        elseif strcmp(mech,'ds')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).n(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        elseif strcmp(mech,'wg')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).p(i) mechanisms.(mech).n(i) mechanisms.(mech).r(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        elseif strcmp(mech,'wf')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).p(i) mechanisms.(mech).r(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        elseif strcmp(mech,'ws')
            [relvar(:,jj),rvarej,ep_pred(:,jj)] = ...
                rel_var([mechanisms.(mech).r(i) mechanisms.(mech).n(i) mechanisms.(mech).E(i) mechanisms.(mech).V(i) mechanisms.(mech).A(i)],Xvec,obsdata,mech);
        end
    end
    denom = rvarej; %(dedot./edot).^2; %err.^2;
    eptot = sum(ep_pred,2).*exp(Xvec);    
    num = (log(obsdata(:,5))-log(eptot)).^2;
    chi2MbyN(i,1) = mean(num./denom);
end

end