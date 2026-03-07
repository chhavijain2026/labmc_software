% construct confidence intervals about mean or median model
% cint = confidence interval (in %)
% expdata = [T P sig d] for wet or [T P sig d Cw] for wet conditions
% flows = list of mechanisms we have assumed.
% mechanisms: assumed flow laws & their structure- cellstr like {'df','gbs'}
% M = number of MCMC models.
function [epsL, epsU] = construct_CI(cint,expdata,flows,mechanisms,M)
    
    R = 8.3144;

    LL = 0.5 - cint/200; % lower limit
    UL = 0.5 + cint/200; % upper limit
    
    % number of flow laws assumed
    nm = numel(flows);
    
    % dry or wet?
    mech = flows{1};
    if strcmp((mech(1)),'d')
        condition = 'dry';
    elseif strcmp(flows{1}(1),'w')
        condition = 'wet';
    else
        error('flow law declaration start with a "d" or "w"');
    end
            
    [dn,cn] = size(expdata);
    
    
    if strcmp(condition,'wet')
        if cn<5
            error('Conditions are insufficient.');
        else
            Cw = expdata(:,5);
        end
    end
    
    epsU = zeros(dn,1);
    epsL = zeros(dn,1);
    
    if strcmp(condition,'dry')
        for j = 1:dn
            T = expdata(j,1);
            P = expdata(j,2);
            sig = expdata(j,3);
            d = expdata(j,4);
            
            ep = zeros(M,1);
            emech = zeros(M,1);
            for i = 1:nm
                mech = flows{i};
                if strcmp(mech,'df')
                    emech = mechanisms.(mech).A.*(sig).*...
                        (d.^-mechanisms.(mech).p).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                elseif strcmp(mech,'ds')
                    emech = mechanisms.(mech).A.*...
                        (sig.^mechanisms.(mech).n).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                elseif strcmp(mech,'dg')
                    emech = mechanisms.(mech).A.*...
                        (d.^-mechanisms.(mech).p).*...
                        (sig.^mechanisms.(mech).n).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                end
                % total strain rate at the given conditions 
                % predicted by all M MCMC models
                ep = ep+emech;
            end
            ep = sort(ep); % order ep
            epsU(j,1) = ep(round(UL*M)); % upper limit prediction at given conditions
            epsL(j,1) = ep(round(LL*M)); % lower limit prediction at given conditions
        end
        
    elseif strcmp(condition,'wet')
        for j = 1:dn
            T = expdata(j,1);
            P = expdata(j,2);
            sig = expdata(j,3);
            d = expdata(j,4);
            Cw = expdata(j,5);
            
            ep = zeros(M,1);
            emech = zeros(M,1);
            for i = 1:nm
                mech = flows{i};
                if strcmp(mech,'wf')
                    emech = mechanisms.(mech).A.*(sig).*...
                        (d.^-mechanisms.(mech).p).*...
                        (Cw.^mechanisms.(mech).r).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                elseif strcmp(mech,'ws')
                    emech = mechanisms.(mech).A.*...
                        (sig.^mechanisms.(mech).n).*...
                        (Cw.^mechanisms.(mech).r).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                elseif strcmp(mech,'wg')
                    emech = mechanisms.(mech).A.*...
                        (d.^-mechanisms.(mech).p).*...
                        (sig.^mechanisms.(mech).n).*...
                        (Cw.^mechanisms.(mech).r).*...
                        exp(-(mechanisms.(mech).E + P*mechanisms.(mech).V)*1e3./(R*T));
                end
                % total strain rate at the given conditions 
                % predicted by all M MCMC models
                ep = ep+emech;
            end
            ep = sort(ep); % order ep
            epsU(j,1) = ep(round(UL*M)); % upper limit prediction at given conditions
            epsL(j,1) = ep(round(LL*M)); % lower limit prediction at given conditions
        end
end