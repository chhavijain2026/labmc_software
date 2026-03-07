% calculate mean & std. dev. of all parameters except 
% scaling coefficients. If a parameter is assumed constant, 
% mean = fixed_val and std. dev. = 0; 
% scaling coeff. are treated separately. 
% 
% mechanisms: assumed flow laws- cellstr like {'df','gbs'}
% meansp, sdevsp = empty 1xkn vectors where 
%                           kn = number of unknown flow-law
%                           parameters, excluding the scaling coeff.
function [meansp, sdevsp, mechanisms] = calc_means_sdevs_params(mechanisms,meansp,sdevsp)

j = 0;
fieldsm = fieldnames(mechanisms);

for i = 1:numel(fieldsm)
    mech = fieldsm{i};
    if ~mechanisms.(mech).active, continue; end
    params = fieldnames(mechanisms.(mech));
    np = numel(params)-3; % ignore 'active','mean','std' fields
    
    for nn=1:np-2 %only consider non-scaling parameters
        pname = params{nn};
        val = mechanisms.(mech).(pname);
        mval = mean(val);
        sval = std(val);
        mechanisms.(mech).means(nn) = mval;
        mechanisms.(mech).sdevs(nn) = sval;
        if (val==mean(val)), continue; end % fixed parameter
        j = j+1;
        meansp(j) = mval; % only for unfixed parameters
        sdevsp(j) = sval;            
    end
end
