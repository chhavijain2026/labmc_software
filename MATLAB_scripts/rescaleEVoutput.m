% MCMC output gives E & V in SI units. 
% convert to kJ/mol and cm^3/mol, respectively
% invertedList: list of inverted flow-law parameters - cellstr like {'df.E','df.A','dg.n','dg.E'} (file order)
% outm: MCMC output
function outm = rescaleEVoutput(outm,invertedList)

J2kJ = 1e-3;
m2cm3 = 1e6;

kn = numel(invertedList);

for i=1:kn
    token = invertedList{i};
    parts = split(token, '.');
    if numel(parts) ~= 2
        error('Bad inverted token "%s". Use format "mech.param", e.g. "df.E".', token);
    end
    %mech  = parts{1};
    param = parts{2};
    if strcmp(param,'E')
        outm(:,3+i) = outm(:,3+i)*J2kJ;
    elseif strcmp(param,'V')
        outm(:,3+i) = outm(:,3+i)*m2cm3;
    end
end
