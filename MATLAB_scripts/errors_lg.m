% function to calculate error bars on observed strain rate: dry or wet
% after calling function rel_var in the main script
% input variables:
% 1. relvar : matrix in which no. of rows = no. of data points, no. of
% columns = no. of flow laws; contains relvar (output of rel_var) from each
% flow law combined into one matrix
% 2. rvarej : column vector containing value (dep/ep).^2 : output of
% rel_var
% 3. ep_pred : matrix in which no. of rows = no. of data points, no. of
% columns = no. of flow laws; contains predicted edot (output of rel_var) 
% from each flow law combined into one matrix

function err = errors_lg(relvar,rvarej,ep_pred)

[r1,c1] = size(relvar);
[r2,c2] = size(rvarej);
[r3,c3] = size(ep_pred);

if (r1~=r2 || r1~=r3 || c1~=c3)
    error('Matrices size do not agree');
end

rvarsj = zeros(r1,1);

if c1==1
    rvarsj = relvar;
else
    ep_tot = zeros(r1,1);
    for i=1:c1
        ep_tot(:,1) = ep_tot(:,1)+ep_pred(:,i);
    end
    for j=1:r1
        tmp = 0;
        for i=1:c1
            tmp = tmp + ((ep_pred(j,i)/ep_tot(j,1))^2).*relvar(j,i);
        end
        rvarsj(j,1) = rvarsj(j,1)+tmp;
    end
end
rvar = rvarej + rvarsj;
err = rvar.^0.5;

return;