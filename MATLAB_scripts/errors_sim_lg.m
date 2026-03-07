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

function err = errors_sim_lg(rvarej,ep_pred)

[r2,c2] = size(rvarej);
[r3,c3] = size(ep_pred);

if (r2~=r3)
    error('Matrix sizes do not agree');
end

ep_tot = zeros(r2,1);
for i=1:c3
     ep_tot(:,1) = ep_tot(:,1)+ep_pred(:,i);
end
rvar = rvarej;% + rvarsj;
err = rvar.^0.5;

return;