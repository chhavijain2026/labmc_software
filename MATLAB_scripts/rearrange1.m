% function to rearrange elements in 'var' vector in ascending order. Related
% vectors are re-organized in accordance to the 'var' vector.

function [var,varargout] = rearrange1(var,varargin)

n=length(var);

if nargin ~= nargout
    error('Error: no. of inputs & no. of outputs not same');
end;
if nargin < 1
    error('Error: min. 1 & max. 8 input arguments necessary');
end;

nOutputs = nargout-1;
% varargout = cell(1,nOutputs);

for index=1:n
    p=index;
    tmp = var(index); 
    for i=index+1:n
        if var(i)<tmp;
            tmp=var(i);
            p=i;
        end;
    end;
    var(p)=var(index);
    var(index)=tmp;
    for k=1:nOutputs
        tmp=varargin{k}(index); 
        varargin{k}(index)=varargin{k}(p); 
        varargin{k}(p)=tmp;
    end;
end;
varargout = varargin;
return;