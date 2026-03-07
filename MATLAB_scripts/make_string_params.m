%  make a string of unknown parameters
% invertedList: list of inverted flow-law parameters - cellstr like {'df.E','df.A','gbs.n','gbs.E'} (file order)
% kn = number of flow-law parameters inverted for (excluding scaling
% coeff.)
% An = number of scaling coeff. inverted for
% Xn = no. of inter-run biases estimated.
function paramstr = make_string_params(invertedList,kn,An,Xn)
    N=kn+Xn+An; % total number of unknown parameters

    % create a string of symbols for the unknown parameters
    paramstr = cell(1,N);
    if kn+An ~=numel(invertedList)
        error('Disagreement between kn+An and number of unknown paramaters declared in invertedList.');
    end
    for i=1:kn+An 
        token = invertedList{i};
        parts = split(token, '.');
        if numel(parts) ~= 2
            error('Bad inverted token "%s". Use format "mech.param", e.g. "df.E".', token);
        end
        mech  = parts{1};
        param = parts{2};
        if strcmp(mech,'df')
            suff = 1;
        elseif strcmp(mech,'ds')
            suff = 3;
        elseif strcmp(mech,'dg')
            suff = 5;
        elseif strcmp(mech,'wf')
            suff = 2;
        elseif strcmp(mech,'ws')
            suff = 4;
        elseif strcmp(mech,'wg')
            suff = 6;
        end
        if strcmp(param,'A')
            paramstr{Xn+i}=['log10' param num2str(suff)];
        else
            paramstr{i}=[param num2str(suff)];
        end
    end
    if Xn>0
        for i=1:Xn
            paramstr{kn+i} = ['X',num2str(i)];
        end
    end
end
