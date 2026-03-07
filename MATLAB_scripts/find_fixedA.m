% find which mechanisms have fixed scaling coeff. 
% mechanisms: assumed flow laws- cellstr like {'df','gbs'}
% invertedList: list of inverted flow-law parameters - cellstr like {'df.E','df.A','gbs.n','gbs.E'} (file order)
% An: number of scaling coefficients estimated.
function fixedA = find_fixedA(mechanisms,invertedList,An)

    %nm = numel(mechanisms);
    np = numel(invertedList);
    fieldsm = fieldnames(mechanisms);
    
    nm = 0;
    for i=1:numel(fieldsm)
        mech = fieldsm{i};
        if mechanisms.(mech).active
            nm = nm+1;
        end
    end
    
    if An==0
        fixedA = ones(1,nm);
    elseif An==nm
        fixedA = zeros(1,nm);
    else
        fixedA = ones(1,nm);

        for j = 1:np
            token = invertedList{j};
            parts = split(token, '.');
            if numel(parts) ~= 2
                error('Bad inverted token "%s". Use format "mech.param", e.g. "df.E".', token);
            end
            mechf = parts{1};
            paramf = parts{2};
            if strcmp(paramf,'A')
                for i=1:nm
                    mech = fieldsm{i};
                    if strcmp(mechf,mech)
                        fixedA(i)=0;
                    end
                end
            end
        end
    end

end
