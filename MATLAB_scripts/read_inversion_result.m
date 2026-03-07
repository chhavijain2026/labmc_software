function model = read_inversion_result(outputFile, mechanisms, invertedList, condition, fixedVals)
% outputFile: MCMC output columns that contain inverted flow-law parameters
% mechanisms: assumed flow laws- cellstr like {'df','dg'}
% invertedList: list of inverted flow-law parameters - cellstr like {'df.E','df.A','dg.n','dg.E'} (file order)
% fixedVals: struct with optional nested structs, e.g.:
%   fixedVals.df.p = 3; fixedVals.df.V = 0;
%   fixedVals.dg.p = 1.4; fixedVals.dg.V = 0;
% conditions = 'dry' or 'wet'

if nargin < 5 || isempty(fixedVals)
    fixedVals = struct();   % create empty struct
end

if strcmp(condition,'dry')
    schema = flowlaw_schema_dry();
else
    schema = flowlaw_schema_wet();
end

% ---------- initialize model with NaNs ----------
allMechs = fieldnames(schema);
for i = 1:numel(allMechs)
    mech = allMechs{i};
    model.(mech) = struct();
    for f = schema.(mech)
        model.(mech).(f{1}) = NaN;
    end
        model.(mech).active = ismember(mech, mechanisms);
        if (model.(mech).active)
            N = numel(fieldnames(model.(mech)));
            model.(mech).means = NaN(1,N-1);
            model.(mech).sdevs = NaN(1,N-1);
        end
end

% model.meta.mechanisms = mechanisms; % mechanisms we assumed
% model.meta.inverted   = invertedList; % parameters we inverted in MCMC
% model.meta.fixed      = fixedVals; % parameters we assumed constant & their values

% ---------- apply fixed values (if provided) ----------
% if nargin >= 5 && ~isempty(fixedVals)
%     mechsFixed = fieldnames(fixedVals);
if ~isempty(fieldnames(fixedVals))
    mechsFixed = fieldnames(fixedVals);
    for i = 1:numel(mechsFixed)
        mech = mechsFixed{i};
        paramsFixed = fieldnames(fixedVals.(mech));
        for j = 1:numel(paramsFixed)
            p = paramsFixed{j};
            model.(mech).(p) = fixedVals.(mech).(p);
        end
    end
end

% ---------- read inversion results from file ----------
% bestfit = readmatrix(outputFile);
bestfit = outputFile;
[nSolutions, nParams] = size(bestfit);

if nParams ~= numel(invertedList)
    error('Mismatch: output file has %d columns but invertedList has %d items.', ...
        nParams, numel(invertedList));
end

% ---------- assign MCMC output values to the right fields ----------
for k = 1:numel(invertedList)
    token = invertedList{k};   % e.g. 'dg.E'
    parts = split(token, '.');
    if numel(parts) ~= 2
        error('Bad inverted token "%s". Use format "mech.param", e.g. "df.E".', token);
    end
    mech  = parts{1};
    param = parts{2};

    % Validate against schema (catches typos)
    if ~isfield(schema, mech) || ~any(strcmp(schema.(mech), param))
        error('Unknown parameter "%s" for mechanism "%s".', param, mech);
    end
    
    model.(mech).(param) = bestfit(:,k);  % assign COLUMN
end

% expand fixed scalars
mechsActive = fieldnames(model);
for i = 1:numel(mechsActive)
    mech = mechsActive{i};
    if ~model.(mech).active, continue; end
    params = fieldnames(model.(mech));
    for j = 1:numel(params)
        p = params{j};
        if strcmp(p,'active'), continue; end
        if strcmp(p,'B'), continue; end
        if strcmp(p,'means'), continue; end
        if strcmp(p,'sdevs'), continue; end
        val = model.(mech).(p);
        if isscalar(val)
            model.(mech).(p) = val * ones(nSolutions,1);
        end
        %model.(mech).mean(j) = mean(val);
    end
end

end
