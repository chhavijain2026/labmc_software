% p = grain-size exp.; n = stress exp.; 
% E = activation energy; V = activation volume; 
% A = scaling coeff.; B = normalized scaling coeff. 
function schema = flowlaw_schema_dry()
% Defines the canonical parameter fields for each mechanism.
schema = struct();
schema.df  = {'p','E','V','A','B'};
schema.ds  = {'n','E','V','A','B'};
schema.dg = {'p','n','E','V','A','B'};
end