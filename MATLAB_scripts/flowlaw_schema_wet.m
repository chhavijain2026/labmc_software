% p = grain-size exp.; n = stress exp.; r = water-content exp.
% E = activation energy; V = activation volume; 
% A = scaling coeff.; B = normalized scaling coeff. 
function schema = flowlaw_schema_wet()
% Defines the canonical parameter fields for each mechanism.
schema = struct();
schema.wf  = {'p','r','E','V','A','B'};
schema.ws  = {'r','n','E','V','A','B'};
schema.wg = {'p','n','r','E','V','A','B'};
end