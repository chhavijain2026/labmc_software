% function to calculate Ai from Bi
% Ei : kJ/mol; Vi : cm^3/mol;

function Ai = AifromBi(Bi,Ei,Vi)

G2Pa = 1e9; % conversion factor from GPa to Pa
cm2m3 = 1e-6; % conversion factor from cm^3 to m^3
kJ2J = 1e3; % conversion factor from kJ to J

% reference state
p0 = 0.3; % GPa
T0 = 1523; % K
R = 8.314; % J/K/mol

Ei=Ei*kJ2J; p0=p0*G2Pa; Vi=Vi*cm2m3;

% invert eqn(20) using input variables to calculate Ai
Ai = (10.^Bi) .* exp((Ei+(p0*Vi))/(R*T0));

return;