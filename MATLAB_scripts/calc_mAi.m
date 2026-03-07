% function to calculate B and mean(B) from equation (20) on pg.7 of KK08
% input parameters:
% Ai : column vector containing estimates of Ai
% Ei : column vector containing estimates of Ei (kJ/mol)
% Vi : column vector containing estimates of Vi (cm3/mol)
% mEi, mVi : mean Ei and mean Vi respectively
% % i corresponds to dry diffusion(1), wet diffusion(2), dry dislocation(3),
                   % wet dislocation(4) 
% Ref.: KK08: Korenaga & Karato (2008), JGR.
function [mAi,sAi,Bi,mBi] = calc_mAi(Ai,Ei,Vi,mEi,mVi)

G2Pa = 1e9; % conversion factor from GPa to Pa
cm2m3 = 1e-6; % conversion factor from cm^3 to m^3
kJ2J = 1e3; % conversion factor from kJ to J

P0 = 0.3; % GPa
T0 = 1523; % K
R = 8.314; % J/K/mol

Ei=Ei*kJ2J; P0=P0*G2Pa; Vi=Vi*cm2m3; mEi=mEi*kJ2J; mVi=mVi*cm2m3;

Bi(:,1) = log10(Ai.*exp(-(Ei+(P0*Vi))/(R*T0))); % eqn(20)
mBi = mean(Bi);

% invert eqn(20) using mean variables to calculate mean(Ai)
mAi = 10^mBi * exp((mEi+(P0*mVi))/(R*T0));

% calculate std. dev. of Ai = std dev of Bi
sBi = (sum((Bi-mBi).^2)/length(Bi))^0.5;
sAi = sBi; % this is standard dev of log10(Ai)

return;





