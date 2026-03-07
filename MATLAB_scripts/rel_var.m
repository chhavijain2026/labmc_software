% function to help calculate error bars on observed strain rate: dry or wet
% d : microns, sig: Mpa, E: kJ/mol, p: GPa, V: cm^3/mol, T: K, R: kJ/K/mol

% input parameters:
% 1. meandata : matrix containing mean model parameters for flow laws
%   p1, E1, V1, A1: mean model parameters for dry diff. creep. 
%   n3, E3, V3, A3: mean model parameters for dry disl. creep.
%   p2, r2, E2, V2, A2: mean model parameters for wet diff. creep. 
%   r4, E4, E4, V4, A4: mean model parameters for wet disl. creep.
%   sigP0, E5, V5, A5, q1, q2: mean model parameters+exponents for Peierls
%   & Peierls2
%   sigP0, E5, A5, q1, q2, A5: mean model parameters+exponents for Peierls3
%   p6, n6, E6, V6, A6: mean model parameters for dry GBS
%   p7, n7, r7, E7, V7, A7: mean model parameters for wet gbs
%   E ->kJ/mol, V ->cm^3/mol, sigP0 ->GPa
% 2. X : mean inter-run bias
% 3. obsdata : matrix of observed data, to be read as follows:
%   T (K), dT(K), P(GPa), dP(GPa), strain rate(1/s), strain rate error (1/s),
%   stress (MPa), stress error (MPa), grain size (microns), 
%   grain size error (microns), C_OH (ppm H/Si), dC_OH (ppm H/Si), run#
%   T: column vector (K); P: column vector (GPa); sig: column vector (MPa)
%   d: column vector (micron); 
% 4. condition = 'df' for dry diff; 'ds' for dry disl; 'wf' for wet diff; 
%                'ws' for wet disl; 'dg' for dry GBS; 'wg' for wet gbs; 
%                'LP1' for Peierls; 'LP2' for Peierls2; 
%               'LP3' for Peierls3; 'LP4' for Peierls4;
% See labmc code manual for definitions of LP1-LP4

% relvar = total relative variance 
% rvarej = rel. variance of only strain rate
function [relvar,rvarej,ep_pred] = rel_var(meandata,X1,obsdata,condition)

R = 8.314*1e-3; % kJ/K/mol
X = zeros(size(X1)); % I don't know why, but I didn't actually consider X=X1

% categorize experimental data:
T = obsdata(:,1);
P = obsdata(:,3);
ep = obsdata(:,5);
sig = obsdata(:,7);
d = obsdata(:,9);
Cw = obsdata(:,11);
% error in experimental data:
dT = obsdata(:,2);
dP = obsdata(:,4);
dep = obsdata(:,6);
dsig = obsdata(:,8);
dd = obsdata(:,10);
dCw = obsdata(:,12);

% calculate rvar(ep(j))
rvarej = (dep./ep).^2;

% % load model parameters
if strcmp(condition,'df') %condition == 1 % dry dif
    p1 = meandata(:,1); E1 = meandata(:,2); V1 = meandata(:,3); A1 = meandata(:,4);
    % calculate strain rate for dry diffusion creep : ep1
    ep_pred = A1 * (d.^(-p1)).*sig.*exp(-(E1+(P*V1))./(R*T)).*exp(X);
    term2 = (p1 * dd./d).^2;
    term3 = (dsig./sig).^2;
    term4 = ((V1./(R*T)).*dP).^2;
    term5 = (((E1+(P*V1))./(R*(T.^2))) .* dT).^2;
    relvar = term2 + term3 + term4 + term5;
    
elseif strcmp(condition,'ds') %condition == 2 % dry dis
    n3 = meandata(:,1); E3 = meandata(:,2); V3 = meandata(:,3); A3 = meandata(:,4);
    % calculate strain rate for dry dislocation creep : ep3
    ep_pred = A3 * (sig.^n3).*exp(-(E3+(P*V3))./(R*T)).*exp(X);
    term7 = (n3 * dsig./sig).^2;
    term8 = ((V3./(R*T)).*dP).^2;
    term9 = (((E3+(P*V3))./(R*(T.^2))) .* dT).^2;
    relvar = term7 + term8 + term9;
    
elseif strcmp(condition,'dg') %condition == 6 %GBS
    p6 = meandata(:,1); n6 = meandata(:,2); E6 = meandata(:,3); V6 = meandata(:,4);
    A6 = meandata(:,5);
    % calculate strain rate for dry gbs : ep2 <- might not be correct
    ep_pred = A6 * (d.^(-p6)).*(sig.^(n6)).*exp(-(E6+(P*V6))./(R*T)).*exp(X);
    term2 = (p6 * dd./d).^2;
    term3 = (n6 * dsig./sig).^2;
    term4 = ((V6./(R*T)).*dP).^2;
    term5 = (((E6+(P*V6))./(R*(T.^2))) .* dT).^2;
    relvar = term2 + term3 + term4 + term5;
    
elseif strcmp(condition,'wf') %condition == 3 % wet dif
    p2 = meandata(:,1); r2 = meandata(:,2); E2 = meandata(:,3); V2 = meandata(:,4);
    A2 = meandata(:,5);
    % calculate strain rate for wet diffusion creep : ep2
    ep_pred = A2 * (d.^(-p2)).*sig.*(Cw.^(r2)).*exp(-(E2+(P*V2))./(R*T)).*exp(X);
    term2 = (p2 * dd./d).^2;
    term3 = (dsig./sig).^2;
    term4 = (r2 * dCw./Cw).^2; 
    term5 = ((V2./(R*T)).*dP).^2;
    term6 = (((E2+(P*V2))./(R*(T.^2))) .* dT).^2;
    relvar = term2 + term3 + term4 + term5 + term6;
    
elseif strcmp(condition,'ws') %condition == 4 % wet dis
    n4 = meandata(:,2); r4 = meandata(:,1); E4 = meandata(:,3); V4 = meandata(:,4);
    A4 = meandata(:,5);
    % calculate strain rate for wet dislocation creep : ep4
    ep_pred = A4 * (sig.^n4).*(Cw.^(r4)).*exp(-(E4+(P*V4))./(R*T)).*exp(X);
    term8 = (n4 * dsig./sig).^2;
    term9 = (r4 * dCw./Cw).^2;
    term10 = ((V4./(R*T)).*dP).^2;
    term11 = (((E4+(P*V4))./(R*(T.^2))) .* dT).^2;
    relvar = term8 + term9 + term10 + term11;
    
elseif strcmp(condition,'wg') %condition == 10 %wet GBS
    p7 = meandata(:,1); n7 = meandata(:,2); r7 = meandata(:,3); 
    E7 = meandata(:,4); V7 = meandata(:,5);
    A7 = meandata(:,end);
    % calculate strain rate for wet gbs : ep2 <- might not be correct
    ep_pred = A7 * (d.^(-p7)).*(sig.^(n7)).*(Cw.^(r7))...
        .*exp(-(E7+(P*V7))./(R*T)).*exp(X);
    term2 = (p7 * dd./d).^2;
    term3 = (n7 * dsig./sig).^2;
    term4 = (r7 * dCw./Cw).^2;
    term5 = ((V7./(R*T)).*dP).^2;
    term6 = (((E7+(P*V7))./(R*(T.^2))) .* dT).^2;
    relvar = term2 + term3 + term4 + term5 + term6;
    
elseif strcmp(condition,'LP1') %condition == 5 % Peierls
    sigP0 = meandata(:,1)*1e3; %MPa
    E5 = meandata(:,2); V5 = meandata(:,3); A5 = meandata(:,4);
    q1 = meandata(:,5); q2 = meandata(:,6);
    G0 = 77.4; %GPa
    dGp = 1.61;
    dGt = -0.013; % GPa/K
    ep_pred = A5 * (sig.^2).*exp(-((E5+P*V5)./(R*T)).*...
        (1-(sig*G0./(sigP0*(G0+dGp*P))).^q1).^q2).*exp(X);
    termi = (sig*G0./(sigP0*(G0+dGp*P))).^q1;
    termii = (1-termi).^(q2-1);
    termiii = termii.*(1-termi);
    termiv = (E5+P*V5)./(R*T);
    term12 = ((dsig./sig).*(2 + q1*q2*termiv.*termi.*termii)).^2;
    term13 = ((dP./P).*((P*V5./(R*T)).*termiii + ...
        q1*q2*(dGp*P./(G0+dGp*P)).*termiv.*termi.*termii)).^2;
    term14 = ((dT./T).*termiv.*termiii).^2;
    relvar = term12 + term13 + term14;
    
elseif strcmp(condition,'LP2') %condition == 7 %Peierls2
    sigP0 = meandata(:,1)*1e3; %MPa
    E5 = meandata(:,2); V5 = meandata(:,3); A5 = meandata(:,4);
    q1 = meandata(:,5); q2 = meandata(:,6);
    G0 = 77.4; %GPa
    dGp = 1.61;
%     dGt = -0.013; % GPa/K
    ep_pred = A5 *exp(-((E5+P*V5)./(R*T)).*...
        (1-(sig*G0./(sigP0*(G0+dGp*P))).^q1).^q2).*exp(X);
    termi = (sig*G0./(sigP0*(G0+dGp*P))).^q1;
    termii = (1-termi).^(q2-1);
    termiii = termii.*(1-termi);
    termiv = (E5+P*V5)./(R*T);
    term12 = ((dsig./sig).*(q1*q2*termiv.*termi.*termii)).^2;
    term13 = ((dP./P).*((P*V5./(R*T)).*termiii + ...
        q1*q2*(dGp*P./(G0+dGp*P)).*termiv.*termi.*termii)).^2;
    term14 = ((dT./T).*termiv.*termiii).^2;
    relvar = term12 + term13 + term14; 
    
elseif strcmp(condition,'LP3') %condition == 8 % Peierls3
    sigP0 = meandata(:,1)*1e3; %MPa
    E5 = meandata(:,2); A5 = meandata(:,3);
    q1 = meandata(:,4); q2 = meandata(:,5);
    ep_pred = A5 * (sig.^2).*exp(-(E5./(R*T)).*...
        (1-(sig./sigP0).^q1).^q2).*exp(X);
    termi = (sig./sigP0).^q1;
    termii = (1-termi).^(q2-1);
    termiii = termii.*(1-termi);
    termiv = (E5)./(R*T);
    term12 = ((dsig./sig).*(2 + q1*q2*termiv.*termi.*termii)).^2;
    term14 = ((dT./T).*termiv.*termiii).^2;
    relvar = term12 + term14;
    
elseif strcmp(condition,'LP4') %condition == 9 % Peierls4
    sigP0 = meandata(:,1)*1e3; %MPa
    E5 = meandata(:,2); A5 = meandata(:,3);
    q1 = meandata(:,4); q2 = meandata(:,5);
    ep_pred = A5 *exp(-(E5./(R*T)).*...
        (1-(sig./sigP0).^q1).^q2).*exp(X);
    termi = (sig./sigP0).^q1;
    termii = (1-termi).^(q2-1);
    termiii = termii.*(1-termi);
    termiv = (E5)./(R*T);
    term12 = ((dsig./sig).*(q1*q2*termiv.*termi.*termii)).^2;
    term14 = ((dT./T).*termiv.*termiii).^2;
    relvar = term12 + term14;
    
end
    
return;