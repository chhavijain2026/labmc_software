% try to draw deformation map using MCMC inversion results
% Use results from case S7fs to -
% 1) make an ensemble of perturbed models 
% 2) calculate 90% & 50% conf. intervals about mean. 
% 3) check if using all models versus using only models with
% misfit<200 produces the same results or not. It doesn't 
% because when we truncate models, 
% the distributions are not fully gaussian anymore. 

clear
clc

kJ2J = 1e3;
cm2m3 = 1e-6;
R=8.314;

% grain growth parameters:
k0 = 3.8e-9; %m2/s
Eg = 160; %kJ/mol

path = '../CaseS7/';
fname = 'out.S7fs.resampled.dat';
dataf = 'dataS7.dat';
obs = load([path dataf]);
params = [{'p1'};{'n3'}]; % list of unknown params
kn = 2; % 2 estimated flow-law parameters
An = 2; % 2 estimated scaling coeff.
Xn = 0; % no inter-run bias


%%%%%%%%%%%%%%%%%%%
% MCMC output file
outm = load([path fname]);

% calculate scaling coefficients
E1 = zeros(size(outm(:,4))); %E1=0
V1 = zeros(size(outm(:,4))); %V1=0
A1 = outm(:,end-1);
[mA1,stdA1,Bi1,mBi1] = calc_mAi(A1,E1,V1,mean(E1),mean(V1));
E3 = zeros(size(outm(:,4))); %E3=0
V3 = zeros(size(outm(:,4))); %V3=0
A3 = outm(:,end);
[mA3,stdA3,Bi3,mBi3] = calc_mAi(A3,E3,V3,mean(E3),mean(V3));
% % inversion results, after running datafit_dry_flowlaw.m
% means = [1.9932 3.0345 -4.7061 -10.8518]; %+/-[0.1585 0.2453 0.1466 0.5235]
% mp1 = means(1); mE1 = 0; mV1 = 0; mn1=1;
% mn3 = means(2); mE3 = 0; mV3 = 0;

% % calculate mean, std and R. 
N=kn+Xn+An;
means = zeros(1,N);
sdevs = zeros(1,N);
sdevs2 = zeros(1,N);
cormat = zeros(N,N);
for i=1:kn+Xn
    means(i) = mean(outm(:,3+i));
    sdevs(i) = std(outm(:,3+i));
    sdevs2(i) = std(outm(:,3+i),1);
end
if An==2
    means(end-1:end) = [mBi1 mBi3];
    sdevs(end-1:end) = [stdA1 stdA3];
    sdevs2(end-1:end) = [stdA1 stdA3];
    Acols = outm(:,end-1:end);
    outm(:,end-1:end) = [Bi1 Bi3];
elseif An==1
    means(end) = [mBi1];
    sdevs(end) = [stdA1];
    sdevs2(end) = [stdA1];
    Acols = outm(:,end);
    outm(:,end) = [Bi1];
end
for i=1:N
    for j=1:N
        Rcor = mean((outm(:,i+3)-means(i)).*(outm(:,j+3)-means(j)))/(sdevs2(i)*sdevs2(j));
        cormat(i,j) = Rcor;
    end
end
outm(:,end-1:end) = Acols;
means(end-1:end) = [log10(mA1) log10(mA3)];

% mean parameters
mp1 = means(1); mE1 = 0; mV1 = 0; mn1=1;
mn3 = means(2); mE3 = 0; mV3 = 0;

% perturb mean model to construct an ensemble of perturbed models 
evpos = []; % column numbers at which E & V are located, if at all
mBi = [mBi1 mBi3];
q_para = model_ensemble(100000,kn,An,evpos,means,sdevs,cormat,mBi);
q_para = [[means(1:2) mBi1 mBi3]; q_para];
M = length(q_para(:,1));

% experimental conditions
Td=obs(:,1);
Pd=obs(:,3);
edotd=obs(:,5);
dedotd=obs(:,6);
sigd=obs(:,7);
dd=obs(:,9);

% declare constant reference conditions
sig0=0.1;
d0=1e3;
P0=Pd(1);
T0=Td(1);
% Earth conditions?
sigE = [0.1 0.3 0.5 0.6 0.8 1 2 200]'; %MPa
dE = [1e-3 2e-3 0.05 0.1 0.2 0.5 0.8 1 2 3 5 8 10]'*1e3; %0.1-10mm in microns
% confidence intervals
UL = 90; %90% confidence intervals
UL2 = 50; %50% confidence intervals
%     UL2 = 20; %20% confidence intervals
% UL2 = 60; %60% confidence intervals

misfit_vec=NaN(1,M);

for m=1:M
    % Diffusion (dry):
    p_dif = q_para(m,1);
    E_dif = mE1; 
    V_dif = mV1; 
    A_dif = AifromBi(q_para(m,3),E_dif,V_dif); % A_dif = 10^q_para(m,end-1);

    % Dislocation (dry):
    n_dis = q_para(m,2); %q_para(m,3);
    E_dis = mE3; %465; %q_para(m,4);
    V_dis = mV3; %15; %0;
    A_dis = AifromBi(q_para(m,4),E_dis,V_dis); % A_dis = 10^q_para(m,end);

    % raw strain rate predictions
    edot_pred_diff = A_dif*(sigd).*(dd.^(-p_dif)).*exp(-(E_dif+Pd*V_dif)*1e3./(R*Td));
    edot_pred_disl = A_dis*(sigd.^n_dis).*exp(-(E_dis+Pd*V_dis)*1e3./(R*Td));
    edot_pred = edot_pred_diff + edot_pred_disl;
    misfit_vec(m)= mean((log(edotd)-log(edot_pred)).^2./(dedotd./edotd).^2);

    % normalize data
    % normalized by randomized model
    ef = A_dif*(sigd).*(d0.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
    es = A_dis*(sigd.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
    eps(:,m) = ef + es;
    if m==1
        ep_norms = (edotd./edot_pred).*eps(:,m);
    end
    %
    ef = A_dif*(sig0).*(dd.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
    es = A_dis*(sig0.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
    epd(:,m) = ef + es;
    if m==1
        ep_normd = (edotd./edot_pred).*epd(:,m); 
    end
        
    % extrapolate to Earth conditions?
    eef = A_dif*(sigE).*(d0.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
    ees = A_dis*(sigE.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
    eeps(:,m) = eef+ees;
    eef = A_dif*(sig0).*(dE.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
    ees = A_dis*(sig0.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
    eepd(:,m) = eef+ees;

end

% draw confidence intervals
M=length(q_para(:,1));
epsm = eps(:,1); % mean model
epdm = epd(:,1); % mean model
eepsm = eeps(:,1); % mean model
eepdm = eepd(:,1); % mean model
epsmv = [eepsm; epsm];
epdmv = [eepdm; epdm];
sigv = [sigE; sigd];
dv = [dE; dd];
epsv = [eeps; eps];
epdv = [eepd; epd];
for j=1:length(sigv)
    epss = sort(epsv(j,:));
    epsmed(j,1) = epss(round(50*M/100)); % median
    %90% CI
    epsUL(j,1) = epss(round((50+UL/2)*M/100));
    epsLL(j,1) = epss(round((50-UL/2)*M/100));
    % 50% CI
    epsUL2(j,1) = epss(round((50+UL2/2)*M/100));
    epsLL2(j,1) = epss(round((50-UL2/2)*M/100));
end
[sigd1,epsmv,epsmed,epsLL,epsUL,epsLL2,epsUL2]=rearrange1(sigv,epsmv,epsmed,epsLL,epsUL,epsLL2,epsUL2);
for j=1:length(dv)
    epds = sort(epdv(j,:));
    epdmed(j,1) = epds(round(50*M/100)); % median
    %90% CI
    epdUL(j,1) = epds(round((50+UL/2)*M/100));
    epdLL(j,1) = epds(round((50-UL/2)*M/100));
    % 50% CI
    epdUL2(j,1) = epds(round((50+UL2/2)*M/100));
    epdLL2(j,1) = epds(round((50-UL2/2)*M/100));
end
[dd1,epdmv,epdmed,epdLL,epdUL,epdLL2,epdUL2]=rearrange1(dv,epdmv,epdmed,epdLL,epdUL,epdLL2,epdUL2);
edifd = mA1*(sig0).*(dd1.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
edisd = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
% etotd = edifd+edisd;
edifs = mA1*(sigd1).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
ediss = mA3*(sigd1.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
% etots = edifs+ediss;
figure
subplot(2,1,1)
hold on
fill_area(sigd1,epsLL,epsUL,[0.75 0.75 0.75],0.3,'90%')
fill_area(sigd1,epsLL2,epsUL2,[0.5 0.5 0.5],0.3,'50%')
plot(sigd,ep_norms,'ko')
hs=plot(sigd1,epsmv,'b',sigd1,epsmed,'b--',sigd1,edifs,'r',sigd1,ediss,'g');
set(hs,'LineWidth',2)
xlabel('stress [MPa]')
ylabel('strain rate [s^{-1}]')
set(gca,'XScale','log','YScale','log','FontSize',16)
box on
hold off
subplot(2,1,2)
hold on
fill_area(dd1,epdLL,epdUL,[0.75 0.75 0.75],0.3,'90%')
fill_area(dd1,epdLL2,epdUL2,[0.5 0.5 0.5],0.3,'50%')
plot(dd,ep_normd,'ko')
hd=plot(dd1,epdmv,'b',dd1,epdmed,'b--',dd1,edifd,'r',dd1,edisd*ones(size(dd1)),'g');
set(hd,'LineWidth',2)
xlabel('grain size [\mu m]')
ylabel('strain rate [s^{-1}]')
set(gca,'XScale','log','YScale','log','FontSize',16)
box on
hold off

nbins=100;
figure
plothist(misfit_vec,100,'\chi^2_M/N',0,[3,1,1])
plothist(q_para(:,1),nbins,'p_1',0,[3,1,2])
plothist(q_para(:,2),nbins,'n_3',0,[3,1,3])
title('Histogram of perturbed models')


% return

%Amongst perturbed models, choose those that produce a very small misfit
    smallchi=find(misfit_vec<3); %0.5); %7); %10);
    qs=q_para(smallchi,:);
    misf=misfit_vec(smallchi);
    [mmisf,pval]=rearrange1(misf,[1:length(misf)]');
    qs = qs(pval,:);
    for m=1:length(qs(:,1))
        A_dif(m,1) = AifromBi(qs(m,end-1),300,10);
        A_dis(m,1) = AifromBi(qs(m,end),465,15);
    end
    qs = [qs A_dif A_dis mmisf'];
    if ~isempty(qs)
        
        for m=1:2%length(qs(:,1))
            % Diffusion (dry):
            p_dif = qs(m,1);
            E_dif = 300; %q_para(m,2);
            V_dif = 10; %0;
            A_dif = AifromBi(qs(m,3),E_dif,V_dif); % A_dif = 10^q_para(m,end-1);

            % Dislocation (dry):
            n_dis = qs(m,2); %q_para(m,3);
            E_dis = 465; %q_para(m,4);
            V_dis = 15; %0;
            A_dis = AifromBi(qs(m,4),E_dis,V_dis); % A_dis = 10^q_para(m,end);

            % raw strain rate predictions
            edot_pred_diff = A_dif*(sigd).*(dd.^(-p_dif)).*exp(-(E_dif+Pd*V_dif)*1e3./(R*Td));
            edot_pred_disl = A_dis*(sigd.^n_dis).*exp(-(E_dis+Pd*V_dis)*1e3./(R*Td));
            edot_pred = edot_pred_diff + edot_pred_disl;
            
            % normalized by randomized model
            ef = A_dif*(sigd).*(d0.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
            es = A_dis*(sigd.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
            eps = ef + es;
            ep_norms = (edotd./edot_pred).*eps;
            %
            ef = A_dif*(sig0).*(dd.^(-p_dif)).*exp(-(E_dif+P0*V_dif)*1e3./(R*T0));
            es = A_dis*(sig0.^n_dis).*exp(-(E_dis+P0*V_dis)*1e3./(R*T0));
            epd = ef + es;
            ep_normd = (edotd./edot_pred).*epd;
            
            % normalized by n=3.2 (no covariance)
            edotf2 = mA1*(sigd).*(dd.^(-mp1)).*exp(-(mE1+Pd*mV1)*1e3./(R*Td));
            edots2 = mA3*(sigd.^n_dis).*exp(-(mE3+Pd*mV3)*1e3./(R*Td));
            edotp2 = edotf2 + edots2;
            %
            ef2 = mA1*(sigd).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
            es2 = mA3*(sigd.^n_dis).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
            eps2 = ef2 + es2;
            ep_norms2 = (edotd./edotp2).*eps2;
            %
            ef2 = mA1*(sig0).*(dd.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
            es2 = mA3*(sig0.^n_dis).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
            epd2 = ef2 + es2;
            ep_normd2 = (edotd./edotp2).*epd2;
            
            if m==1
                pred_datas = [sigd ep_norms eps];
                pred_datad = [dd ep_normd epd];
            elseif m==2
                pred_datas = [pred_datas ep_norms eps ep_norms2 eps2];
                pred_datad = [pred_datad ep_normd epd ep_normd2 epd2];
            end

            figure
            subplot(2,2,1)
            plot(sigd,ep_norms,'ko',sigd,eps,'b')
            xlabel('\sigma [MPa]')
            ylabel('strain rate [s^{-1}]')
            set(gca,'FontSize',16,'XScale','log','YScale','log')
            subplot(2,2,2)
            plot(dd,ep_normd,'ko',dd,epd,'b')
            xlabel('d [\mu m]')
            ylabel('strain rate [s^{-1}]')
            set(gca,'FontSize',16,'XScale','log','YScale','log')
            subplot(2,2,3)
            plot(sigd,ep_norms2,'ko',sigd,eps2,'b')
            xlabel('\sigma [MPa]')
            ylabel('strain rate [s^{-1}]')
            set(gca,'FontSize',16,'XScale','log','YScale','log')
            subplot(2,2,4)
            plot(dd,ep_normd2,'ko',dd,epd2,'b')
            xlabel('d [\mu m]')
            ylabel('strain rate [s^{-1}]')
            set(gca,'FontSize',16,'XScale','log','YScale','log')
        end

    end



return
