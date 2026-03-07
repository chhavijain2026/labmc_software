% Analyze data S9 for only a single & composite flow laws (wet)
% - gbs, diffusion, or dislocation.
% Include inter-run bias where assumed.
% 1. Resample output files at intervals determined from ACF
% 2. Combine resampled parallel runs
% 3. Determine scaling coefficients from MCMC output
% 4. Calculate mean, std. dev., correlations. Write in a file
% 5. Draw pdf of each unknown parameter.
% 6. Draw normalized strain rate graphs to test data fit
% GBS creep - invert for p6, n6, r6, E6, V6.
% Diffusion creep - invert for p2, r2, E2, V2. 
% Dislocation creep - invert for r4, n4, E4, V4.

% Provide the following information- 
% 1. path = path to input & output data files (e.g., './')
% 2. dataf = name of input data file <filename.dat>
% 3. outf = name of MCMC output file (without '.dat')
% 4. flows = array of the identifiers for each of the   
%     flow laws assumed, provided  in the same order 
%     in which they were declared in the MCMC file: 
%       'wf' - wet diffusion
%       'ws' - wet dislocation
%       'wg' - wet GBS
%     E.g., If parameter file was declared with -Pc<> -Pg<>,
%           then An = 2 and flow = {'wf', 'wg'};
% 5. An = the number of scaling coefficients inverted for
% 6. kn = total number of flow-law parameters inverted for 
%     (excluding the scaling co-efficients).
% 7. invertparams = An array of parameters being inverted for, in the same
%     sequence in which they appear in the MCMC output file
%     {'wf.p','wf.r','wf.E','wf.V','wf.A'} - wet diffusion creep
%     {'ws.r','ws.n','ws.E','ws.V','ws.A'} - dislocation creep
%     {'wg.p','wg.r','wg.n','wg.E','wg.V','wg.A'} - gbs creep
%          E.g., If parameter file was declared with :
%               -Pc0/0/1/3/-1/1/10/10/0/0 -Pg<0/0/1/3/0.2/0.2/1/5/0/0/0/0,
%          then flow = {'wf', 'wg'}; An = 2; kn = 4; invertparams = {'wf.p', 'wf.r', 'wg.p', 'wg.n', 'wf.A', 'wg.A'};
% 8. nfixed = no. of fixed flow-law parameters. E.g., nfixed = 5 in the prev. example.
% 9. fixed.(<mech>).<para> = fixed_val : 
%     parameter para of mechanism mech whose value 
%     is assumed to be contant = fixed_val.
%     Call individually for each such parameter (para = p,n,r,E,V,A) 
%     for each mechanism (mech = wf,ws,wg) if nfixed>1
%           E.g., nfixed = 5; fixed.wf.E=10; fixed.wf.V=10; fixed.wg.r=0.2; fixed.wg.E=0; fixed.wg.V=0;
%     NOTE: if A is fixed, provide actual value of A and not log10(A).
% 10. XXn = no. of runs in the original data set
% 11. Xn = number of inter-run biases estimated.
% 12. parallelruns = Number of parallel simulations conducted
% 13. resample = Resampling interval for the MCMC output based on the 
%     autocorrelation function (ACF) 
% 14. printresults = enable/disable writing the results in files.
%           0: disable, 1: enable
% 15. (OPTIONAL) para0: vector of "expected values" of each unknown parameter (not
%       scaling coeff.), if known. comment if not known.
% 16. (OPTIONAL) paraX0: Same as para0 but for inter-run biases.

clc
clear

J2kJ = 1e-3;
m2cm3 = 1e6;
R = 8.3144;
cols={[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [1 1 0]; [0 1 1]; ...
    [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]};
grey = [0.75 0.75 0.75];
grey2 = [0.5 0.5 0.5];
ci = 1; % plot confidence intervals or not (90% & 50%)
printresults = 0; % write means etc. or not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case S9fsX = Wet data, f+s 
path = '../CaseS9/';
dataf = 'dataS9.dat'; 
outf = 'out.S9fsX.'; 
flows = {'ws','wf'}; 
An = 2;
kn = 8; %p2,r2,E2,V2,r4,n4,E4,V4
invertparams ={'ws.r','ws.n','ws.E','ws.V','wf.p','wf.r','wf.E','wf.V','ws.A','wf.A'}; % list of unknown params
nfixed = 0;
Xn = 12; % 
XXn = 12; % 
parallelruns = 2; 
resample = 6000; 
para0 = [1.2 3.5 480 15 2 1 330 5 1e2 5E3];
paraX0 = [-1 -0.5 0 0.5 1 -0.7 0.7 -0.3 0.3 0 0.2 -0.2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dry or wet?
if strcmp(flows{1}(1),'d')
    %condition = 'dry';
    error('This script is only for wet flow laws.');
elseif strcmp(flows{1}(1),'w')
    condition = 'wet';
else
    error('flow laws must start be "d" or "w".');
end

resampledout = [outf 'resampled.'];
covfile = [outf 'R_cov.txt']; % file in which covariance matrix is printed
mname = [outf 'means.']; %
nname = outf; %

N=kn+Xn+An; % total number of unknown parameters
% if An>0
    meansA=NaN(1,An);
    sdevsA = NaN(1,An);
% end
% if kn>0
    meansp=NaN(1,kn);
    sdevsp = NaN(1,kn);
% end
% if Xn>0
    meansX=NaN(1,Xn);
    sdevsX = NaN(1,Xn);
% end
% means = zeros(1,N);
% sdevs = zeros(1,N);
covmat = zeros(N,N);

% create a string of symbols for the unknown parameters
paramstr = make_string_params(invertparams,kn,An,Xn);

% parallel runs- output files
for i=0:parallelruns-1
    fname = strcat(outf,num2str(i),'.dat');
    files{i+1,:} = fname;
end

% % % % % % % % % % % % % % % 
% comment out this block if the resampled MCMC output has been saved. 
% Simply load that saved file instead. 
% % % % % % % % % % % % % % % 
% resample MCMC outputs from all parallel runs & compile into one
outm = combine_data(path,files,resample); % resampled & compiled data
% remove models for which the scaling coefficients are nan, inf, or 0
for a=0:An-1
    ff = find(isnan(outm(:,end-a)) | isinf(outm(:,end-a)) | outm(:,end-a)==0);
    outm(ff,:) = [];
end
% remove models for which the chi2 is unnaturally large
f4 = find(outm(:,2)>1e15);
outm(f4,:) = [];
clearvars ff f4
if printresults ==1
    % write resampled file
    fn=fopen([path outf 'resampled.dat'],'w');
    fspec = '%d %8.6e %8.6e';
    for i=4:length(outm(1,:))
        fspec = [fspec ' %8.6e'];
    end
    fspec = [fspec '\n'];
    fprintf(fn,fspec,outm');
    fclose(fn);
end
% % % % % % % % % % % % % % % % % % % 

% % % % uncomment if the resampled MCMC output has been saved.
% outm=load([path outf 'resampled.dat']);
% % % uncomment if the resampled MCMC output has been saved.

Nout = length(outm(:,1));
Nmech = numel(flows);

% convert J/mol to  kJ/mol &  m^3/mol to cm^3/mol
if (kn+An)~=numel(invertparams), error('kn+An~=length(invertparams.'); end
outm = rescaleEVoutput(outm,invertparams);

% construct & initialize the assumed rheological model
% read MCMC output and assign columns to the appropriatte variables
% assign assumed values to the fixed variables, if any
if nfixed > 0
    model = read_inversion_result(outm(:,[4:kn+3 4+kn+Xn:end]), flows, invertparams, condition, fixed);
else
    model = read_inversion_result(outm(:,[4:kn+3 4+kn+Xn:end]), flows, invertparams, condition);
end

% % determine which scaling coefficients are fixed
fixedA = find_fixedA(model,invertparams(kn+1:end),An);
Acols = outm(:,end-An+1:end); 
Bcols = NaN(size(Acols)); %NaN(Nout,An);

% calculate mean & std. dev. of all parameters 
% if a parameter is assumed constant, mean = fixed_val and std. dev. = 0;
% scaling coeff. are treated differently from non-scaling parameters
[meansp, sdevsp, model] = calc_means_sdevs_params(model,meansp,sdevsp);
% calculate mean & std. dev. of the unknown scaling coefficients 
% A = scaling coeff.; B = normalized scaling coeff.
% Based on Equation (20) on pg.7 of Korenaga & Karato (2008).  
[meansA, sdevsA, meansB, Bcols,model] = calc_means_sdevs_scaling(model,fixedA,meansA,sdevsA,Bcols);

% calculate mean inter-run biases
if Xn>0
    for i=1:Xn
        meansX(1,i) = mean(outm(:,3+kn+i));
        sdevsX(1,i) = std(outm(:,3+kn+i));
    end
end

% % Calculate parameter correlations, R. 
% % Write them in a text file
means = [meansp meansX meansB];
sdevs = [sdevsp sdevsX sdevsA];
outm(:,end-An+1:end) = Bcols;
for i=1:N
    for j=1:N
        Rcor = mean((outm(:,i+3)-means(i)).*(outm(:,j+3)-means(j)))/(sdevs(i)*sdevs(j));
        covmat(i,j) = Rcor;
    end
end
outm(:,end-An+1:end) = Acols;
means(end-An+1:end) = log10(meansA);
if printresults==1
    fn=fopen([path mname 'dat'],'w');
    fspec = '%8.6e %8.6e\n';
    fprintf(fn,fspec,[means' sdevs']');
    fclose(fn);

    fn = fopen([path covfile],'w'); %,'a');
    fspec = '%3s';
    for i=1:kn+Xn+An
        fspec = [fspec '\t\t%6.3f'];
    end
    fspec = [fspec ' \t\t%6.3f +/- %6.3f\n'];
    fprintf(fn,'%s\n',['Inversion (' outf '.dat) of data (' dataf '): correlations, means, std']);
    fprintf(fn, '\t\t\t');
    for i=1:N
        fprintf(fn,'%s\t\t',paramstr{i});
    end
    fprintf(fn,'%s +/- %s\n','mean','std');
    for i=1:N
        fprintf(fn,fspec,paramstr{i},covmat(i,:),means(i),sdevs(i));  
    end
    fclose(fn);
    clearvars Rcor
end

% draw normalized datafit figures
% read observed data
obsdata = load([path dataf]);
[dn,colobs] = size(obsdata);
T=obsdata(:,1);
dT=obsdata(:,2);
P=obsdata(:,3);
dP=obsdata(:,4);
edot=obsdata(:,5);
dedot=obsdata(:,6);
sig=obsdata(:,7);
dsig=obsdata(:,8);
d=obsdata(:,9);
dd=obsdata(:,10);
if strcmp(condition,'wet')
    Cw=obsdata(:,11);
    dCw=obsdata(:,12);
end
runs=obsdata(:,end);     
uniqruns = unique(runs,'stable');
if Xn==0
    runindex = ones(dn,1);
else
    runindex = NaN(dn,1);
    for j=1:length(uniqruns)
        f=find(runs==uniqruns(j));
        runindex(f,1) = j*ones(length(f),1);
    end
end
if colobs<13
    obsdata = [obsdata(:,1:10) zeros(dn,2) runindex];
else
    obsdata = [obsdata(:,1:12) runindex];
end
if Xn~=XXn %Xn==0 && XXn>0
    obsdatax = load([path dataf2]);
    runsorig=obsdatax(:,end);
    clearvars obsdatax
    uniqruns = unique(runsorig,'stable');
    runsx = NaN(dn,1);
    for j=1:length(uniqruns)
        f=find(runsorig==uniqruns(j));
        runsx(f,1) = j*ones(length(f),1);
    end
elseif Xn>0
    runsx=runs;
end

%recalculate chi^2_M/N 
chi2M_vec = calc_chi2MbyN(model,obsdata,outm,kn);
% out_dry2 = [outm(:,1:2) chi2M_vec outm(:,4:end)];
% % fn=fopen([path fname 'chi2M.dat'],'w');
% % fspec = '%d %8.6e %8.6e';
% % for i=4:length(out_dry2(1,:))
% %     fspec = [fspec ' %8.6e'];
% % end
% % fspec = [fspec '\n'];
% % fprintf(fn,fspec,out_dry2');
% % fclose(fn);
% % clearvars out_dry2 fn fspec

nbins = 50;
% plot histogram of estimated flow-law parameters
titlestr = [condition 'flow law: ' flows{1}];
if Nmech>1
    for i=2:Nmech
        titlestr = [titlestr ' + ' flows{i}];
    end
end
figure
title({['Inversion (' outf ') : flow parameters'],titlestr},'FontSize',14)
for i=-1:kn
    if i==-1
        plothist(outm(:,2),nbins,'simple \chi^2_J/N_p',0,[3,4,1]);
    elseif i==0
        plothist(chi2M_vec,nbins,'\chi^2_M/N',0,[3,4,2]);
    else 
        plothist(outm(:,3+i),nbins,paramstr{i},1,[3,4,i+2]);
        if exist('para0')
            hold on
            plot([para0(i) para0(i)],[0 500],'b--')
            hold off
        end
    end
end
% plot histogram of inter-run bias
if Xn>0
    figure
    title({['Inversion (' outf ') : X'],titlestr},'FontSize',14)
    for i=1:Xn
        plothist(outm(:,3+kn+i),nbins,paramstr{kn+i},1,[3,4,i]); %['X_' num2str(i)],1,[3,3,i]);
        if exist('paraX0')
            hold on
            plot([paraX0(i) paraX0(i)],[0 500],'b--')
            hold off
        end
    end
end

% correct observed strain rates for inter-run bias
Xmat = zeros(dn,1);
edotX = zeros(dn,1);
for i=1:dn
    if Xn>0
        Xmat(i) = meansX(runindex(i)); %means(kn+runindex(i));
    else
        Xmat(i) = 0;
    end
    edotX(i)=edot(i)/exp(Xmat(i));
end

% mean parameters: 
if model.wf.active
    mp1 = model.wf.means(1); mr1 = model.wf.means(2); 
    mE1 = model.wf.means(3); 
    mV1 = model.wf.means(4); mA1 = model.wf.means(5);
end
if model.ws.active
    mr3 = model.ws.means(1); mn3 = model.ws.means(2); 
    mE3 = model.ws.means(3); 
    mV3 = model.ws.means(4); mA3 = model.ws.means(5);
end
if model.wg.active
    mp5 = model.wg.means(1); mn5 = model.wg.means(2); 
    mr5 = model.wg.means(3); mE5 = model.wg.means(4); 
    mV5 = model.wg.means(5); mA5= model.wg.means(6);
end
    
% strain rates predictions and error bars at experimental conditions
% using mean model parameters
edot_pred_diff = zeros(dn,1);
edot_pred_disl = zeros(dn,1);
edot_pred_gbs = zeros(dn,1);
relvar = zeros(dn,1);
ep_pred = zeros(dn,1);
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:5),Xmat,obsdata,mech);
        edot_pred_diff = ep_pred(:,i);
    end
    if strcmp(mech,'ws')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:5),Xmat,obsdata,mech);
        edot_pred_disl = ep_pred(:,i);
    end
    if strcmp(mech,'wg')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:6),Xmat,obsdata,mech);
        edot_pred_gbs = ep_pred(:,i);
    end
end
edot_pred = edot_pred_diff + edot_pred_disl + edot_pred_gbs;
err = errors_sim_lg(rvarej,ep_pred);
lerr = edotX - exp(log(edotX)-err);
uerr = exp(log(edotX)+err)-edotX;

% reference & extended conditions
if min(sig)==max(sig)
    sig0 = min(sig);
else
    sig0 = 80;
    sig3 = linspace(min(sig)/2,max(sig)*2,100)';
end
if min(d)==max(d)
    d0 = min(d);
else
    d0 = 8; % microns
    d3 = linspace(min(d)/2,max(d)*2,100)'; 
end
if min(T)==max(T)
    T0 = min(T);
else
    T0 = 1523; %1573; %K
    T3 = linspace(min(T)-100,max(T)+100,100)';
end
if min(P)==max(P)
    P0 = min(P);
else
    P0 = 1.5; % GPa
    P3 = linspace(0.1,max(P)+1,100)'; 
end
if min(Cw)==max(Cw)
    Cw0 = min(Cw);
else
    Cw0 = 500; % GPa
    Cw3 = linspace(min(Cw)/2,max(Cw)*2,100)'; 
end


% plot data fits

hfig=figure;

ep_pred_diff = zeros(dn,1);
ep_pred_disl = zeros(dn,1);
ep_pred_gbs = zeros(dn,1);
epf = zeros(length(sig3),1);
eps = zeros(length(sig3),1);
epg = zeros(length(sig3),1);

% normalize at reference T0, P0, d0
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        ep_pred_diff = mA1*(sig).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
        epf = mA1*(sig3).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
    elseif strcmp(mech,'ws')
        ep_pred_disl = mA3*(sig.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*(Cw0.^mr3);
        eps = mA3*(sig3.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*(Cw0.^mr3);
    elseif strcmp(mech,'wg')
        ep_pred_gbs = mA5*(sig.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
        epg = mA5*(sig3.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
    end
end
ep_preds = ep_pred_disl + ep_pred_diff + ep_pred_gbs;
epsex = epf + eps + epg;
% normalize observed strain rate
ep_norms = (edotX./edot_pred).*ep_preds;
errs = err;
uerrs = exp(log(ep_norms)+errs)-ep_norms;
lerrs = ep_norms - exp(log(ep_norms)-errs);
if ci==1
% calculate confidence intervals
[epsL, epsU] = construct_CI(90,[T0*ones(size(sig3)) P0*ones(size(sig3)) sig3 d0*ones(size(sig3)) Cw0*ones(size(sig3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(sig3)) P0*ones(size(sig3)) sig3 d0*ones(size(sig3)) Cw0*ones(size(sig3))],flows,model,Nout); % 50% conf. intervals
end
subplot(2,3,1)
hold on
if ci==1
    fill_area(sig3,epsL,epsU,grey,0.3,'90% CI');
    fill_area(sig3,epsL2,epsU2,grey2,0.3,'50% CI')
end
if Xn==0 && XXn==0
    errorbar(sig,ep_norms,lerrs,uerrs,'ko','DisplayName','data')
else %if XXn>0
    for xx=1:XXn
        hh=errorbar(sig(runsx==xx),ep_norms(runsx==xx),lerrs(runsx==xx),uerrs(runsx==xx),'o','DisplayName',['run' num2str(xx)]);
        %set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh=plot(sig3,epsex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.wf.active, plot(sig3,epf,'r--','DisplayName','diffusion'); end
    if model.ws.active, plot(sig3,eps,'g--','DisplayName','dislocation'); end
    if model.wg.active, plot(sig3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log','XScale','log')
xlabel('\sigma [MPa]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable stress, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(d0) '\mu m, ' ...
    num2str(T0) 'K' num2str(Cw0) 'ppm H/Si']},'FontSize',14)
hold off
% normdatatable = [runs,sig,ep_norms,lerrs,uerrs,ep_preds,ep_pred_diff,ep_pred_disl,ep_pred_gbs];
% preddatatable = [sig3 epsex epf eps epg];


if min(d)<max(d)
epf = zeros(length(d3),1);
eps = zeros(length(d3),1);
epg = zeros(length(d3),1);
% normalize at reference T0, P0, sig0
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        ep_pred_diff = mA1*(sig0).*(d.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
        epf = mA1*(sig0).*(d3.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
    elseif strcmp(mech,'ws')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*ones(size(d)).*(Cw0.^mr3);
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*ones(size(d3)).*(Cw0.^mr3);
    elseif strcmp(mech,'wg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
        epg = mA5*(sig0.^mn5).*(d3.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
    end
end
ep_predd = ep_pred_disl + ep_pred_diff + ep_pred_gbs;
epdex = epf + eps + epg;
% normalize observed strain rate
ep_normd = (edotX./edot_pred).*ep_predd;
% % calculate error 
errd = err;
uerrd = exp(log(ep_normd)+errd)-ep_normd;
lerrd = ep_normd - exp(log(ep_normd)-errd);
if ci==1
% calculate confidence intervals
[epsL, epsU] = construct_CI(90,[T0*ones(size(d3)) P0*ones(size(d3)) sig0*ones(size(d3)) d3 Cw0*ones(size(d3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(d3)) P0*ones(size(d3)) sig0*ones(size(d3)) d3 Cw0*ones(size(d3))],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,3,2)
hold on
if ci==1
    fill_area(d3,epsL,epsU,grey,0.3,'90% CI')
    fill_area(d3,epsL2,epsU2,grey2,0.3,'50% CI')
end
if Xn==0 && XXn==0
    errorbar(d,ep_normd,lerrd,uerrd,'ko','DisplayName','data')
else %if XXn>0
    for xx=1:XXn
        hh=errorbar(d(runsx==xx),ep_normd(runsx==xx),lerrd(runsx==xx),uerrd(runsx==xx),'o','DisplayName',['run' num2str(xx)]);
        %set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(d3,epdex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.wf.active, plot(d3,epf,'r--','DisplayName','diffusion'); end
    if model.ws.active, plot(d3,eps,'g--','DisplayName','dislocation'); end
    if model.wg.active, plot(d3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log','XScale','log')
xlabel('d [\mum]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable grain size, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(T0) 'K, ' ...
    num2str(sig0) 'MPa' num2str(Cw0) 'ppm H/Si']},'FontSize',14)
hold off
% normdatatable = [runs,d,ep_normd,lerrd,uerrd,ep_predd,ep_pred_diff,ep_pred_disl];
% preddatatable = [d3 epdex epf eps epg];
end

if min(T)<max(T)
epf = zeros(length(T3),1);
eps = zeros(length(T3),1);
epg = zeros(length(T3),1);
% normalize at reference P0, sig0, d0
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        ep_pred_diff = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T)).*(Cw0.^mr1);
        epf = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T3)).*(Cw0.^mr1);
    elseif strcmp(mech,'ws')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T)).*(Cw0.^mr3);
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T3)).*(Cw0.^mr3);
    elseif strcmp(mech,'wg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T)).*(Cw0.^mr5);
        epg = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T3)).*(Cw0.^mr5);
    end
end
ep_predT = ep_pred_disl + ep_pred_diff + ep_pred_gbs;
epTex = epf + eps + epg;
% normalize observed strain rate
ep_normT = (edotX./edot_pred).*ep_predT;
% % calculate error 
errT = err;
uerrT = exp(log(ep_normT)+errT)-ep_normT;
lerrT = ep_normT-exp(log(ep_normT)-errT);
if ci==1
% calculate confidence intervals
[epsL, epsU] = construct_CI(90,[T3 P0*ones(size(T3)) sig0*ones(size(T3)) d0*ones(size(T3)) Cw0*ones(size(T3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T3 P0*ones(size(T3)) sig0*ones(size(T3)) d0*ones(size(T3)) Cw0*ones(size(T3))],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,3,3)
hold on
if ci==1
    fill_area(1e3./T3,epsL,epsU,grey,0.3,'90% CI')
    fill_area(1e3./T3,epsL2,epsU2,grey2,0.3,'50% CI')
end
if Xn == 0 && XXn==0
    errorbar(1e3./T,ep_normT,lerrT,uerrT,'ko','DisplayName','data')
else %if XXn>0
    for xx=1:XXn
        hh=errorbar(1e3./T(runsx==xx),ep_normT(runsx==xx),lerrT(runsx==xx),uerrT(runsx==xx),'o','DisplayName',['run' num2str(xx)]);
        %set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(1e3./T3,epTex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.wf.active, plot(1e3./T3,epf,'r--','DisplayName','diffusion'); end
    if model.ws.active, plot(1e3./T3,eps,'g--','DisplayName','dislocation'); end
    if model.wg.active, plot(1e3./T3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log')
xlabel('1000/T [K^{-1}]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable T, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(d0) '\mu m, ' ...
    num2str(sig0) 'MPa' num2str(Cw0) 'ppm H/Si']},'FontSize',14)
hold off
% normdatatable = [runs,T,ep_normT,lerrT,uerrT,ep_predT,ep_pred_diff,ep_pred_disl];
% preddatatable = [T3 epTex epf eps epg];
end

if min(P)<max(P)
epf = zeros(length(P3),1);
eps = zeros(length(P3),1);
epg = zeros(length(P3),1);
% normalize at reference T0, sig0, d0
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        ep_pred_diff = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
        epf = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P3*mV1)*1e3./(R*T0)).*(Cw0.^mr1);
    elseif strcmp(mech,'ws')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P*mV3)*1e3./(R*T0)).*(Cw0.^mr3);
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P3*mV3)*1e3./(R*T0)).*(Cw0.^mr3);
    elseif strcmp(mech,'wg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
        epg = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P3*mV5)*1e3./(R*T0)).*(Cw0.^mr5);
    end
end
ep_predP = ep_pred_disl + ep_pred_diff + ep_pred_gbs;
epPex = epf + eps + epg;
% normalize observed strain rate
ep_normP = (edotX./edot_pred).*ep_predP;
% calculate error 
errP = err;
uerrP = exp(log(ep_normP)+errP)-ep_normP;
lerrP = ep_normP-exp(log(ep_normP)-errP);
if ci==1
% calculate confidence intervals
[epsL, epsU] = construct_CI(90,[T0*ones(size(P3)) P3 sig0*ones(size(P3)) d0*ones(size(P3)) Cw0*ones(size(P3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(P3)) P3 sig0*ones(size(P3)) d0*ones(size(P3)) Cw0*ones(size(P3))],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,3,4)
hold on
if ci==1
    fill_area(P3,epsL,epsU,grey,0.3,'90% CI')
    fill_area(P3,epsL2,epsU2,grey2,0.3,'50% CI')
end
if XXn==0
    errorbar(P,ep_normP,lerrP,uerrP,'ko','DisplayName','data')
else %if XXn>0
    for xx=1:XXn
        hh=errorbar(P(runsx==xx),ep_normP(runsx==xx),lerrP(runsx==xx),uerrP(runsx==xx),'o','DisplayName',['run' num2str(xx)]);
        %set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(P3,epPex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.wf.active, plot(P3,epf,'r--','DisplayName','diffusion'); end
    if model.ws.active, plot(P3,eps,'g--','DisplayName','dislocation'); end
    if model.wg.active, plot(P3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log')
xlabel('P [GPa]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend 
title({['Inversion at variable P, Xn=' num2str(Xn)],...
    [num2str(T0) 'K, ' num2str(d0) '\mu m, ' ...
    num2str(sig0) 'MPa' num2str(Cw0) 'ppm H/Si']},'FontSize',14)
hold off
% normdatatable = [runs,P,ep_normP,lerrP,uerrP,ep_predP,ep_pred_diff,ep_pred_disl];
% preddatatable = [P3 epPex epf eps epg];
end


if min(Cw)<max(Cw)
epf = zeros(length(Cw3),1);
eps = zeros(length(Cw3),1);
epg = zeros(length(Cw3),1);
% normalize at reference T0, sig0, d0, P0
for i=1:Nmech
    mech = flows{i};
    if strcmp(mech,'wf')
        ep_pred_diff = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw.^mr1);
        epf = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0)).*(Cw3.^mr1);
    elseif strcmp(mech,'ws')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*(Cw.^mr3);
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*(Cw3.^mr3);
    elseif strcmp(mech,'wg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw.^mr5);
        epg = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0)).*(Cw3.^mr5);
    end
end
ep_predw = ep_pred_disl + ep_pred_diff + ep_pred_gbs;
epwex = epf + eps + epg;
% normalize observed strain rate
ep_normw = (edotX./edot_pred).*ep_predw;
% calculate error 
errw = err;
uerrw = exp(log(ep_normw)+errw)-ep_normw;
lerrw = ep_normw-exp(log(ep_normw)-errw);
if ci==1
% calculate confidence intervals
[epsL, epsU] = construct_CI(90,[T0*ones(size(P3)) P0*ones(size(Cw3)) sig0*ones(size(P3)) d0*ones(size(P3)) Cw3],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(P3)) P0*ones(size(Cw3)) sig0*ones(size(P3)) d0*ones(size(P3)) Cw3],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,3,5)
hold on
if ci==1
    fill_area(Cw3,epsL,epsU,grey,0.3,'90% CI')
    fill_area(Cw3,epsL2,epsU2,grey2,0.3,'50% CI')
end
if XXn==0
    errorbar(Cw,ep_normw,lerrw,uerrw,'ko','DisplayName','data')
else %if XXn>0
    for xx=1:XXn
        hh=errorbar(Cw(runsx==xx),ep_normw(runsx==xx),lerrw(runsx==xx),uerrw(runsx==xx),'o','DisplayName',['run' num2str(xx)]);
        %set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(Cw3,epwex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.wf.active, plot(Cw3,epf,'r--','DisplayName','diffusion'); end
    if model.ws.active, plot(Cw3,eps,'g--','DisplayName','dislocation'); end
    if model.wg.active, plot(Cw3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log')
xlabel('C_w [ppm H/Si]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend 
title({['Inversion at variable P, Xn=' num2str(Xn)],...
    [num2str(T0) 'K, ' num2str(P0) 'GPa, ' num2str(d0) '\mu m, ' ...
    num2str(sig0) 'MPa']},'FontSize',14)
hold off
% normdatatable = [runs,P,ep_normP,lerrP,uerrP,ep_predP,ep_pred_diff,ep_pred_disl];
% preddatatable = [P3 epPex epf eps epg];
end