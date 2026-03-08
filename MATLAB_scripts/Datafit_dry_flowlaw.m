% Analyze data S1-S8 for only a single & composite flow law (dry)
% - gbs, diffusion, or dislocation.
% Include inter-run bias where assumed.
% 1. Resample output files at intervals determined from ACF
% 2. Combine resampled parallel runs
% 3. Determine scaling coefficients from MCMC output
% 4. Calculate mean, std. dev., correlations. Write in a file
% 5. Draw pdf of each unknown parameter.
% 6. Draw normalized strain rate graphs to test data fit
% GBS creep - invert for n5, p5, E5, V5.
% Diffusion creep - invert for p1, E1, V1. 
% Dislocation creep - invert for n3, E3, V3.

% Provide the following information- 
% 1. path = path to input & output data files (e.g., './')
% 2. dataf = name of input data file <filename.dat>
% 3. outf = name of MCMC output file (without '.dat')
% 4. flows = array of the identifiers for each of the   
%     flow laws assumed, provided  in the same order 
%     in which they were declared in the MCMC file: 
%       'df' - dry diffusion
%       'ds' - dry dislocation
%       'dg' - dry GBS
%     E.g., If parameter file was declared with -Pa<> -Pf<>,
%           then An = 2 and flow = {'df', 'dg'};
% 5. An = the number of scaling coefficients inverted for
% 6. kn = total number of flow-law parameters inverted for 
%     (excluding the scaling co-efficients).
% 7. invertparams = An array of parameters being inverted for, in the same
%     sequence in which they appear in the MCMC output file
%     {'df.p','df.E','df.V','df.A'} - diffusion creep
%     {'ds.n','ds.E','ds.V','ds.A'} - dry dislocation creep
%     {'dg.p','dg.n','dg.E','dg.V','dg.A'} - dry gbs creep
%          E.g., If parameter file was declared with :
%               -Pa0/0/1/3/10/10/0/0 -Pf<0/0/1/3/1/5/0/0/0/0,
%          then flow = {'df', 'dg'}; An = 2; kn = 3; invertparams = {'df.p', 'dg.p', 'dg.n', 'df.A', 'dg.A'};
% 8. nfixed = no. of fixed flow-law parameters. E.g., nfixed=4 in the prev example
% 9. fixed.(<mech>).<para> = fixed_val : 
%     parameter para of mechanism mech whose value 
%     is assumed to be contant = fixed_val.
%     Call individually for each such parameter (para = p,n,E,V,A) 
%     for each mechanism (mech = df,ds,dg) if nfixed>1
%           E.g., nfixed = 4; fixed.df.E=10; fixed.df.V=10; fixed.dg.E=0; fixed.dg.V=0;
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
printresults = 1; % write means etc. or not

idir = 'input/';
odir = 'output/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % choose a Case and associated input parameters
% % Uncomment only that case and comment the others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% % (Case S1g)
% path = '../CaseS1/';
% dataf = 'dataS1.noX.dat'; 
% outf = 'out.S1g.'; % 
% flows = {'dg'}; 
% An = 1;
% kn = 2; % p5,n5
% invertparams ={'dg.p','dg.n','dg.A'}; % list of unknown params
% nfixed = 2;
% fixed.dg.E = 0;
% fixed.dg.V = 0;
% XXn = 3;
% Xn=0;
% dataf2 = 'dataS1.dat'; % original data grouping
% parallelruns = 3; 
% resample = 100;
% % if 'expected values' of unknown parameters is available, declare it here
% % in the same seq. in which the unknown parameters were declared in
% % invertparams
% para0 = [2 5e5*exp(-(300+0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% % RESULTS:
% % % means = [1.8805 0.9907 -4.6963] +/- [0.108 0.095 0.0989]
% % % 2s.d. includes the true values of the parameters. 24 data points.
% % % % Expected A1=5e5*exp(-(300+0.1*10)*1e3/(R*1523))=2.3748e-05; log(A1)=-4.6244 

% % Use same data set, and invert for diffusion creep. (Case S1f)
% path = '../CaseS1/';
% dataf = 'dataS1.noX.dat'; 
% outf = 'out.S1f.';
% flows = {'df'}; 
% An = 1;
% kn = 1; %p1
% invertparams ={'df.p','df.A'}; % list of unknown params
% nfixed = 2;
% fixed.df.E = 0;
% fixed.df.V = 0;
% XXn=3;
% Xn = 0;
% dataf2 = 'dataS1.dat'; % original data grouping
% parallelruns = 3; % done
% resample = 100;
% % if 'expected values' of unknown parameters is available, declare it here
% % in the same seq. in which the unknown parameters were declared in
% % invertparams
% para0 = [2 5e5*exp(-(300+0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% % RESULTS: 
% % means = [1.8892 -4.6990] +/- [0.0879 0.0871]
% % As expected, p1 is within 2s.d. of mean. 
% % Expected A1=5e5*exp(-(300+0.1*10)*1e3/(R*1523))=2.3748e-05; log(A1)=-4.6244

% % % Data are constructed ensuring that they are in the diffusion creep regime
% % Invert for diff. (Case S2f)
% path = '../CaseS2/';
% dataf = 'dataS2.noX.dat'; 
% outf = 'out.S2f.'; 
% flows = {'df'}; 
% An = 1;
% kn = 3; %p1,E1,V1
% invertparams ={'df.p','df.E','df.V','df.A'}; % list of unknown params
% nfixed = 0;
% XXn=8;
% Xn = 0;
% dataf2 = 'dataS2.dat'; % original data grouping
% parallelruns = 3; % done
% resample = 100;
% % if 'expected values' of unknown parameters is available, declare it here
% % in the same seq. in which the unknown parameters were declared in
% % invertparams
% para0 = [2 300 10 5E5];
% % RESULTS:
% % p1 is outside 1s.d. and at 2sd. rest are ok.
% % means = [1.948 300.950 9.884 5.6975]; std=[0.0447 5.8362 0.4667 0.0437]
% % 3 runs: [1.9477 300.9784 9.8834 5.6983]; sd=[0.0441 5.7910 0.4647 0.0431]

% % Modify data set S1 by adding inter-run bias (Case S3)
% % Invert data set for diffusion creep but without considering
% % X (Case S3f)
% path = '../CaseS3/';
% dataf = 'dataS3.noX.dat'; 
% outf = 'out.S3f.'; 
% flows = {'df'}; 
% An = 1;
% kn = 1; %p1
% invertparams ={'df.p','df.A'}; % list of unknown params
% nfixed = 2;
% fixed.df.E = 0;
% fixed.df.V = 0;
% Xn = 0; % we assumed no inter-run bias
% XXn = 3; % original data set has 3 inter-run biases or data-groups
% dataf2 = 'dataS3.dat'; % original data grouping
% parallelruns = 1; % done
% resample = 100;
% para0 = [2 5e5*exp(-(300+0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% % RESULTS:
% % means = [2.2282 -4.3853] +/- [0.1006 0.0984]
% % Whereas p1 is just beyond 2s.d. of mean, by not considering 
% % inter-run biases, we get v large misfits.
% % Expected A1=5e5*exp(-(300+0.1*10)*1e3/(R*1523))=2.3748e-05; log(A1)=-4.6244

% % Invert dataset S3 for diffusion creep with X (Case S3fX).
% path = '../CaseS3/';
% dataf = 'dataS3.dat'; 
% outf = 'out.S3fX.'; 
% flows = {'df'}; 
% An = 1;
% kn = 1; %p1
% invertparams ={'df.p','df.A'}; % list of unknown params
% nfixed = 2;
% fixed.df.E = 0;
% fixed.df.V = 0;
% Xn = 3; % we assumed no inter-run bias
% XXn = 3; % original data set has 3 inter-run biases or data-groups
% parallelruns = 1; % done
% resample = 100;
% para0 = [2 5e5*exp(-(300+0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% paraX0 = [-1 0 1]; %%[-1 -0.5 0 0.5 1];
% % RESULTS: 
% % means = [1.879 -1.0125 -0.0149 1.0275 -4.7084] +/- [0.0959 0.0489 0.0503 0.0488 0.0947]
% % Whereas p1 is within 2s.d. of mean, the inter-run biases are also recovered.
% % Expected A1=5e5*exp(-(300+0.1*10)*1e3/(R*1523))=2.3748e-05; log(A1)=-4.6244

% % Modify data set S1 such that each run is at a different temperature. 
% % Invert for diffusion creep p1, E1 without considering X. (Case S4f)
% path = '../CaseS4/';
% dataf = 'dataS4.noX.dat'; 
% outf = 'out.S4f.'; 
% flows = {'df'}; 
% An = 1;
% kn = 2; %p1,E1
% invertparams ={'df.p','df.E','df.A'}; % list of unknown params
% nfixed = 1;
% fixed.df.V = 0;
% Xn = 0; % we assumed no inter-run bias
% XXn = 3; % original data set has 3 inter-run biases or data-groups
% parallelruns = 1; % done
% resample = 100;
% dataf2 = 'dataS4.dat'; 
% para0 = [2 300 5e5*exp(-(0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% % RESULTS:
% % means = [2.0097 406.38 9.37] +/- [0.0997 8.488 0.0938]
% % Both p1 and E1 are constrained v rightly, and p1 is fully recovered.
% % However, E1 is higher than expected, and data show large misfits. 

% % Invert dataset S4 for diffusion creep p1, E1 considering X. (Case S4fX-nb)
% path = '../CaseS4/';
% dataf = 'dataS4.dat'; 
% outf = 'out.S4fX.'; 
% flows = {'df'}; 
% An = 1;
% kn = 2; %p1,E1
% invertparams ={'df.p','df.E','df.A'}; % list of unknown params
% nfixed = 1;
% fixed.df.V = 0;
% Xn = 3; % 
% XXn = 3; % original data set has 3 inter-run biases or data-groups
% parallelruns = 1; % done
% resample = 100;
% para0 = [2 300 5e5*exp(-(0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% paraX0 = [-1 0 1]; %%[-1 -0.5 0 0.5 1];
% % RESULTS:
% % % means = [1.807 259.83 -1.043 -0.255 1.298 4.148] +/- [0.0873 321.95 0.0875 1.714 1.641  0.088]
% % % Whereas p1 is just beyond 2s.d. of mean, E1 is not constrained. 1s.d. > mean &
% % % we see a large peak at E1=0. a priori bounds were 0<E1<1000. Widen
% % % bounds.

% % Repeat above inversion with even wider bounds on E1.(Case S4fX-wb)
% path = '../CaseS4/';
% dataf = 'dataS4.dat'; 
% outf = 'out.S4fX-wb.'; 
% flows = {'df'}; 
% An = 1;
% kn = 2; %p1,E1
% invertparams ={'df.p','df.E','df.A'}; % list of unknown params
% nfixed = 1;
% fixed.df.V = 0;
% Xn = 3; % 
% XXn = 3; % original data set has 3 inter-run biases or data-groups
% parallelruns = 1; % done
% resample = 100;
% para0 = [2 300 5e5*exp(-(0.1*10)*1e3/(R*1523))]; %[2 300 10 5E5];
% paraX0 = [-1 0 1]; %%[-1 -0.5 0 0.5 1];
% % In the paper's figures, I mistakenly removed the MCMC 
% % output with chi2full>20, so we don't see a second peak
% % at E=-1000 (lower a priori bound), which indicates that 
% % data prefer still lower values of E. Means are slightly different too. 
% % RESULTS: 
% % means = [1.8088 178.47 -1.0235 -0.688 1.7116 1.3509] +/- [0.0864 559.08 0.1354 2.9757 2.8496 0.0975]
% % Whereas p1 is just beyond 2s.d. of mean, E1 is poorly constrained. 1s.d. > mean &
% % we see a large peak at E1~0. a priori bounds were -1000<E1<2000. Main peak is
% % within the a priori bounds. 

% % Invert data S5 for p3 & n3 for dislocation creep using
% % the flow law for gbs creep (i.e., w/o assuming p3=0).
% % Therefore, P=0.4 GPa here. E_i=V_i=0 
% % (Case S5g).
% path = '../CaseS5/';
% dataf = 'dataS5.dat';
% outf = 'out.S5g.';
% flows = {'dg'}; 
% An = 1;
% kn = 2; %p5,n5
% invertparams ={'dg.p','dg.n','dg.A'}; % list of unknown params
% nfixed = 2;
% fixed.dg.E = 0;
% fixed.dg.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 3; 
% resample = 200;
% para0 = [0 3 2.1e5*exp(-(465+0.1*15)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % As expected, p3=0. So, no diffusion creep. 
% % 3 runs: means=[0.0062 2.8547 -10.4636]; sd=[0.1753 0.2261 0.3141]
% % Expected A3=2.1e5*exp(-(465+0.4*15)*1e3/(R*1523))=1.4738e-11; log10(A3)=-10.8316

% % Invert the data S5 for n3 for dislocation creep using
% % the flow law for dislocation creep (i.e., assuming p1=0).  (Case S5s)
% path = '../CaseS5/';
% dataf = 'dataS5.dat'; 
% outf = 'out.S5s.'; 
% flows = {'ds'}; 
% An = 1;
% kn = 1; %n3
% invertparams ={'ds.n','ds.A'}; % list of unknown params
% nfixed = 2;
% fixed.ds.E = 0;
% fixed.ds.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 3; 
% resample = 100;
% para0 = [0 3 2.1e5*exp(-(465+0.1*15)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % expected values lie at 2s.d. of mean.
% % 3 runs: means=[2.8465 -10.4528]; sd=[0.0679 0.1633]
% % Expected A3=2.1e5*exp(-(465+0.4*15)*1e3/(R*1523))=1.4738e-11; log10(A3)=-10.8316

% % Invert dataset S6 for disl. without fixing any params.(Case S6s)
% path = '../CaseS6/';
% dataf = 'dataS6.noX.dat'; 
% outf = 'out.S6s.'; 
% flows = {'ds'}; 
% An = 1;
% kn = 3; %n3,E3,V3
% invertparams ={'ds.n','ds.E','ds.V','ds.A'}; % list of unknown params
% nfixed = 0;
% Xn = 0; % 
% XXn = 8; % 
% dataf2 = 'dataS6.dat'; 
% parallelruns = 3; 
% resample = 100;
% para0 = [0 3 2.1e5*exp(-(465+0.1*15)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % true values lie within 2s.d. of means.
% % mean = [2.9573 460.6663 14.8674 5.2802]; sd = [0.0364 2.2107 0.1644 0.0829]

% % Use data set S7 over a range of d & sigma, constructed 
% % for composite rheology to invert for diffusion creep. (Case S7f)
% path = '../CaseS7/';
% dataf = 'dataS7.dat'; 
% outf = 'out.S7f.'; 
% flows = {'df'}; 
% An = 1;
% kn = 1; %p1
% invertparams ={'df.p','df.A'}; % list of unknown params
% nfixed = 2;
% fixed.df.E = 0;
% fixed.df.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 1; 
% resample = 100;
% para0 = [2 5e5*exp(-(300+0.3*10)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % means = [0.6487 -5.8891] +/- [0.0681 0.0808]

% % Use data set S7 over a range of d & sigma, constructed 
% % for composite rheology to invert for dislocation creep. (Case S7s)
% path = '../CaseS7/';
% dataf = 'dataS7.dat'; 
% outf = 'out.S7s.'; 
% flows = {'ds'}; 
% An = 1;
% kn = 1; %n3
% invertparams ={'ds.n','ds.A'}; % list of unknown params
% nfixed = 2;
% fixed.ds.E = 0;
% fixed.ds.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 1; 
% resample = 100;
% para0 = [3 2.1E5*exp(-(465+0.3*15)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % means = [0.7644 -6.2450] +/- [0.0380 0.0642]

% % Use data set S7 over a range of d & sigma, constructed 
% % for composite rheology toinvert for gbs creep. (Case S7g)
% path = '../CaseS7/';
% dataf = 'dataS7.dat'; 
% outf = 'out.S7g.'; 
% flows = {'dg'}; 
% An = 1;
% kn = 2; %p5,n5
% invertparams ={'dg.p','dg.n','dg.A'}; % list of unknown params
% nfixed = 2;
% fixed.dg.E = 0;
% fixed.dg.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 1; 
% resample = 200;
% para0 = [2 3]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % means = [2.9487 2.2871 -5.3667] +/- [0.0927 0.0724 0.1217]

% % Use data set S7 over a range of d & sigma, constructed 
% % for composite rheology to invert for diff+disl creep. (Case S7fs)
% % Given time constraints, I only ran till iteration# = 500,000
% path = '../CaseS7/';
% dataf = 'dataS7.dat'; 
% outf = 'out.S7fs.'; 
% flows = {'df','ds'}; 
% An = 2;
% kn = 2; %p1,n3
% invertparams ={'df.p','ds.n','df.A','ds.A'}; % list of unknown params
% nfixed = 4;
% fixed.df.E = 0;
% fixed.df.V = 0;
% fixed.ds.E = 0;
% fixed.ds.V = 0;
% Xn = 0; % 
% XXn = 0; % 
% parallelruns = 3; 
% resample = 100;
% para0 = [2 3 5e5*exp(-(300+0.3*10)*1e3/(R*1523)) 2.1e5*exp(-(465+0.3*15)*1e3/(R*1523))]; %[2 300 10 3 465 15 5E5 2.1e5];
% % RESULTS: 
% % Thankfully, excellent fit!
% % 1 run: means=[1.9975 3.0309 -4.7019 -10.8442]+/-[0.1622 0.2480 0.1496 0.5300]
% % all runs: [1.9932 3.0345 -4.7061 -10.8518]+/-[0.1585 0.2453 0.1466 0.5235]
% % % Expected: [2 3 5e5*exp(-(300+0.3*10)*1e3/(8.314*1523))
% % 2.1e5*exp(-(465+0.3*15)*1e3/(8.314*1523))]
% % = [2 3 -4.6935 -10.7809]

% Case S8fsX
path = '../CaseS8/';
dataf = 'dataS8.dat'; 
% outf = 'out.S8fsX.'; # no sticky parameters
outf = 'out.S8fsX.sticky.'; %assuming sticky -G1/3 -G2/5 -G3/6, resample=5000-7000
flows = {'df','ds'}; 
An = 2;
kn = 6; %p1,E1,V1,n3,E3,V3
invertparams ={'df.p','df.E','df.V','ds.n','ds.E','ds.V','df.A','ds.A'}; % list of unknown params
nfixed = 0;
Xn = 5; % 
XXn = 5; % 
parallelruns = 3; 
resample = 8000; % 5000
para0 = [2 300 10 3 465 15 5E5 2.1e5];
paraX0 = [-1 -0.5 0 0.5 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Main script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dry or wet?
if strcmp(flows{1}(1),'d')
    condition = 'dry';
elseif strcmp(flows{1}(1),'w')
    %condition = 'wet';
    error('This script is only for dry flow laws.');
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
outm = combine_data([path odir],files,resample); % resampled & compiled data
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
    fn=fopen([path odir outf 'resampled.dat'],'w');
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
% outm=load([path odir outf 'resampled.dat']);
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
obsdata = load([path idir dataf]);
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
    obsdatax = load([path idir dataf2]);
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
        plothist(outm(:,2),nbins,'simple \chi^2_J/N_p',0,[3,3,1]);
    elseif i==0
        plothist(chi2M_vec,nbins,'\chi^2_M/N',0,[3,3,2]);
    else 
        plothist(outm(:,3+i),nbins,paramstr{i},1,[3,3,i+2]);
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
        plothist(outm(:,3+kn+i),nbins,paramstr{kn+i},1,[3,3,i]); %['X_' num2str(i)],1,[3,3,i]);
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
if model.df.active
    mp1 = model.df.means(1); mE1 = model.df.means(2);
    mV1 = model.df.means(3); mA1 = model.df.means(4);
end
if model.ds.active
    mn3 = model.ds.means(1); mE3 = model.ds.means(2);
    mV3 = model.ds.means(3); mA3 = model.ds.means(4);
end
if model.dg.active
    mp5 = model.dg.means(1); mn5 = model.dg.means(2); 
    mE5 = model.dg.means(2);
    mV5 = model.dg.means(3); mA5= model.dg.means(5);
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
    if strcmp(mech,'df')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:4),Xmat,obsdata,mech);
        edot_pred_diff = ep_pred(:,i);
    end
    if strcmp(mech,'ds')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:4),Xmat,obsdata,mech);
        edot_pred_disl = ep_pred(:,i);
    end
    if strcmp(mech,'dg')
        [relvar(:,i),rvarej,ep_pred(:,i)] = rel_var(model.(mech).means(1:5),Xmat,obsdata,mech);
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
    % sig0 = 50; %100; % MPa
    sig0 = 20;
    sig3 = linspace(min(sig)/2,max(sig)*2,100)';
end
if min(d)==max(d)
    d0 = min(d);
else
    %d0 = 12; % micron
    d0 = 10; % microns
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
%     P0 = 0.1; %1; % GPa
%     P0 = 0.3; %1; % GPa
    P0 = 0.2; % GPa
    P3 = linspace(0.1,max(P)+1,100)'; 
end
T0 = 1523; % K
sig0 = 100; % MPa
d0 = 15; %12; % micron
P0 = 0.3; %1; % GPa

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
    if strcmp(mech,'df')
        ep_pred_diff = mA1*(sig).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
        epf = mA1*(sig3).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
    elseif strcmp(mech,'ds')
        ep_pred_disl = mA3*(sig.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
        eps = mA3*(sig3.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0));
    elseif strcmp(mech,'dg')
        ep_pred_gbs = mA5*(sig.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0));
        epg = mA5*(sig3.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0));
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
[epsL, epsU] = construct_CI(90,[T0*ones(size(sig3)) P0*ones(size(sig3)) sig3 d0*ones(size(sig3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(sig3)) P0*ones(size(sig3)) sig3 d0*ones(size(sig3))],flows,model,Nout); % 50% conf. intervals
end
subplot(2,2,1)
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
        set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh=plot(sig3,epsex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.df.active, plot(sig3,epf,'r--','DisplayName','diffusion'); end
    if model.ds.active, plot(sig3,eps,'g--','DisplayName','dislocation'); end
    if model.dg.active, plot(sig3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log','XScale','log')
xlabel('\sigma [MPa]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable stress, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(d0) '\mu m, ' ...
    num2str(T0) 'K']},'FontSize',14)
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
    if strcmp(mech,'df')
        ep_pred_diff = mA1*(sig0).*(d.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
        epf = mA1*(sig0).*(d3.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T0));
    elseif strcmp(mech,'ds')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*ones(size(d));
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T0)).*ones(size(d3));
    elseif strcmp(mech,'dg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0));
        epg = mA5*(sig0.^mn5).*(d3.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T0));
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
[epsL, epsU] = construct_CI(90,[T0*ones(size(d3)) P0*ones(size(d3)) sig0*ones(size(d3)) d3],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(d3)) P0*ones(size(d3)) sig0*ones(size(d3)) d3],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,2,2)
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
        set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(d3,epdex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.df.active, plot(d3,epf,'r--','DisplayName','diffusion'); end
    if model.ds.active, plot(d3,eps,'g--','DisplayName','dislocation'); end
    if model.dg.active, plot(d3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log','XScale','log')
xlabel('d [\mum]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable grain size, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(T0) 'K, ' ...
    num2str(sig0) 'MPa']},'FontSize',14)
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
    if strcmp(mech,'df')
        ep_pred_diff = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T));
        epf = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P0*mV1)*1e3./(R*T3));
    elseif strcmp(mech,'ds')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T));
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P0*mV3)*1e3./(R*T3));
    elseif strcmp(mech,'dg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T));
        epg = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P0*mV5)*1e3./(R*T3));
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
[epsL, epsU] = construct_CI(90,[T3 P0*ones(size(T3)) sig0*ones(size(T3)) d0*ones(size(T3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T3 P0*ones(size(T3)) sig0*ones(size(T3)) d0*ones(size(T3))],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,2,3)
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
        set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(1e3./T3,epTex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.df.active, plot(1e3./T3,epf,'r--','DisplayName','diffusion'); end
    if model.ds.active, plot(1e3./T3,eps,'g--','DisplayName','dislocation'); end
    if model.dg.active, plot(1e3./T3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log')
xlabel('1000/T [K^{-1}]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend
title({['Inversion at variable T, Xn=' num2str(Xn)],...
    [num2str(P0) 'GPa, ' num2str(d0) '\mu m, ' ...
    num2str(sig0) 'MPa']},'FontSize',14)
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
    if strcmp(mech,'df')
        ep_pred_diff = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P*mV1)*1e3./(R*T0));
        epf = mA1*(sig0).*(d0.^(-mp1)).*exp(-(mE1+P3*mV1)*1e3./(R*T0));
    elseif strcmp(mech,'ds')
        ep_pred_disl = mA3*(sig0.^mn3).*exp(-(mE3+P*mV3)*1e3./(R*T0));
        eps = mA3*(sig0.^mn3).*exp(-(mE3+P3*mV3)*1e3./(R*T0));
    elseif strcmp(mech,'dg')
        ep_pred_gbs = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P*mV5)*1e3./(R*T0));
        epg = mA5*(sig0.^mn5).*(d0.^(-mp5)).*exp(-(mE5+P3*mV5)*1e3./(R*T0));
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
[epsL, epsU] = construct_CI(90,[T0*ones(size(P3)) P3 sig0*ones(size(P3)) d0*ones(size(P3))],flows,model,Nout); % 90% conf. intervals
[epsL2, epsU2] = construct_CI(50,[T0*ones(size(P3)) P3 sig0*ones(size(P3)) d0*ones(size(P3))],flows,model,Nout); % 50% conf. intervals
end
figure(hfig)
subplot(2,2,4)
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
        set(hh,'MarkerEdgeColor',cols{xx})
    end
end
hh = plot(P3,epPex,'b','DisplayName','prediction');
set(hh,'LineWidth',1)
if Nmech>1
    if model.df.active, plot(P3,epf,'r--','DisplayName','diffusion'); end
    if model.ds.active, plot(P3,eps,'g--','DisplayName','dislocation'); end
    if model.dg.active, plot(P3,epg,'m--','DisplayName','dislocation'); end
end
set(gca,'YScale','log')
xlabel('P [GPa]','FontSize',14)
ylabel('\epsilon''','FontSize',14)
legend 
title({['Inversion at variable P, Xn=' num2str(Xn)],...
    [num2str(T0) 'K, ' num2str(d0) '\mu m, ' ...
    num2str(sig0) 'MPa']},'FontSize',14)
hold off
% normdatatable = [runs,P,ep_normP,lerrP,uerrP,ep_predP,ep_pred_diff,ep_pred_disl];
% preddatatable = [P3 epPex epf eps epg];
end
