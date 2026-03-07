% % calculate autocorrelation function for seed=0 run 
% for dry (or wet) MCMC inversion output
% 
% % format of MCMC output file :
% % iteration#, simple_chi2, full_chi2, p1, E1 (J/mol), V1 (m3/mol), ...
% % n3 E3 V3 A1 A3

% path = (string) path to MCMC output file
% para = (cell) vector of parameters names inverted in MCMC
% kn = number of unknown flow-law parameters determined 
% by MCMC (ignoring the scaling coeff. & inter-run biases)
% fname = name of the MCMC output file run with seed=0
% tstr = plot title string

clear
clc

path = '../CaseS1/'; 
para = [{'p5'};{'n5'}];
kn=2;
fname = ['out.S1g.0'];
tstr = 'S1g';

path = '../CaseS1/'; 
para = [{'p1'}];
kn=1;
fname = ['out.S1f.0'];
tstr = 'S1f';

path = '../CaseS2/'; 
para = [{'p1'};{'E1'};{'V1'}];
kn=3;
fname = ['out.S2f.0'];
tstr = 'S2f';

path = '../CaseS3/'; 
para = [{'p1'}];
kn=1;
fname = ['out.S3f.0'];
tstr = 'S3f';

path = '../CaseS3/'; 
para = [{'p1'}];
kn=1;
fname = ['out.S3fX.0'];
tstr = 'S3fX';

path = '../CaseS4/'; 
para = [{'p1'};{'E1'}];
kn=2;
fname = ['out.S4f.noX.0'];
tstr = 'S4f (no X)';

path = '../CaseS4/'; 
para = [{'p1'};{'E1'}];
kn=2;
fname = ['out.S4fX.0'];
tstr = 'S4fX';

path = '../CaseS4/'; 
para = [{'p1'};{'E1'}];
kn=2;
fname = ['out.S4fX-wb.0'];
tstr = 'S4fX (wb)';

path = '../CaseS5/'; 
para = [{'p5'};{'n5'}];
kn=2;
fname = ['out.S5g.0'];
tstr = 'S5g';

path = '../CaseS5/'; 
para = [{'n3'}];
kn=1;
fname = ['out.S5s.0'];
tstr = 'S5s';

path = '../CaseS6/'; 
para = [{'n3'};{'E3'};{'V3'}];
kn=3;
fname = ['out.S6s.0'];
tstr = 'S6s';

path = '../CaseS7/'; 
para = [{'p1'}];
kn=1;
fname = ['out.S7f.0'];
tstr = 'S7f';

path = '../CaseS7/'; 
para = [{'n3'}];
kn=1;
fname = ['out.S7s.0'];
tstr = 'S7s';

path = '../CaseS7/'; 
para = [{'p5'};{'n5'}];
kn=2;
fname = ['out.S7g.0'];
tstr = 'S7g';

path = '../CaseS7/'; 
para = [{'p1'};{'n3'}];
kn=2;
fname = ['out.S7fs.0'];
tstr = 'S7fs';

path = '../CaseS8/'; 
para = [{'p1'};{'E1'};{'V1'};{'n3'};{'E3'};{'V3'}];
kn=6;
%     fname = ['out.S8fsX.0']; % w/o sticky parameters
%     tstr = 'S8fsX';
fname = ['out.S8fsX.sticky.0'];
tstr = 'S8fsX (sticky)'; % w/ sticky parameters 
%  does not converge too well. lag ~ 8e3. 
                                 
path = '../CaseS9/'; 
para = [{'r4'};{'n4'};{'E4'};{'V4'};{'p2'};{'r2'};{'E2'};{'V2'}];
kn=8;
fname = ['out.S9fsX.0'];
tstr = 'S9fsX\_wet'; % resample at 6000.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 10e3; %2e3; <-- no. of lags
%
outm = load([path fname '.dat']);
outm(:,[1:3])=[];

ACF = zeros(M,kn+1);
for i=1:kn
    mean1 = mean(outm(:,i));
    var1 = (outm(:,i)-mean1)'*(outm(:,i)-mean1);
    for j=0:M
        ACF(j+1,i+1)=sum(((outm(j+1:end,i)-mean1).*...
            (outm(1:end-j,i)-mean1))/var1);
        if i==1
            ACF(j+1,i) = j;
        end
    end
end

figure
hold on
for kk=1:kn
    subplot(ceil(kn/3),3,kk)
    plot(ACF(:,1),ACF(:,kk+1),'b')
    legend(para{kk})
    xlabel('lags','FontSize',14)
    ylabel('ACF','FontSize',14)
    if kk==1
        title(['ACF of inversion (' tstr ')'],'FontSize',14)
    end
end
hold off

% % write ACF output
% fn=fopen([path 'ACF.' fname '.dat'],'w');
% fspec = '%d';
% for i=2:length(ACF(1,:))
%   fspec = [fspec ' %8.6f'];
% end
% fspec = [fspec '\n'];
% fprintf(fn,fspec,ACF');
% fclose(fn);
