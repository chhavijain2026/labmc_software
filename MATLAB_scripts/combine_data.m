% combine different parallel runs after resampling them
% input : 
% files = list of names of all .dat/.mat files of parallel runs;
% resample = resampling interval according to ACF

function data = combine_data(path,file,resample)

n = length(file); % no. of parallel runs being combined

k=0;
for i=1:n
    D=load([path file{i,:}]);
    if i==1
        % for 0th run, output is printing at every step
        % We want to remove the first 10^4 outputs
        remrows = 10^4;
        resample2 = round(resample/100)*100;
    else
         % for run number > 0, output is printing after every 100 steps
         % We want to remove the first 10^4 outputs = first 100 rows.
        remrows = 100;
        resample2 = round(resample/100);
    end
    D(1:remrows,:)=[];  % ignore the first 10^4 models
   % resample at intervals of 'resample' to ensure independent data pts
    for x = 1:length(D(:,1))
        if mod(x,resample2)==0
            k=k+1;
            % combine all parallel runs
            data(k,:) = D(x,:);
        end
    end
    clearvars D
end
    
