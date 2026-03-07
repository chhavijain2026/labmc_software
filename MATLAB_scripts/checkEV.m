% function, incorporated in perturb_model3.m, meant to check if 
% all parameters identified by 'evpos', such as Ei and/or Vi, are > 0
% in all perturbed models. If parameters at positions marked by evpos are 
% <0 in any model, that row is deleted from the q_pert ensemble. 
% Size of ensemble < M in this case. Finally, more 
% models are added to make size of ensemble = M again.

function q_pert = checkEV(q_pert,M,Km,Am,evpos,eigval,eigvec,meanm,stdm,mBi)

k=0;
i=1;
% i1=1;
while i<length(q_pert(:,1))
    for i1=1:length(evpos)
       if q_pert(i,evpos(i1))<0 %|| q_pert(i,i1+1)<0 || q_pert(i,i2)<0 || q_pert(i,i2+1)<0
          k=k+1;
          f(k)=i;
          break;
       end
    end
    i=i+1;
end
if exist('f','var')
    q_pert(f,:)=[];
    lq = M-k; % new no. of rows in q_pert matrix. # of rows should be = M
    q2 = perturb_model(5e3,Km,Am,evpos,eigval,eigvec,meanm,stdm,mBi,0);
    i=1;
    p=0;
    while p<k
        flag = 0;
        for i1=1:length(evpos)
           if q2(i,evpos(i1))<0 
              flag=1;
              break;
           end
        end
        if flag==0
              p=p+1;
              q_pert(lq+p,:) = q2(i,:);
        end
%        flag=1;
%        for m=1:Km-2
%            if q2(i,m)<0
%                flag=0;
%                break;
%            end;
%        end;
%        if flag==1
%            p=p+1;
%            q_pert(lq+p,:) = q2(i,:);
%        end;
       i=i+1;
    end
end

return;
    
