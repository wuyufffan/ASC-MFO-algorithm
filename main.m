  
clear all 
clc
func_number=8;
runs=30;
SearchAgents_no=50;
max_nfes = 50000;

% kj=func_number*(runs+1);


k=1;
% kl=zeros(30,1);
Max_iteration = max_nfes/SearchAgents_no;
Algorithm_name={'SCA_MFO_02'};
for k = 1:1
    cg_curve = [];
    FUNC = Algorithm_name{k};
    fprintf('Algorithm =\t %d\n',k);
    for i=1:8
            fprintf('Problem =\t %d\n',i);
            Function_name=i;
            [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
            for j=1:runs
                fprintf('run =\t %d\n',j);
                tic;
                [Best_score1,Best_pos1,Convergence_curve]=feval(FUNC,SearchAgents_no,max_nfes,lb,ub,dim,fobj);
                t(1,j) = toc;
                gbest_val(1,j) = Best_score1;
                gbest_pos(j,:) = Best_pos1;
                cg_curve(j,:)=Convergence_curve(1,:);

            end

            gbest_all_val{1,i} =gbest_val;
            gbest_all_pos{1,i} =gbest_pos;
            cg_curve_all{1,i} = cg_curve;
            time_all{1,i} = t;    
    end
    
    if k==1
        save('SCA_MFO.mat','gbest_all_val','gbest_all_pos','cg_curve_all','time_all');
    end

end

