%  Adaptive Sine-cosine Moth-Flame Optimization Algorithm£¨ASC-MFO algorithm£©                                                           
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2019a                                                                
%                                                                                                     
%  programmer: YuFan Wu
%  e-Mail: 1095978552@qq.com
function [Best_score1,Best_pos1,Convergence_curve]=SCA_MFO_02(NP,Max_FES,lb,ub,dim,fobj)  
format long; 
display('SCA_MFO_02 is optimizing your problem');
rand('seed', sum(100 * clock));
N = NP;
Xmin = lb;
Xmax = ub;
D = dim;
dimension = D;
popsize=NP;
% counter_g1=0;
% counter_g2=0;


X = repmat(Xmin, popsize, 1) + rand(popsize, D) .* (repmat(Xmax-Xmin, popsize, 1));
fitcount =0; 
 maxFES = Max_FES;
 Max_iteration = maxFES/NP;
 
numg1 = 25;

     
for i=1:size(X,1)
    Val_X(1,i)=fobj(X(i,:));
    fitcount = fitcount + 1;
end
     
[gbest_fitness,g_index]=min(Val_X);
gbest_pos=X(g_index,:);

pbest_pos=X;
pbest_fitness=Val_X;

k = 0;
a =2;
while k<=Max_iteration && fitcount<=Max_FES
    k=k+1;
%     mean_val_x = sum(Val_X)/NP;      
%     diversity = sum((Val_X-mean_val_x).^2)/NP;
%     D_all(k) = diversity;
    r1=a-k*((a)/Max_iteration);
        for i =1:numg1
            for d = 1:D
                    r2=(2*pi)*rand();
                    r3=rand;
                    r4=rand();
                    if r4<0.5
                % Eq. (3.1)
                        X(i,d)= X(i,d)+(r1*sin(r2)*abs(r3*gbest_pos(1,d)-X(i,d)));
                    else
                % Eq. (3.2)
                        X(i,d)= X(i,d)+(r1*cos(r2)*abs(r3*gbest_pos(1,d)-X(i,d)));
                    end
            end
                    X(i,:)=boundConstraint_absorb(X(i,:), Xmin, Xmax);
                    Val_X(1,i)=fobj(X(i,:));
                    fitcount = fitcount + 1;
                    if fitcount>=Max_FES
                        break;
                    end
                    if  Val_X(1,i)<pbest_fitness(1,i)
                        pbest_fitness(1,i) = Val_X(1,i);
                        pbest_pos(i,:) =X(i,:);
%                         counter_g1 = counter_g1 + 1;
                        if Val_X(1,i)<gbest_fitness
                            gbest_fitness = Val_X(1,i);
                            gbest_pos = X(i,:);
                        end
                    end
        end
%         mean_val_g1= sum(Val_X(1,1:numg1))/numg1;
%         diversity_g1 = sum((Val_X(1,1:numg1)-mean_val_g1).^2)/numg1;
%         D_g1(k) = diversity_g1;
        [~,index] = sort(pbest_fitness);
        flame=(pbest_pos(index(1:(NP-numg1)),:));
        for i=numg1+1:NP
            a=-1+k*((-1)/Max_iteration);
            distance_to_flame_temp=abs(X(i,:)-flame(i-numg1,:));
            bb=1;
            t=(a-1)*rand(1,dimension)+1;
            X(i,:)=distance_to_flame_temp.*exp(bb.*t).*cos(t.*2*pi)+flame(i-numg1,:);
% end


                
            X(i,:)=boundConstraint_absorb(X(i,:), Xmin, Xmax);

            Val_X(1,i)=fobj(X(i,:));
            fitcount = fitcount + 1;
%         mean_val_g2= sum(Val_X(1,numg1+1:NP))/(NP-numg1);
%         diversity_g2 = sum((Val_X(1,numg1+1:NP)-mean_val_g2).^2)/(NP-numg1);
%         D_g2(k) = diversity_g2;
        
            if fitcount>=Max_FES
                break;
            end

            if  Val_X(1,i)<pbest_fitness(1,i)
                    pbest_fitness(1,i) = Val_X(1,i);
                    pbest_pos(i,:) =X(i,:);
%                     counter_g2 = counter_g2+1;
                    if Val_X(1,i)<gbest_fitness
                        gbest_fitness = Val_X(1,i);
                        gbest_pos = X(i,:);
                    end
            end
        end
%         if (counter_g1/numg1)>= counter_g2/(NP - numg1)
%             numg1 = numg1 + 1;
%         else
%             numg1 = numg1 - 1;
%         end
%         if numg1 < 5
%             numg1 = 5;
%         end
%         if numg1 >45
%             numg1 = 45;
%         end
%         counter_g1 = 0;
%         counter_g2 = 0;
%         if mod(fitcount,10)==0
%             display(['At FES ', num2str(fitcount), ' the best fitness is ', num2str(gbest_fitness)]);
%         end
% X1 = X(:,2);
% X2 = X(:,3);
% x = 0.22728;
% y = 55.4612;
% figure(1);
% plot(X1,X2,'go','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
% hold on
% plot(x,y,'cs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',10);
% xlabel('Isd1','FontSize',15);ylabel('Rsh','FontSize',15);
% hold on
% axis([Xmin(1,2) Xmax(1,2) Xmin(1,3) Xmax(1,3)]);
% grid on 
% T = legend('particle','gbest');
% set(T,'Fontsize',12); 
% hold off
% grid off
Convergence_curve(1,k) = gbest_fitness;
end



Best_score1 = gbest_fitness;
Best_pos1 = gbest_pos;
end
            
            
            
            
            
        
    
    