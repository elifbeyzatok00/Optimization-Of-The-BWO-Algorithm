%___________________________________________________________________%
%  Smell Agent Optimzation (SAO) source codes version 1.0.1               %
%                                                                   %
%  Developed in MATLAB R2020b                                 %
%                                                                   %
%  Author and programmer: Salawudeen Ahmed Tijani                        %
%                                                                   %
%         e-Mail: atsalawudeen@unijos.edu.ng                           %
%                      tasalawudeen@abu.edu.ng             %
%                                                                   %

%   Main paper: Salawudeen A. T., Mu'azu M. B., Sha'aban Y. A., Adedokun E.A. %
%               A Novel Smell Agent Optimization: An Extensive CEC Study and Engineering Application        %
%               DOI: 10.1016/j.knosys.2021.107486              %
%                                                                   %
%___________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% myCost = @YourCostFunction
% dim = number of your variables
% Max_Iter = maximum number of generations
% nMole = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run SAO: [Smell_Object,Best_mole,Convergence]=SAO(nMole,Max_Iter,lb,ub,dim,fobj);
%__________________________________________

clear all
close all
clc
format long
nMole=50; % Number of Smell Molecules (Search Agent)

F_name='F1'; % Selecte Benchmark Function

Max_Iter=5000; % Maximum numbef of iterations
% run50;
% Load Function Details
[lb,ub,dim,myCost]=Select_Function(F_name);
% for k=1:run
    [Smell_Object,Best_mole,Convergence]=SAO(nMole,Max_Iter,lb,ub,dim,myCost);
%     Smell_Object(k,:)=Smell_Object;
% end

BestCost=Smell_Object

figure('Position',[500 500 660 290])
%Draw objective space
subplot(1,2,2);
semilogy(Convergence,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('SAO')
%Draw function in hyperspace
subplot(1,2,1);
func_plot(F_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_i , x_j )'])

%% Uncomment the lines Bellow to obtained the statistical analysis
% Best=min(Smell_Object)
%  Worst=max(Smell_Object)
% Average=mean(Smell_Object)
% STD=std(Smell_Object)







