%__________________________________________________________________
%  Dandelion Optimizer
%  Developed in MATLAB R2018a
%
%  programmer: Shijie Zhao and Tianran Zhang   
%  E-mail: zhaoshijie@lntu.edu.cn
%          ztr20010118@126.com
%  The code is based on the following papers.
%  Shijie Zhao, Tianran Zhang, Shilin Ma, and Miao Chen 
%  Dandelion Optimizer: A nature-inspired metaheuristic algorithm for
%  engineering applications.
%  Engineering Applications of Artificial Intelligence
%  DOI:10.1016/j.engappai.2022.105075
%
%__________________________________________________________________

clear all 
clc

N=30; % Number of search agents
Max_iter=500; % Maximum numbef of iterations
F_name='F12'; % Name of the test function

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(F_name);

tic;
[Bestfitness,Bestposition,Convergencecurve]=DO(N,Max_iter,lb,ub,dim,fobj);
Run_time=toc;
semilogy(1:Max_iter,Convergencecurve,'color','r','linewidth',2.5);
title('Convergence curve');
xlabel('Iteration');
ylabel('Best score obtained so far')
display(['The running time is:', num2str(Run_time)]);
display(['The best fitness is:', num2str(Bestfitness)]);
display(['The best position is: ', num2str(Bestposition)]);

