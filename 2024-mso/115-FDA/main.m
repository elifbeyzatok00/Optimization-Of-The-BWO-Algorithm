%___________________________________________________________________%
%  Flow Direction Algorithm (FDA): source codes version 1.0         %
%                                                                   %
%  Developed in MATLAB R2017b                                       %
%                                                                   %
%  Authors:H Karami, M Valikhan Anaraki, S Farzin, S. Mirjalili     % 
% programmers: H Karami, M Valikhan Anaraki, S Farzin, S. Mirjalili % 
%                                                                   %
%         e-Mails: h.karami@semnan.ac.ir                            %
%                 mvalikhan@semnan.ac.ir                            %
%                 saeed.farzin@semnan.ac.ir                         %
%                 ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %   
%   Main paper: H Karami, M Valikhan Anaraki, S Farzin, S. Mirjalili%
%               Flow Direction Algorithm (FDA):                     %
%               A Novel Optimization Approach for                   %
%               Solving Optimization Problems,                      %
%               Computers & Industrial Engineering                  %
%                                                                   %
%               DOI: https://doi.org/10.1016/j.cie.2021.107224      %
%                                                                   %
%___________________________________________________________________%
%
% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% alpha = number of flows
% beta= number of neighborhoods
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
%
%
% To run FDA: [Best_fitness,BestX,FDA_cg_curve]=FDA(Max_iteration,lb,ub,dim,fobj,alpha,beta);
%__________________________________________

clear all 
clc

alpha=50; % Number of flows
beta=1; %Number of neighborhoods
Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=200; % Maximum numbef of iterations

% Load details of the selected benchmark function
[dim,fobj,ub, lb]  = Select_Functions(Function_name);

[Best_fitness,BestX,FDA_cg_curve]=FDA(Max_iteration,lb,ub,dim,fobj,alpha,beta);


figure('Position',[500 500 660 290])
%Draw objective space
semilogy(FDA_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('FDA')

display(['The best solution obtained by DFA is : ', num2str(BestX)]);
display(['The best optimal value of the objective funciton found by DFA is : ', num2str(Best_fitness)]);

        



