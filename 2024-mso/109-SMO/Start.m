% =========================================================================
%     Starling murmuration optimizer: A novel bio-inspired algorithm
%                   for global and engineering optimization
%        Computer Methods in Applied Mechanics and Engineering
%                    Volume 392, 15 March 2022, 114616
%   https://www.sciencedirect.com/science/article/abs/pii/S0045782522000330?via%3Dihub 
%    ---------------------------------------------------------------------
%  Authors: HodaZamani, Mohammad H.Nadimi-Shahraki*, Amir H.Gandomi
%  *Corresponding author: Faculty of Computer Engineering, Najafabad 
%                Branch, Islamic Azad University, Najafabad, Iran.
%  e-Mails: nadimi@ieee.org, zamanie_hoda@ymail.com, 
%           nadimi.mh@gmail.com, a.h.gandomi@gmail.com
% =========================================================================
clc 
clear
tic 
fprintf('================================================================\n');
fprintf(' Starling Murmuration Optimizer: A Novel Bio-inspired Algorithm \n');
fprintf('           for global and engineering optimization \n');
fprintf(' ---------------------------------------------------------------\n');
warning off
format long
%% Problem defination
fhd = @cec17_func;
Problem_Size  = 10;    % Dimensions D  = 2, 10, 30, 50, 100
lu = [-100 * ones(1, Problem_Size); 100 * ones(1, Problem_Size)];
Alg_run = 30;
%% SMO parameter setting
PopSize    = 100;          % Maximum number of population
MaxIt      = Problem_Size *10000/PopSize;   % Maximum number of iterations
k  = 10;               % Number of flocks (range k = 3, 5, 10, 15, 20) 
%% Start
for Func = [1,3:30]
    
    Optimum = Func * 100.0;
    fprintf('\n-------------------------------------------------------\n');
    fprintf('Function = %d, Dimension size = %d\n', Func, Problem_Size);
    for Run =1: Alg_run
        [Covergence, Best_pos,Best_Fit]= SMO(PopSize,Problem_Size,MaxIt,k,lu,Func);
        fprintf('run = %d, Best_Fit= %d\n', Run, Best_Fit);
        SMO_Results(Func).Covergence(Run,:) = Covergence;
        SMO_Results(Func).Pos(Run,:) = Best_pos;
        SMO_Results(Func).Fit(Run,:) = Best_Fit;
    end
end


