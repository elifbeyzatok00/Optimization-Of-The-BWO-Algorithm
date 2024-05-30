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
% This function initialize the first population of search agents
function [flow_x,fitness_flow]=initialization(alpha,dim,ub,lb)

for i=1:alpha
    flow_x(i,:)=lb+rand(1,dim).*(ub-lb);%position of each flow
end