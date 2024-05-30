
%___________________________________________________________________%
%  Golden Jackal Optimization  (GJO)            
%  Developed in MATLAB R2018a                                       
%  Authors: Nitish Chopra and Muhammad Mohsin Ansari
%  Programmer: Nitish Chopra                            

%   Main paper: Chopra, Nitish, and Muhammad Mohsin Ansari. "Golden Jackal Optimization: A
%              Novel Nature-Inspired Optimizer for Engineering Applications." 
%              Expert Systems with Applications (2022): 116924. 
%
%               DOI: https://doi.org/10.1016/j.eswa.2022.116924             
%                                                                   %
%___________________________________________________________________%
%% TESTING GJO ON ENGINEERING DESIGN
SearchAgents_no=50; % Number of search agents
Max_iteration=3000; % Maximum numbef of iterations

Function_name='F1';
% 'F1'Tension/compression spring design; 
% 'F2' %Pressure vessel design
% 'F3' %Welded beam design 
% 'F4' %Speed Reducer design 
% 'F5'  Gear train design problem
% 'F6' Three-bar truss design problem

%% Load details of the selected engineering design problem
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
runn=5;% maximum number of re-run of GJO
cost=zeros(runn,1);pos=zeros(runn,4);
for i=1:runn
    disp(['Run no: ',num2str(i)]);
 [Male_Jackal_score,Male_Jackal_pos,GJO_cg_curve]=GJO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
cost(i,:)=Male_Jackal_score;
end
mean_cost=mean(cost);
min_cost=min(cost);
max_cost=max(cost);

disp(['best value GJO:  ',num2str(min_cost,10),'  Mean: ', num2str(mean_cost),'  Max: ', num2str(max_cost)]);
    

