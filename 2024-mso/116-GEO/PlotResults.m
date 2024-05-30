
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Golden Eagle Optimizer (GEO) source codes version 1.0
%  
%  Developed in:	MATLAB 9.6 (R2019a)
%  
%  Programmer:		Abdolkarim Mohammadi-Balani
%  
%  Original paper:	Abdolkarim Mohammadi-Balani, Mahmoud Dehghan Nayeri, 
%					Adel Azar, Mohammadreza Taghizadeh-Yazdi, 
%					Golden Eagle Optimizer: A nature-inspired 
%					metaheuristic algorithm, Computers & Industrial Engineering.
%
%                  https://doi.org/10.1016/j.cie.2020.107050               
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotResults (fun,lb,ub, FunctionNumber,ConvergenceCurve) 

% ----------------

close ('all');
figure ('Position', [469,200,887,395]);

% ----------------

subplot (1,2, 1)
NumPoints = 100;
[SurfX,SurfY] = meshgrid (linspace(lb(1),ub(1),NumPoints),linspace(lb(2),ub(2),NumPoints));
SurfZ = reshape(fun([SurfX(:),SurfY(:)]),NumPoints,NumPoints);
surf (SurfX,SurfY,SurfZ, 'EdgeAlpha',0.2);
box ('on');
title ('Landscape');

% ----------------

subplot (1,2, 2)
stairs (ConvergenceCurve, 'r', 'LineWidth',2);
grid ('on');
ax = gca;
ax.YScale = 'log';
axis tight;
title ('Convergence curve');

% ----------------

sgtitle (sprintf('F%d',FunctionNumber), 'FontWeight','Bold');
