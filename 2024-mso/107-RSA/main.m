%_______________________________________________________________________________________%
%  Reptile Search Algroithm (RSA) source codes demo version 1.0                         %
%                                                                                       %
%  Developed in MATLAB R2015a (7.13)                                                    %
%                                                                                       %
%  Author and programmer: Laith Abualigah                                               %
%                                                                                       %
%         e-Mail: Aligah.2020@gmail.com                                                 %
%       Homepage:                                                                       %
%         1- https://scholar.google.com/citations?user=39g8fyoAAAAJ&hl=en               %
%         2- https://www.researchgate.net/profile/Laith_Abualigah                       %
%_______________________________________________________________________________________%
%  Main paper:            Reptile Search Algorithm (RSA):                               %
%                  A novel nature-inspired meta-heuristic optimizer                     %                                                                       %
%_______________________________________________________________________________________%

clear all 
clc

Solution_no=10;  %Number of search solution
F_name='F10';     %Name of the test function
T=500;           %Maximum number of iterations

    
[LB,UB,Dim,F_obj]=Get_F(F_name); %Give details of the underlying benchmark functions
[Best_F,Best_P,Conv]=RSA(Solution_no,T,LB,UB,Dim,F_obj); % Call Reptile Search Algorithm (RSA)



figure('Position',[454   445   694   297]);
subplot(1,2,1);
func_plot(F_name);     % Function plot
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])
subplot(1,2,2);       % Convergance plot
plot(Conv,'LineWidth',3)
xlabel('Iteration#');
ylabel('Best fitness function');
legend('RSA');



display(['The best-obtained solution by RSA is : ', num2str(Best_P)]);  
display(['The best optimal value of the objective funciton found by RSA is : ', num2str(Best_F)]);  