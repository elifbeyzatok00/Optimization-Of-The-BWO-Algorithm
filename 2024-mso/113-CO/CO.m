function [bestSolution, bestFitness, iteration]=CO(fhd, dimension, maxIteration, fNumber)

config;

tic
for fnm = 1 % Test functions
    for run = 1 : 1
lb=lbArray;
ub=ubArray;
D=dimension;
        %% Problem definition (Algorithm 1, L#1)
%         if fnm <= 9
%             Function_name = f_name(fnm,:);
%         else
%             Function_name = f_name1(fnm-9,:);
%         end
%         [lb,ub,D,fobj] = Get_Functions_details(Function_name); % The shifted functions' information
%         D = 100;                   % Number of Decision Variables
%         if length(lb) == 1
%             ub = ub.*ones(1,D);   % Lower Bound of Decision Variables
%             lb = lb.*ones(1,D);   % Upper Bound of Decision Variables
%         end
        
        n = 6;                   % Population Size
        m = 2;                    % Number of search agenets in a group
        
        %% Generate initial population of cheetahs (Algorithm 1, L#2)
        empty_individual.Position = [];
        empty_individual.Cost = [];
        BestSol.Cost = inf;
        pop = repmat(empty_individual,n,1);
        
        for i=1:n
            pop(i).Position = lb+rand(1,D).*(ub-lb);
%             pop(i).Cost = fobj(pop(i).Position);
            pop(i).Cost=testFunction(pop(i).Position', fhd, fNumber);
            if pop(i).Cost < BestSol.Cost
                BestSol = pop(i); % Initial leader position
            end
        end
        
        %% Initialization (Algorithm 1, L#3)
        pop1 = pop;               % Population's initial home position
        BestCost = [];            % Leader fittnes value in a current hunting period
        X_best = BestSol;         % Prey solution sofar
        Globest = BestCost;       % Prey fittnes value sofar
        
        %% Initial parameters
        t = 0;                    % Hunting time counter (Algorithm 1, L#4)
        it = 1;                   % Iteration counter(Algorithm 1, L#5)
        MaxIt =maxIteration;           % Maximum number of iterations (Algorithm 1, L#6)
        T = ceil(D/10)*60;        % Hunting time (Algorithm 1, L#7)
        FEs = 0;                  % Counter for function evaluations (FEs)
        %% CO Main Loop
        while FEs <= MaxIt % Algorithm 1, L#8
            %  m = 1+randi (ceil(n/2)); 
            i0 = randi(n,1,m);    % select a random member of cheetahs (Algorithm 1, L#9)
            for k = 1 : m % Algorithm 1, L#10
                i = i0(k);
                
                % neighbor agent selection (Algorithm 1, L#11)
                if k == length(i0)
                    a = i0(k-1);
                else
                    a = i0(k+1);
                end
                
                X = pop(i).Position;    % The current position of i-th cheetah
                X1 = pop(a).Position;   % The neighbor position
                Xb = BestSol.Position;  % The leader position
                Xbest = X_best.Position;% The pery position
                             
                kk=0;
                % Uncomment the follwing statements, it may improve the performance of CO
                                if i<=2 && t>2 && t>ceil(0.2*T+1) && abs(BestCost(t-2)-BestCost(t-ceil(0.2*T+1)))<=0.0001*Globest(t-1)
                                    X = X_best.Position;
                                    kk = 0;
                                elseif i == 3
                                    X = BestSol.Position;
                                    kk = -0.1*rand*t/T;
                                else
                                    kk = 0.25;
                                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                if mod(it,100)==0 || it==1
                   xd = randperm(numel(X));
                end
                Z = X;
                
                %% Algorithm 1, L#12
                for j = xd % select arbitrary set of arrangements
                    %% Algorithm 1, L#13
                    r_Hat = randn;         % Randomization paameter, Equation (1)
                    r1 = rand;
                    if k == 1              % The leader's step length (it is assumed that k==1 is associated to the leade number)
                        alpha = 0.0001*t/T.*(ub(j)-lb(j)); % Step length, Equation (1) %This can be updated by desired equation 
                    else                   % The members' step length
                        alpha = 0.0001*t/T*abs(Xb(j)-X(j))+0.001.*round(double(rand>0.9));%member step length, Equation (1)%This can be updated by desired equation
                    end
                    
                    r = randn;
                    r_Check = abs(r).^exp(r/2).*sin(2*pi*r); % Turning factor, Equation (3)%This can be updated by desired equation
                    beta = X1(j)-X(j);     % Interaction factor, Equation (3)
                    
                    h0 = exp(2-2*t/T);
                    
                    H = abs(2*r1*h0-h0);
                    
                    %% Algorithm 1, L#14
                    
                    r2 = rand;
                    r3 = kk+rand;
                    
                    %% Strategy selection mechanism
                    if r2 <= r3              % Algorithm 1, L#15
                        r4 = 3*rand;         % Algorithm 1, L#16
                        if H > r4            % Algorithm 1, L#17
                            Z(j) = X(j)+r_Hat.^-1.*alpha;    % Search, Equation(1) (Algorithm 1, L#18)
                        else
                            Z(j) = Xbest(j)+r_Check.*beta;    % Attack, Equation(3) (Algorithm 1, L#20)
                        end
                    else
                        Z(j) = X(j);         % Sit&wait, Equation(2) (Algorithm 1, L#23)
                    end
                end
                %% Update the solutions of member i (Algorithm 1, L#26)
                % Check the limits
                xx1=find(Z<lb);
                Z(xx1)=lb(xx1)+rand(1,numel(xx1)).*(ub(xx1)-lb(xx1));
                xx1=find(Z>ub);
                Z(xx1)=lb(xx1)+rand(1,numel(xx1)).*(ub(xx1)-lb(xx1));
                
                % Evaluate the new position
                NewSol.Position = Z;
%                 NewSol.Cost = fobj(NewSol.Position);
                NewSol.Cost=testFunction(NewSol.Position', fhd, fNumber);
                if NewSol.Cost < pop(i).Cost
                    pop(i) = NewSol;
                    if pop(i).Cost < BestSol.Cost
                        BestSol = pop(i);
                    end
                end
                FEs = FEs+1;
            end
            
            t = t+1; % (Algorithm 1, L#28)
            
            %% Leave the prey and go back home (Algorithm 1, L#29)
            if t>T && t-round(T)-1>=1 && t>2
                if  abs(BestCost(t-1)-BestCost(t-round(T)-1))<=abs(0.01*BestCost(t-1))
                    
                    % Change the leader position (Algorithm 1, L#30)
                    best = X_best.Position;
                    j0=randi(D,1,ceil(D/10*rand));
                    best(j0) = lb(j0)+rand(1,length(j0)).*(ub(j0)-lb(j0));
%                     BestSol.Cost = fobj(best);
                    BestSol.Cost=testFunction(best', fhd, fNumber);
                    BestSol.Position = best; % Leader's new position 
                    FEs = FEs+1;
                    
                    i0 = randi(n,1,round(1*n));
                    % Go back home, (Algorithm 1, L#30)
                    pop(i0(n-m+1:n)) = pop1(i0(1:m)); % Some members back their initial positions 
                    
                    pop(i) = X_best; % Substitude the member i by the prey (Algorithm 1, L#31)
                    
                    t = 1; % Reset the hunting time (Algorithm 1, L#32)
                end
            end
            
            it = it +1; % Algorithm 1, L#34
            
            %% Update the prey (global best) position (Algorithm 1, L#35)
            if BestSol.Cost<X_best.Cost
                X_best=BestSol;
            end
            BestCost(t)=BestSol.Cost;
            Globest(1,t)=X_best.Cost;
            
            %% Display
            if mod(it,500)==0
%             disp([' FEs>> ' num2str(FEs) '   BestCost = ' num2str(Globest(t))]);
            end
           
        end
        
%         fit(run) = Globest(end) ;
%         mean_f = mean (fit);
%         std_f = std (fit);
%         clc
        
    end
end
% toc
bestSolution=X_best.Cost;
bestFitness= BestSol.Cost;
iteration=FEs;
end
