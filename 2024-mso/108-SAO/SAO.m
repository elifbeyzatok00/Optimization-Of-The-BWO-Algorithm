function [bestSolution, bestFitness, iteration]=SAO(fhd, dimension, maxIteration, fNumber)

config;

nMole=50; % Number of Smell Molecules (Search Agent)
lb=lbArray;
ub=ubArray;
dim=dimension;
Max_iter=maxIteration; % Maximum numbef of iterations
Agent_Pos=zeros(1,dim);
Agent_Fit=inf; %



olf=0.9;
K=0.6;
T=0.95;
M=0.9;
Step=0.02;

%Create the initial position of smell molecules
moles_Pos=initialization(nMole,dim,ub,lb);
BestScore=inf;
Converge_curve=zeros(1,Max_iter);


iter=1;% 

% Main loop
while iter<Max_iter
    for i=1:size(moles_Pos,1)        
%Make Sure smell molecules remains in the search space.
        Clip_ub=moles_Pos(i,:)>ub;
        Clip_lb=moles_Pos(i,:)<lb;
        moles_Pos(i,:)=(moles_Pos(i,:).*(~(Clip_ub+Clip_lb)))+ub.*Clip_ub+lb.*Clip_lb;                      
        % Calculate objective function for each molecules
%         fitness=fobj(moles_Pos(i,:)); 
        fitness=testFunction(moles_Pos(i,:)', fhd, fNumber);
    iter=iter+1;    

        % Agent Fitness
        if fitness<Agent_Fit 
            Agent_Fit=fitness; % Update Agent fitness
            Agent_Pos=moles_Pos(i,:);
        end
    end  
    % Update the Position of molecules
    for i=1:size(moles_Pos,1)
        for j=1:size(moles_Pos,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r3=rand();       
            r4=rand();
            r5=rand();
            Sniff_mole(i,j)=moles_Pos(i,j)+r1*sqrt(3*K*T/M); %Sniffing Mode
        end
%         fitness=fobj(Sniff_mole(i,:));
        fitness=testFunction(Sniff_mole(i,:)', fhd, fNumber);
    iter=iter+1;    

        [~,Index]=min(fitness);
        Agent_Pos=Sniff_mole(:,Index);
        [~,Indes]=max(fitness);
        Worst_Pos=Sniff_mole(:,Indes);
        if fitness<BestScore
            BestScore=fitness;
            moles_Pos(i,:)=Sniff_mole(i,:);
        end
        %Trailing Mode       
        for j=1:size(moles_Pos,2)     
            Trail_mole(i,j)=moles_Pos(i,j)+r3*olf*(moles_Pos(i,j)-Agent_Pos(i,1))...
                -r4*olf*(moles_Pos(i,j)-Worst_Pos(i,1)); %Traili Mode
        end
%         fitness=fobj(Trail_mole(i,:));
        fitness=testFunction(Trail_mole(i,:)', fhd, fNumber);
    iter=iter+1;    

        if fitness<BestScore
            BestScore=fitness;
            moles_Pos(i,:)=Trail_mole(i,:);
        end
         %Random Mode
        for j=1:size(moles_Pos,2)     
            Random_mole(i,j)=moles_Pos(i,j)+r5*Step; 
        end
%         fitness=fobj(Random_mole(i,:));
        fitness=testFunction(Random_mole(i,:)', fhd, fNumber);
    iter=iter+1;    

        if fitness<BestScore
            BestScore=fitness;
            moles_Pos(i,:)=Random_mole(i,:);
        end                  
    end
%     Agent_Fit=BestScore;
Agent_Pos=moles_Pos(1:dim);
    Converge_curve(iter)=Agent_Fit;
bestSolution=Agent_Pos ;
bestFitness=Agent_Fit;
iteration=iter;
end



