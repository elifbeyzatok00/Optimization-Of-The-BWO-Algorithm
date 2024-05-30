function [bestSolution, bestFitness, iteration]=edo(fhd, dimension, maxIteration, fNumber)

config;
    
N = 30;
Max_iter = ceil(maxIteration / N);
LB = lbArray;
UB = ubArray;
Dim = dimension;

BestSol=zeros(1,Dim);   % the optimal solution
BestFitness = inf;      % the optimal Fitness value
Xwinners=initialization(N,Dim,UB,LB);    % Initialise the population X using Eq. (12)
Fitness=zeros(1,size(Xwinners,1));
for i=1:N
    Fitness(1,i) = testFunction(Xwinners(i,:)', fhd, fNumber);
    if Fitness(1,i) < BestFitness
        BestFitness=Fitness(1,i);
        BestSol = Xwinners(i,:);
    end
end
Memoryless=Xwinners;%old population
iter=1;
while iter<Max_iter+1
    V=zeros(N,Dim);
    cgcurve(iter)=BestFitness;
    
    %rank the solutions according to the fitness function sort fitness by indices
    [Fitness,sorted_indices]=sort(Fitness);
    temp_Xwinners=Xwinners;
    Xwinners=temp_Xwinners(sorted_indices,:);
    
    %EDO parameters
    d=(1-iter/Max_iter);  %Eq.(23)
    f= 2*rand-1; %Eq. (17)
    a=f^10;      %Eq.(15)
    b=f^5;       %Eq.(16)
    c=d*f;       %Eq. (22)
    sum=(Xwinners(1,:)+Xwinners(2,:)+Xwinners(3,:));
    X_guide=(sum/3);   % Eq.(13) calculates the guiding solution using the average of the best two solutions (x(1,:) x(2,:) x(3,:)  mean(X(i,:))
    
    for  i=1:N
        alpha=rand;
        if alpha<0.5
            
            %---------------- Exploiation -------------------------------------------------------------
            if Memoryless(i,:)==Xwinners(i,:)
                Mu=(X_guide+Memoryless(i,:))/2.0;    %  Eq. (19) is the mean or expected value of an exponentially distributed random vaiable X
                ExP_rate=1./Mu;    % Eq.(18) calculates the rate of the exponential distribution vector lambda
                variance=1./ExP_rate.^2;    %Eq. (4) calculates the vaiance vector
                V(i,:)=a.*(Memoryless(i,:)-variance)+b.*X_guide;   %Eq. (14) first branch
            else
                Mu=(X_guide+Memoryless(i,:))/2.0;    %Eq. (19)is the mean or expected value of an exponentially distributed random vaiable X
                ExP_rate=1./Mu;    %Eq.(18) calculates the rate of the exponential distribution vector lambda
                variance=1./ExP_rate.^2;    %Eq. (4) calculates the vaiance vector
                phi=rand;%random numbers between 0 and 1
                V(i,:)=b.*(Memoryless(i,:)-variance)+log (phi).*Xwinners(i,:);% Eq. (14) second branch
            end
        else
            %----------------------------Exploration----------------------------------------
            M=mean(Xwinners);
            s=randperm(N);
            D1=M-Xwinners(s(1),:); %Eq.(26)
            D2=M-Xwinners(s(2),:); %Eq.(27)
            Z1=M-D1+D2; %Eq.(24)
            Z2=M-D2+D1; %Eq.(25)
            V(i,:)=(Xwinners(i,:)+(c.*Z1+(1-c).*Z2)-M); %Eq.(20)
        end
        %check bounds of V
        F_UB=V(i,:)>UB;
        F_LB=V(i,:)<LB;
        V(i,:)=(V(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
    end
    
    for i=1:N
        for j=1:Dim
            Memoryless(i,j)=V(i,j);
        end
    end
     
    %generate the next population
    V_Fitness=zeros(1,size(Xwinners,1));
    for i=1:N
        V_Fitness(1,i) = testFunction(V(i,:)', fhd, fNumber);
        if V_Fitness(1,i)<Fitness(1,i)
            Xwinners(i,:)=V(i,:);
            Fitness(1,i)=V_Fitness(1,i);
            if Fitness(1,i) < BestFitness
                BestFitness=Fitness(1,i);
                BestSol = Xwinners(i,:);
            end
        end
    end
    iter=iter+1;
end
bestSolution = BestSol;
bestFitness = BestFitness;
iteration = 1;
end

function Positions=initialization(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end