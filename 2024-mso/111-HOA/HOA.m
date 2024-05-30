function [bestSolution, bestFitness, iteration]=HOA(fhd, dimension, maxIteration, fNumber)

config;

%% Problem Definition

global nVar;
nVar=dimension;           % Number of Decision Variables

VarSize=[1 nVar];  % Size of Decision Variables Matrix

VarMin=lbArray;        % Lower Bound of Variables
VarMax=ubArray;         % Upper Bound of Variables

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
%% Algorithm Parameters

MaxIt=maxIteration;                 % Maximum Number of Iterations
nHourse=50;                % Number of Mean Points

%% Initialization
tic

empty_Horse.Position=[];
empty_Horse.Cost=[];
empty_Horse.Velocity=[];
empty_Horse.Best.Position=[];
empty_Horse.Best.Cost=[];

Hourse=repmat(empty_Horse,nHourse,1);

GlobalBest.Cost=inf;

CostPositionCounter=zeros(nHourse,2+nVar);

for i=1:nHourse
    
    % Initialize Position
    Hourse(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    Hourse(i).Velocity=zeros(VarSize);
    
    % Evaluation
%     Hourse(i).Cost=CostFunction(Hourse(i).Position);
    Hourse(i).Cost=testFunction(Hourse(i).Position', fhd, fNumber);

    % Update Personal Best
    Hourse(i).Best.Position=Hourse(i).Position;
    Hourse(i).Best.Cost=Hourse(i).Cost;
    
    % Update Global Best
    if Hourse(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=Hourse(i).Best;
        
    end
    
    CostPositionCounter(i,:)=[i Hourse(i).Best.Cost Hourse(i).Best.Position];
end

BestCost=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% HOA Main Loop
tic

w=1;
phiD=0.02;
phiI=0.02;


g_Alpha=1.50;       % Grazing
d_Alpha=0.5;        % Defense Mechanism
h_Alpha=1.5;        % Hierarchy

g_Beta=1.50;       % Grazing
h_Beta=0.9;        % Hierarchy 
s_Beta=0.20;       % Sociability
d_Beta=0.20;       % Defense Mechanism

g_Gamma=1.50;       % Grazing 
h_Gamma=0.50;       % Hierarchy 
s_Gamma=0.10;    % Sociability 
i_Gamma=0.30;       % Imitation
d_Gamma=0.10;       % Defense Mechanism 
r_Gamma=0.05;       % Random (Wandering and Curiosity)

g_Delta=1.50;       % Grazing
r_Delta=0.10;       % Random (Wandering and Curiosity) 
it=1;
while it<MaxIt
    
    CostPositionCounter=sortrows(CostPositionCounter,2);
    MeanPosition=mean(CostPositionCounter(1:nHourse,2));
    BadPosition=mean(CostPositionCounter((1-phiD)*nHourse:nHourse,3));
    GoodPosition=mean(CostPositionCounter(1:phiI*nHourse,3));
    
    for i=1:nHourse
        
        CC=find(CostPositionCounter==i); %CC is Cost Counter
        
        % Update Velocity
        if CC<=0.1*nHourse

           Hourse(i).Velocity = +h_Alpha*rand(VarSize).*(GlobalBest.Position-Hourse(i).Position)...
                                -d_Alpha*rand(VarSize).*(Hourse(i).Position)...
                                +g_Alpha*(0.95+0.1*rand)*(Hourse(i).Best.Position-Hourse(i).Position);
               
        elseif CC<=0.3*nHourse
           Hourse(i).Velocity = s_Beta*rand(VarSize).*(MeanPosition-Hourse(i).Position)...
                               -d_Beta*rand(VarSize).*(BadPosition-Hourse(i).Position)...
                               +h_Beta*rand(VarSize).*(GlobalBest.Position-Hourse(i).Position)...
                               +g_Beta*(0.95+0.1*rand)*(Hourse(i).Best.Position-Hourse(i).Position);
        
        elseif CC<=0.6*nHourse
           Hourse(i).Velocity = s_Gamma*rand(VarSize).*(MeanPosition-Hourse(i).Position)...
                               +r_Gamma*rand(VarSize).*(Hourse(i).Position)...
                               -d_Gamma*rand(VarSize).*(BadPosition-Hourse(i).Position)...
                               +h_Gamma*rand(VarSize).*(GlobalBest.Position-Hourse(i).Position)...
                               +i_Gamma*rand(VarSize).*(GoodPosition-Hourse(i).Position)...
                               +g_Gamma*(0.95+0.1*rand)*(Hourse(i).Best.Position-Hourse(i).Position);
        
        else
           Hourse(i).Velocity = +r_Delta*rand(VarSize).*(Hourse(i).Position)...
                                +g_Delta*(0.95+0.1*rand)*(Hourse(i).Best.Position-Hourse(i).Position);
        end
        
        % Apply Velocity Limits
        Hourse(i).Velocity = max(Hourse(i).Velocity,VelMin);
        Hourse(i).Velocity = min(Hourse(i).Velocity,VelMax);
        
        % Update Position
        Hourse(i).Position = Hourse(i).Position + Hourse(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(Hourse(i).Position<VarMin | Hourse(i).Position>VarMax);
        Hourse(i).Velocity(IsOutside)=-Hourse(i).Velocity(IsOutside);
        
        % Apply Position Limits
        Hourse(i).Position = max(Hourse(i).Position,VarMin);
        Hourse(i).Position = min(Hourse(i).Position,VarMax);
        
        % Evaluation
%         Hourse(i).Cost = CostFunction(Hourse(i).Position);
        Hourse(i).Cost=testFunction(Hourse(i).Position', fhd, fNumber);
          it=it+1;
        % Update Personal Best
        if Hourse(i).Cost<Hourse(i).Best.Cost
            
            Hourse(i).Best.Position=Hourse(i).Position;
            Hourse(i).Best.Cost=Hourse(i).Cost;
            
            % Update Global Best
            if Hourse(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=Hourse(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    nfe(it)=it;
    
%     disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost = ' num2str(BestCost(it))]);

% cD1=cD1*w;
% cD2=cD2*w; cS2=cS2*w; cR2=cR2*w; cI2=cI2*w;
% cD3=cD3*w; cS3=cS3*w; cR3=cR3*w; cI3=cI3*w;
% cR4=cR4*w; cI4=cI4*w;

d_Alpha=d_Alpha*w; g_Alpha=g_Alpha*w;
d_Beta=d_Beta*w; s_Beta=s_Beta*w; g_Beta=g_Beta*w;
d_Gamma=d_Gamma*w; s_Gamma=s_Gamma*w; r_Gamma=r_Gamma*w; i_Gamma=i_Gamma*w;
g_Gamma=g_Gamma*w;
r_Delta=r_Delta*w; g_Delta=g_Delta*w;

end

%% Results
% disp([ ' Best Position = '  num2str(GlobalBest.Position)])
% disp([ ' Best Cost = '  num2str(GlobalBest.Cost)])
% disp([ ' Time = '  num2str(toc)])
% plot(nfe,BestCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');

bestSolution=GlobalBest.Cost;
bestFitness= Hourse(i).Cost ;
iteration=it;
end
