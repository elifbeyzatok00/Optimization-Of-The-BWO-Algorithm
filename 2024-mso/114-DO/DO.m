function [bestSolution, bestFitness, iteration]=DO(fhd, dimension, maxIteration, fNumber)

config;

LB=lbArray;
UB=ubArray;
% N=30; % Number of search agents
Popsize=30;
Maxiteration=maxIteration; % Maximum numbef of iterations
Dim=dimension;
tic;

dandelions=initialization(Popsize,Dim,UB,LB);
dandelionsFitness = zeros(1,Popsize);
Convergence_curve=zeros(1,Maxiteration);
for i=1:Popsize
%     dandelionsFitness(1,i)=Fobj(dandelions(i,:));
    dandelionsFitness(1,i)=testFunction(dandelions(i,:)', fhd, fNumber);
   
end
% Calculate the fitness values of initial dandelions.
[~,sorted_indexes]=sort(dandelionsFitness);
Best_position=dandelions(sorted_indexes(1),:);
Best_fitness = dandelionsFitness(sorted_indexes(1));

 Maxiteration=Maxiteration/Popsize;

for t=1:Maxiteration
    
    %% Rising stage
    beta=randn(Popsize,Dim);
    alpha=rand()*((1/Maxiteration^2)*t^2-2/Maxiteration*t+1); % eq.(8) in this paper
    a=-1/(Maxiteration^2-2*Maxiteration+1);
    b=-2*a;
    c=1-a-b;
    k=1-rand()*(c+a*t^2+b*t); % eq.(11) in this paper
    if randn()<1.5
        for i=1:Popsize
            lamb=abs(randn(1,Dim));
            theta=(2*rand()-1)*pi;
            row=1/exp(theta);
            vx=row*cos(theta);
            vy=row*sin(theta);
            NEW=rand(1,Dim).*(UB-LB)+LB;
            dandelions_1(i,:)=dandelions(i,:)+alpha.*vx.*vy.*lognpdf(lamb,0,1).*(NEW(1,:)-dandelions(i,:)); % eq.(5) in this paper
        end
    else
        for i=1:Popsize
            dandelions_1(i,:)=dandelions(i,:).*k; % eq.(10) in this paper
            
        end
    end
    dandelions=dandelions_1;
    % Check boundries
    for i=1:Popsize
        dandelions(i) = max(dandelions(i, :),LB);
        dandelions(i) = min(dandelions(i, :),UB);
    end
    
    %% Decline stage
    dandelions_mean=sum(dandelions,1)/Popsize; % eq.(14) in this paper
    for i=1:Popsize
        for j=1:Dim
            dandelions_2(i,j)=dandelions(i,j)-beta(i,j)*alpha*(dandelions_mean(1,j)-beta(i,j)*alpha*dandelions(i,j)); % eq.(13) in this paper
        end
    end
    dandelions=dandelions_2;
    % Check boundries
    dandelions = max(dandelions,LB);
    dandelions = min(dandelions,UB);
    
    %% Landing stage
    Step_length=levy(Popsize,Dim,1.5);
    Elite=repmat(Best_position,Popsize,1);
    for i=1:Popsize
        for j=1:Dim
            dandelions_3(i,j)=Elite(i,j)+Step_length(i,j)*alpha*(Elite(i,j)-dandelions(i,j)*(2*t/Maxiteration)); % eq.(15) in this paper
        end
    end
    dandelions=dandelions_3;
    % Check boundries
    dandelions = max(dandelions,LB);
    dandelions = min(dandelions,UB);
    
    %%
    % Calculated all dandelion seeds' fitness values
    for i=1:Popsize
%         dandelionsFitness(1,i)=Fobj(dandelions(i,:));
        dandelionsFitness(1,i)=testFunction(dandelions(i,:)', fhd, fNumber);
   
    end
    
    
    % Arrange dandelion seeds from good to bad according to fitness values
    [~,sorted_indexes]=sort(dandelionsFitness);
    dandelions=dandelions(sorted_indexes(1:Popsize),:);
    SortfitbestN = dandelionsFitness(sorted_indexes(1:Popsize));
    
    %Update the optimal dandelion seed
    if SortfitbestN(1)<Best_fitness
        Best_position=dandelions(1,:);
        Best_fitness=SortfitbestN(1);
    end
end
t=t+1;
% time = toc;
bestSolution= Best_position;
bestFitness=Best_fitness;
iteration=t;
 

end


% ___________________________________
function [z] = levy(n,m,beta)
% beta is set to 1.5 in this paper
num = gamma(1+beta)*sin(pi*beta/2);
den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma_u = (num/den)^(1/beta);
u = random('Normal',0,sigma_u,n,m);
v = random('Normal',0,1,n,m);
z =u./(abs(v).^(1/beta));
end


