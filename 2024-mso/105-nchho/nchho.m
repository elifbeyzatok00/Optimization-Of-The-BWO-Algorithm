function [bestSolution, bestFitness, iteration]=nchho(fhd, dimension, maxIteration, fNumber)

config;
N=100;
T=ceil(maxIteration/(N*2.5));
lb=lbArray;
ub=ubArray;
dim=dimension;
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=0;
%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
CNVG=zeros(1,T);
t=0; % Loop counter
while t<T
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness=testFunction(X(i,:)', fhd, fNumber);
        % Update the location of Rabbit
        if fitness>Rabbit_Energy
            Rabbit_Energy=fitness;
            Rabbit_Location=X(i,:);
        end
    end
    
    E1=abs(2*(1-(t/T))-2); % factor to show the decreaing energy of rabbit
    a1 = 4;              % Initial chaotic map parameter configuration
    teta = 0.7;         % Initial chaotic map parameter configuration
    % Update the location of Harris' hawks
    for i=1:size(X,1)
         for ii=1:4
          Cm(1,ii) = abs((a1/4)*sin(pi*teta));
          teta = Cm(1,ii);
        end
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if q<0.5
                % perch based on other family members
                 X(i,:)=X_rand-Cm(1,1)*abs(X_rand-2*Cm(1,2)*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
               X(i,:)=(Rabbit_Location(1,:)-mean(X))-Cm(1,3)*((ub-lb)*Cm(1,4)+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                w1=2*exp(-(8*t/T)^2);         % Non-linear control Parameter
                Jump_strength=2*(1-rand());
                X1=w1*Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                
                if testFunction(X1', fhd, fNumber) > testFunction(X(i,:)', fhd, fNumber)
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=w1*Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                    if testFunction(X2', fhd, fNumber) > testFunction(X(i,:)', fhd, fNumber)
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                w1=2*exp(-(8*t/T)^2);
                Jump_strength=2*(1-rand());
                X1=w1*Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                
                if testFunction(X1', fhd, fNumber) > testFunction(X(i,:)', fhd, fNumber)
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=w1*Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                    if testFunction(X2', fhd, fNumber) > testFunction(X(i,:)', fhd, fNumber)
                        X(i,:)=X2;
                    end
                end
            end
            %%
        end
    end
    t=t+1;
    CNVG(t)=Rabbit_Energy;
end
bestSolution=Rabbit_Location;
bestFitness=Rabbit_Energy;
iteration=1;
end
% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

function [X]=initialization(N,dim,up,down)
if size(up,2)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,2)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
end