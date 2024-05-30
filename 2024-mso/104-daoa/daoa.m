function [bestSolution, bestFitness, iteration]=daoa(fhd, dimension, maxIteration, fNumber)

config;

Dim=dimension;                          
LB = lbArray;       
UB = ubArray;           
N=5;                           
M_Iter=ceil(maxIteration/N);
Best_P=zeros(1,Dim);
Best_FF=inf;
%% Initialize the positions of solution
X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));          % (fitness values)
Ffun_new=zeros(1,size(Xnew,1));   % (fitness values)
C_Iter=1;
Mu=0.001;
alpha=25;
for i=1:size(X,1)   
    Ffun(1,i)=testFunction(X(i,:)', fhd, fNumber);
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
%% Main Loop
while C_Iter<M_Iter+1                        % Main loop
    
   %DAF=((C_Iter/M_Iter+1)^((-1)*alpha));    %DAF1
    DAF=(M_Iter+1/C_Iter)^((alpha));        %DAF2
    
    DCS=.99*((1-(C_Iter/M_Iter)^(0.5)));     % DCS
    
    
    %Upate the Position of solutions
    
    for i=1:size(X,1)                        % if each of the UB and LB has a just value
        for j=1:size(X,2)
            r1=rand();
            if (size(LB,2)==1)
                if r1<DAF
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=(Best_P(1,j)/(DCS+eps)*((UB-LB)*Mu+LB));
                    else
                        Xnew(i,j)=(Best_P(1,j)*DCS*((UB-LB)*Mu+LB));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=(Best_P(1,j)-DCS*((UB-LB)*Mu+LB));
                    else
                        Xnew(i,j)=(Best_P(1,j)+DCS*((UB-LB)*Mu+LB));
                    end
                end
            end
            
            
            if (size(LB,2)~=1)                          % if each of the UB and LB has more than one value
                r1=rand();
                if r1<DAF
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=((Best_P(1,j)/(DCS+eps)*((UB(j)-LB(j))*Mu+LB(j))));
                    else
                        Xnew(i,j)=((Best_P(1,j)*DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=((Best_P(1,j)-DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    else
                        Xnew(i,j)=((Best_P(1,j)+DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    end
                end
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB;                    % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB;                    % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
        
        Ffun_new(1,i)=testFunction(Xnew(i,:)', fhd, fNumber);
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
            Best_FF=Ffun(1,i);
            Best_P=X(i,:);
        end
        
    end
    
    C_Iter=C_Iter+1;  % incremental iteration
    
end
bestSolution=Best_P;
bestFitness=Best_FF;
iteration=1;

end

function X=initialization(N,Dim,UB,LB)
B_no= size(UB,2); % numnber of boundaries
if B_no==1
    X=rand(N,Dim).*(UB-LB)+LB;
end
% If each variable has a different lb and ub
if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
end
end
