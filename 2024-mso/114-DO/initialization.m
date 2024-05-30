function Positions=initialization(Popsize,Dim,UB,LB)

Boundary_no= size(UB,2); % numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% number for both UB and LB
if Boundary_no==1
     Positions=rand(Popsize,Dim).*(UB-LB)+LB;
end

% If each variable has a different LB and UB
if Boundary_no>1
    for i=1:Dim
        ub_i=UB(i);
        lb_i=LB(i);
        Positions(:,i)=rand(Popsize,1).*(ub_i-lb_i)+lb_i;      
    end
end