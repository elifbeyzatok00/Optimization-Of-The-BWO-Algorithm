
function [lb,ub,dim,fobj] = Get_Functions_details(F)

switch F
    case 'F1' % Tension/compression spring design
        fobj = @F1;        lb=[0.05,0.25,2];        ub=[2,1.3,15];        dim=3;   
    case 'F2' %Pressure vessel design
        fobj = @F2;        lb=[0,0,10,10];        ub=[99,99,200,200];        dim=4;
    case 'F3' %Welded beam design
        fobj = @F3;        lb=[0.1,0.1,0.1,0.1];        ub=[2,10,10,2];        dim=4;
    case 'F4' % speed reducer problem
        fobj = @F4;  lb=[2.6,0.7,17,7.3,7.3,2.9,5.0];   ub=[3.6,0.8,28,8.3,8.3,3.9,5.5];  dim=7;
    case 'F5' % Gear train design problem
        fobj = @F5;  lb=[12,12,12,12]; ub=[60,60,60,60]; dim=4;
   
    case 'F6' % Three bar truss design
        fobj = @F6; lb=[0,0]; ub=[1,1]; dim=2;
   
end
end
% F1
function o = F1(x) % Tension/compression spring design
cost=(x(3)+2)*x(2)*(x(1)^2);
g(1)=1-((x(3)*(x(2)^3))/(71785*(x(1)^4)));   
gtmp=(4*x(2)^2-x(1)*x(2))/(12566*(x(2)*x(1)^3-x(1)^4));
g(2)=gtmp+1/(5108*x(1)^2)-1;                      
g(3)=1-((140.45*x(1))/((x(2)^2)*x(3)));        
g(4)=((x(1)+x(2))/1.5)-1;            
% Apply all inequality constraints as a penalty function 
lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
o=cost+ Z;
% Apply all equality constraints (when geq=[], length->0)
% for k=1:length(geq),
%    Z=Z+lam*geq(k)^2*geteqH(geq(k));
% end
end

% F2

function o = F2(x) % Pressure vessel design
cost=0.6224*x(1)*x(3)*x(4)+1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3);
g(1)=-x(1)+0.0193*x(3);                                  
g(2)=-x(2)+0.00954*x(3);                                  
g(3)=-pi*x(3)^2*x(4)-(4/3)*pi*x(3)^3+1296000;  
g(4)=x(4)-240;       
lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
o=cost+ Z;

end

% F3   Welded beam design
function o = F3(x)
cost=1.10471*(x(1)^2)*x(2)+0.04811*x(3)*x(4)*(14.0+x(2));
% Inequality constraints
Q=6000*(14+x(2)/2);
D=sqrt(x(2)^2/4+(x(1)+x(3))^2/4);
J=2*(x(1)*x(2)*sqrt(2)*(x(2)^2/12+(x(1)+x(3))^2/4));
alpha=6000/(sqrt(2)*x(1)*x(2));
beta=Q*D/J;
tau=sqrt(alpha^2+2*alpha*beta*x(2)/(2*D)+beta^2);
sigma=504000/(x(4)*x(3)^2);
delta=65856000/(30*10^6*x(4)*x(3)^3);
F=4.013*(30*10^6)/196*sqrt(x(3)^2*x(4)^6/36)*(1-x(3)*sqrt(30/48)/28);

g(1)=tau-13600;                  
g(2)=sigma-30000;             
g(3)=x(1)-x(4);                    
g(4)=0.10471*x(1)^2+0.04811*x(3)*x(4)*(14+x(2))-5.0;  
g(5)=0.125-x(1);                
g(6)=delta-0.25;                  
g(7)=6000-F;                    
lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
o=cost+ Z;
end

function o = F4(x) % speed reducer problem
cost=0.7854*x(1)*(x(2)^2)*(3.3333*(x(3)^2)+ 14.9334*x(3)- 43.0934)-1.508*x(1)*((x(6)^2)+(x(7)^2))...
    +7.4777*((x(6)^3)+(x(7)^3))+0.7854*(x(4)*(x(6)^2)+x(5)*(x(7)^2));
g(1)=(27/(x(1)*(x(2)^2)*x(3)))-1; 
g(2)=(397.5/(x(1)*(x(2)^2)*(x(3)^2)))-1; 
g(3)=(1.93*(x(4)^3)/(x(2)*(x(6)^4)*x(3)))-1;   
g(4)=(1.93*(x(5)^3)/(x(2)*(x(7)^4)*x(3)))-1; 
g(5)=(((((745*x(4)/(x(2)*x(3)))^2)+16.9*(10^6))^(1/2))/(110*(x(6)^3)))-1;
g(6)=(((((745*x(5)/(x(2)*x(3)))^2)+157.5*(10^6))^(1/2))/(85*(x(7)^3)))-1;
g(7)=(x(2)*x(3)/40)-1;
g(8)=((5*x(2))/x(1))-1;
g(9)=(x(1)/(12*x(2)))-1;
g(10)=((1.5*x(6)+1.9)/x(4))-1;
g(11)=((1.1*x(7)+1.9)/x(5))-1;
lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
o=cost+ Z;
end
function o = F5(x) % Gear train design problem
cost=((1/6.931)-((x(3)*x(2))/(x(1)*x(4))))^2;  
%lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
o=cost;
end

function o = F6(x) % %% Three-bar truss design problem
    cost= (2.*sqrt(2).*x(1)+x(2))*100;
    g(1) = (sqrt(2).*x(1)+x(2))./(sqrt(2).*x(1).^2+2.*x(1).*x(2))*2-2;
    g(2) = x(2)./(sqrt(2).*x(1).^2+2.*x(1).*x(2))*2-2;
    g(3) = 1./(sqrt(2).*x(2)+x(1))*2-2;
    lam=10^15; Z=0;     for k=1:length(g),    Z=Z+ lam*g(k)^2*getH(g(k));   end
    o=cost+ Z;
end

function H=getH(g)
if g<=0 
    H=0; 
else
    H=1; 
end
end
