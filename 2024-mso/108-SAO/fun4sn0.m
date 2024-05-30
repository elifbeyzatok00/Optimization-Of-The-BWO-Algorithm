function e=fun4sn0(b)

global u;

th=pi/180;
u=[1 2 35*th 70*th]; % Desired Vector

%%%%%%%%%%%%%%%
% Swapping vector b
%%%%%%%%%%%%%%%
[~, ix] = sort(u); 
[~, ix1(ix)] = sort(b);  
b = b(ix1);


[R,C]=size(b);
P=C/2;
M=2*C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate observed vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%
xo=zeros(1,M);
for k=1:M 
    for i=1:P 
        xo(1,k)=xo(1,k)+u(i)*exp(-1i*(k-1)*pi*cos(u(P+i)));
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Estimated vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xe=zeros(1,M);
for k=1:M
        for i=1:P
            xe(1,k)=xe(1,k)+b(i)*exp(-1i*(k-1)*pi*cos(b(P+i)));
        end 
end

%%%%%%%%%%%%%%%%%%%%
% Mean Square Error
%%%%%%%%%%%%%%%%%%%%
abc=0.0;
for m1=1:M
  abc=abc+(abs(xo(1,m1)-xe(1,m1))).^2; %This is square of the absolute of the difference of the two values i.e. |(observed-estimated)|^2
end 
abc=abc/M; %This is the mean squared error value

e=abc;



