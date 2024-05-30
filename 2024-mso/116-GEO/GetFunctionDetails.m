
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Golden Eagle Optimizer (GEO) source codes version 1.0
%
%  Developed in:	MATLAB 9.6 (R2019a)
%
%  Programmer:		Abdolkarim Mohammadi-Balani
%
%  Original paper:	Abdolkarim Mohammadi-Balani, Mahmoud Dehghan Nayeri,
%					Adel Azar, Mohammadreza Taghizadeh-Yazdi,
%					Golden Eagle Optimizer: A nature-inspired
%					metaheuristic algorithm, Computers & Industrial Engineering.
%
%                  https://doi.org/10.1016/j.cie.2020.107050
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fun,nvars,lb,ub] = ...
	GetFunctionDetails (FunctionNumber)

switch FunctionNumber
	case 1
		fun   = @ F1;
		nvars = 2;
		lb    = -4.5 .* ones (1,nvars);
		ub    =  4.5 .* ones (1,nvars);
	case 2
		fun   = @ F2;
		nvars = 2;
		lb    = -10 .* ones (1,nvars);
		ub    =  10 .* ones (1,nvars);
	case 3
		fun   = @ F3;
		nvars = 2;
		lb    = -5 .* ones (1,nvars);
		ub    =  5 .* ones (1,nvars);
	case 4
		fun   = @ F4;
		nvars = 30;
		lb    = -1 .* ones (1,nvars);
		ub    =  1 .* ones (1,nvars);
	case 5
		fun   = @ F5;
		nvars = 30;
		lb    = -5 .* ones (1,nvars);
		ub    =  5 .* ones (1,nvars);
	case 6
		fun   = @ F6;
		nvars = 30;
		lb    = -100 .* ones (1,nvars);
		ub    =  100 .* ones (1,nvars);
	case 7
		fun   = @ F7;
		nvars = 30;
		lb    = -5.12 .* ones (1,nvars);
		ub    =  5.12 .* ones (1,nvars);
	case 8
		fun   = @ F8;
		nvars = 2;
		lb    = -5.2 .* ones (1,nvars);
		ub    =  5.2 .* ones (1,nvars);
	case 9
		fun   = @ F9;
		nvars = 2;
		lb    = -512 .* ones (1,nvars);
		ub    =  512 .* ones (1,nvars);
	case 10
		fun   = @ F10;
		nvars = 2;
		lb    = -5 .* ones (1,nvars);
		ub    =  5 .* ones (1,nvars);
	case 11
		fun   = @ F11;
		nvars = 2;
		lb    = -10 .* ones (1,nvars);
		ub    =  10 .* ones (1,nvars);
	case 12
		fun   = @ F12;
		nvars = 30;
		lb    = -32 .* ones (1,nvars);
		ub    =  32 .* ones (1,nvars);
	case 13
		fun   = @ F13;
		nvars = 30;
		lb    = -600 .* ones (1,nvars);
		ub    =  600 .* ones (1,nvars);
	case 14
		fun   = @ F14;
		nvars = 30;
		lb    = -2 .* ones (1,nvars);
		ub    =  2 .* ones (1,nvars);
	case 15
		fun   = @ F15;
		nvars = 10;
		lb    = zeros (1,nvars);
		ub    = pi .* ones (1,nvars);
	case 16
		fun   = @ F16;
		nvars = 30;
		lb    = -50 .* ones (1,nvars);
		ub    =  50 .* ones (1,nvars);
	case 17
		fun   = @ F17;
		nvars = 30;
		lb    = -50 .* ones (1,nvars);
		ub    = 50 .* ones (1,nvars);
	case 18
		fun   = @ F18;
		nvars = 30;
		lb    = -50 .* ones (1,nvars);
		ub    =  50 .* ones (1,nvars);
	case 19
		fun   = @ F19;
		nvars = 30;
		lb    = -500 .* ones (1,nvars);
		ub    =  500 .* ones (1,nvars);
	case 20
		fun   = @ F20;
		nvars = 30;
		lb    = -5.12 .* ones (1,nvars);
		ub    =  5.12 .* ones (1,nvars);
	case 21
		fun   = @ F21;
		nvars = 30;
		lb    =  -5 .* ones (1,nvars);
		ub    =  10 .* ones (1,nvars);
	case 22
		fun   = @ F22;
		nvars = 30;
		lb    = -100 .* ones (1,nvars);
		ub    =  100 .* ones (1,nvars);
	case 23
		fun   = @ F23;
		nvars = 30;
		lb    = -10 .* ones (1,nvars);
		ub    =  10 .* ones (1,nvars);
	otherwise
		fun   = @ F1;
		nvars = 2;
		lb    = -4.5 .* ones (1,nvars);
		ub    =  4.5 .* ones (1,nvars);
end

end

function f = F1 (x)
f1 = (1.5   - x(:,1) + x(:,1) .* x(:,2)   ).^2;
f2 = (2.25  - x(:,1) + x(:,1) .* x(:,2).^2).^2;
f3 = (2.625 - x(:,1) + x(:,1) .* x(:,2).^3).^2;
f =  f1 + f2 + f3;
end

function f = F2 (x)
f =   0.26 .* sum(x.^2, 2) ...
	- 0.48 .* prod(x, 2);
end

function f = F3 (x)
f =   2 .* x(:,1).^2 ...
	- 1.05 .* x(:,1).^4 ...
	+ (x(:,1).^6) ./ 6 ...
	+ x(:,1) .* x(:,2) ...
	+ x(:,2).^2;
end

function f = F4 (x)
f = -exp(-0.5 .* sum(x.^2, 2));
end

function f = F5 (x)
d = 2;
alpha = 0.1;
f = x(:,1) + d .* (sum(x(:,2:end).^2, 2).^alpha);
end

function f = F6 (x)
f = sum(x.^2, 2);
end

function f = F7 (x)
f = sum((x+0.5).^2, 2);
end

function f = F8 (x)
f = -  (1 + cos(12 .* sqrt(sum(x.^2, 2)))) ...
	./ (0.5 .* sum(x.^2, 2) + 2);
end

function f = F9 (x)
f = sum(-(x(:,2)+47).*sin(sqrt(abs(x(:,2)+x(:,1)/2+47)))-x(:,1).*sin(sqrt(abs(x(:,1)-(x(:,2)+47)))), 2);
end

function f = F10 (x)
f =   (x(:,1).^2 + x(:,2) - 11) .^ 2 ...
	+ (x(:,1) + x(:,2).^2 - 7) .^ 2;
end

function f = F11 (x)
f1 = (sin(3.*pi.*x(:,1))).^2;
f2 = ((x(:,1)-1).^2) .* (1 + (sin(3.*pi.*x(:,2))).^2);
f3 = ((x(:,2)-1).^2) .* (1 + (sin(2.*pi.*x(:,2))).^2);
f = f1 + f2 + f3;
end

function f = F12 (x)
nvars = size(x, 2);
f1 = exp(-0.2 .* sqrt(sum(x.^2,2)./nvars));
f2 = exp(sum(cos(2.*pi.*x),2)/nvars);
f = -20 .* f1 - f2 + 20 + exp(1);
end

function f = F13 (x)
nvars = size(x, 2);
t = 1:nvars;
f1 = sum((x-100).^2, 2);
f2 = prod(cos((x-100)./sqrt(t)), 2);
f = (1/4000) .* f1 - f2 + 1;
end

function f = F14 (x)
nvars = size(x, 2);
alpha = 1/8;
v = vecnorm(x,2,2).^2;
f = ((v-nvars).^2) .^ alpha ...
	+ (1/nvars) .* (0.5.*v + sum(x,2)) ...
	+ 0.5;
end

function f = F15 (x)
nvars = size(x, 2);
f1 = (1:nvars) .* (x.^2) / pi ;
f = - sum(sin(x) .* sin(f1).^20, 2);
end

function f = F16 (x)
a = 10;
k = 100;
m = 4;
nvars = size(x,2);
f1 = 10 .* (sin(pi.*funY(x(:,1)))).^2;
f2 = sum ((funY(x(:,1:end-1))-1).^2 .* (1+10.*sin(pi.*funY(x(:,2:end))).^2), 2);
f3 = (funY(x(:,end))-1) .^2 ;
f = (pi./nvars) .* (f1+f2+f3) + sum(funU(x,a,k,m), 2);
end

function f = F17 (x)
a = 5;
k = 100;
m = 4;
f1 = sin(3.*pi.*x(:,1)) .^ 2;
f2 = sum(((x(:,1:end-1)-1).^2) .* (1+sin(3.*pi.*x(:,2:end)).^2), 2);
f3 = ((x(:,end)-1).^2) .* (1+sin(2.*pi.*x(:,end)).^2);
f = 0.1 .* (f1+f2+f3) + sum(funU(x,a,k,m), 2);
end

function f = F18 (x)
f =   1 ...
	+ sum(sin(x).^2, 2) ...
	- 0.1 .* exp(-sum(x.^2, 2));
end

function f = F19 (x)
nvars = size(x, 2);
f = sum((x.^2 - (1:nvars)).^2,2);
end

function f = F20 (x)
f = sum(x.^2 - 10 .* cos(2.*pi.*x) + 10, 2);
end

function f = F21 (x)
f1 = zeros(size(x));
for i1 = 1:size(x, 2)-1
	f1(:,i1) = 100 .* (x(:,i1+1)-x(:,i1).^2).^2 + (1-x(:,i1)).^2;
end
f = sum(f1, 2);
end

function f = F22 (x)
f = 1 - cos(2 .* pi .* sqrt(sum(x.^2,2))) + 0.1 .* sqrt(sum(x.^2, 2));
end

function f = F23 (x)
f =    (sum((sin(x)).^2, 2) - exp(-sum(x.^2, 2))) ...
	.* exp(-sum(sin(sqrt(abs(x))).^2, 2));
end

function u = funU (x,a,k,m)
u = k .* ((x-a).^m) .* (x>a) + k.*((-x-a).^m) .* (x<(-a));
end

function y = funY (x)
y = 1 + (x+1) ./ 4;
end
