
% Draw the benchmarck functions in 3D

function func_plot(func_name)

[lb,ub,dim,fobj]=Select_Function(func_name);

switch func_name 
    case 'F1' 
        x=-1:0.1:2; y=x;%-1:0.1:1;%x; %[-100,100]
    case 'F2' 
        x=-4.5:0.1:4.5; y=x; %[-10,10]    
    case 'F3' 
        x=-2*pi:0.1:2*pi; y=x; %[-100,100]
    case 'F4' 
        x=-100:0.1:100; y=x; %[-100,100]
    case 'F5' 
        x=-10:0.1:10; y=x; %[-5,5]
    case 'F6' 
        x=-5:0.1:10; y=0:0.1:15;%x; %[-5,5]
    case 'F7' 
        x=-5:0.1:15;  y=x;  %[-1,1]
    case 'F8' 
        x=-10:0.1:10;y=x; %[-500,500]
    case 'F9' 
        x=-15:0.1:-3;y=-3:0.1:3; %[-5,5]    
    case 'F10' 
        x=-5:0.1:5; y=x;%[-500,500]
    case 'F11' 
        x=-30:0.10:30; y=x;%[-0.5,0.5]
    case 'F12' 
        x=-20:0.1:20; y=x;%[-pi,pi]
    case 'F13' 
        x=-100:0.1:100; y=x;%[-3,1]
    case 'F14' 
        x=-10:0.1:10; y=x;%[-100,100]
    case 'F15' 
        x=-10:0.1:10; y=x;%[-5,5]
    case 'F16' 
        x=0:0.01:pi; y=x;%[-5,5]
    case 'F17' 
        x=-10:0.1:10; y=x;%[-5,5]
    case 'F18' 
        x=-100:0.1:100; y=x;%[-5,5]
    case 'F19' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F20' 
        x=0.9:0.1:9; y=x;%[-5,5]        
    case 'F21' 
        x=-1:0.1:1; y=x;%[-5,5]
    case 'F22' 
        x=-1:0.01:1; y=x;%[-5,5]     
    case 'F23' 
        x=0:0.1:pi; y=x;%[-5,5] 
    case 'F24' 
        x=-1:0.1:1;
        y=x;%[-5,5]  
     case 'F25' 
        x=-100:0.1:100; y=x;%[-5,5]  
    case 'F26'
        x=0:0.1:pi; y=x;%[-5,5]  
    case 'F27'
        x=-10:0.1:10; y=x;%[-5,5]  
     case 'F28'
        x=-32:0.1:32; y=x;%[-5,5]  
     case 'F29'
        x=-1:0.1:4; y=x;%[-5,5]  
    case 'F30'
        x=-5.12:0.1:5.12; y=x;%[-5,5]  
    case 'F31'
        x=-100:0.5:100; y=x;%[-5,5] 
     case 'F32'
        x=0:0.1:1;y=x;
      case 'F33'
        x=-1.28:0.1:1.28; y=x;%[-5,5]  
    case 'F34'
        x=-5.12:0.1:5.12; y=x;%[-5,5]  
    case 'F35'
        x=-30:0.1:30; y=x;%[-5,5]  
    case 'F36'
        x=-100:0.1:100; y=x;%[-5,5]  
    case 'F37'
        x=-100:0.1:100; y=x;%[-5,5] 
     case 'F38'
        x=-10:0.1:10; y=x;%[-5,5] 
end    

L=length(x);
f=[];

for i=1:L
    for j=1:L
        if strcmp(func_name,'F15')==0 && strcmp(func_name,'F19')==0 && strcmp(func_name,'F20')==0 && strcmp(func_name,'F21')==0 && strcmp(func_name,'F22')==0 && strcmp(func_name,'F23')==0
            f(i,j)=fobj([x(i),y(j)]);
        end
        if strcmp(func_name,'F15')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end
        if strcmp(func_name,'F19')==1
            f(i,j)=fobj([x(i),y(j),0]);
        end
        if strcmp(func_name,'F20')==1
            f(i,j)=fobj([x(i),y(j),0,0,0,0]);
        end       
        if strcmp(func_name,'F21')==1 || strcmp(func_name,'F22')==1 ||strcmp(func_name,'F23')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end          
    end
end

surfc(x,y,f,'LineStyle','none');

end

