clc;
clear;


format long;
Run_No=1;
fhd=@cec22_test_func;


D=20;
for fun=1:12
    eval(['load input_data/shift_data_' num2str(fun) '.txt']);
    eval(['O=shift_data_' num2str(fun) '(1,1:D);']);
    O = rand(1,10);
    f(fun,:)= feval(fhd, O', fun);
end
disp(f)
    


