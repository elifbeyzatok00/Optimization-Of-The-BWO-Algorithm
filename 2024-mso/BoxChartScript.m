clear;
clc;
algorithms = {'sos', 'fdb_sos'};
filename = 'test.xlsx'; run = 5;
algorithmsNumber = length(algorithms); functionsNumber = 12; experimentNumber = 1; 
solutionR = zeros(algorithmsNumber, functionsNumber * experimentNumber, run);
solution = zeros(algorithmsNumber * experimentNumber, functionsNumber, run);
for i = 1 : algorithmsNumber
    solutionR(i,:,:) = xlsread(filename, char(algorithms(i)));
    solutionR(i,6,:) = zeros(1, run);
end

for i = 1 : functionsNumber
    for j = 1 : run
        for k = 1 : algorithmsNumber * experimentNumber
            m = mod(k, algorithmsNumber); if(m == 0), m = algorithmsNumber; end
            n = (i - 1) * experimentNumber + ceil(k / algorithmsNumber);
            solution(k, i, j) = solutionR(m, n, j);
        end
    end
end

fn = 5;
en = 2;
values = [];
for m = 1  :algorithmsNumber
    values = [values ;  squeeze(solution(algorithmsNumber * (en - 1) + m, fn, :))'];
end
boxplot(values', {'AGDE','Case-1','Case-2','Case-3'});
xlabel('Algorithms');
ylabel('Function Error Value');

% for i = 1 : functionsNumber
%     for k=1:experimentNumber
%         values = [];
%         for m=1:algorithmsNumber
%             values = [values ;  squeeze(solution(algorithmsNumber*(k-1)+m,i,:))'];
%         end
%         boxplot(values', {'SFS','Case-1','Case-2','Case-3'});
%         xlabel('Algorithms');
%         ylabel('Function Error Value');
%     end
% end