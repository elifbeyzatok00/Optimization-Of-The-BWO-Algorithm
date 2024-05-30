clear;
clc;
algorithms = {'bwo', 'bwo_case1', 'bwo_case2', 'bwo_case3', 'bwo_case4', 'bwo_case5', 'bwo_case6', 'bwo_case7'};
algorithmsNumber = length(algorithms); 
functionsNumber = 12; experimentNumber = 1; 
experimentsName = {'cec2022'};
filename = 'test4.xlsx'; run = 21;
solution = zeros(algorithmsNumber * experimentNumber, functionsNumber, run);
solutionR = zeros(algorithmsNumber, functionsNumber * experimentNumber, run);
friedmanR = zeros(functionsNumber * experimentNumber, algorithmsNumber);
result = zeros(functionsNumber  * algorithmsNumber, 6 * experimentNumber);
firstWilcon = zeros(1, run);
secondWilcon = zeros(1, run);
sumWilcon = zeros(algorithmsNumber*experimentNumber, 3);

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

for i = 1 : functionsNumber
    for k=1:experimentNumber
        for m=1:algorithmsNumber
            result(algorithmsNumber*(i-1)+m,6*k-5) = min(solution(algorithmsNumber*(k-1)+m,i,:));
            result(algorithmsNumber*(i-1)+m,6*k-4) = mean(solution(algorithmsNumber*(k-1)+m,i,:));
            result(algorithmsNumber*(i-1)+m,6*k-3) = std(solution(algorithmsNumber*(k-1)+m,i,:));
            result(algorithmsNumber*(i-1)+m,6*k-2) = median(solution(algorithmsNumber*(k-1)+m,i,:));
            result(algorithmsNumber*(i-1)+m,6*k-1) = max(solution(algorithmsNumber*(k-1)+m,i,:));
        end
        for m=1:algorithmsNumber-1
            for n=1:run
                firstWilcon(n) = solution(algorithmsNumber*(k-1)+1,i,n);
                secondWilcon(n) = solution(algorithmsNumber*(k-1)+m+1,i,n);
            end
            [p,h, stats] = ranksum(firstWilcon, secondWilcon);
            wilcon = 0;
            disp(stats);
            if(h)
                if(p < 0.05)
                    if(stats.zval < 0)
                        wilcon = 2; % firstWilcon daha iyi
                    else
                        wilcon = 1; % secondWilcon daha iyi
                    end  
                end
            end
            result(algorithmsNumber*(i-1)+m+1,6*k) = wilcon;
        end
    end
end

for i=1:experimentNumber
    for j=1:algorithmsNumber
        pW = 0; eW = 0; nW = 0;
        for k=1:functionsNumber
            vW = result(algorithmsNumber * (k - 1) + j, 6 * i);
            if(vW == 0)
                eW = eW + 1;
            elseif(vW == 1)
                pW = pW + 1;
            else
                nW = nW + 1;
            end         
        end
        sumWilcon(algorithmsNumber * (i - 1) + j, 1) = pW;
        sumWilcon(algorithmsNumber * (i - 1) + j, 2) = eW;
        sumWilcon(algorithmsNumber * (i - 1) + j, 3) = nW;
    end
end

for i = 1 : functionsNumber
    for k=1:experimentNumber
        compareResult = zeros(run, algorithmsNumber);
        for m=1:algorithmsNumber
            compareResult(:, m) = solution(algorithmsNumber*(k-1)+m,i,:);
        end
        [~,~,stats] = friedman(compareResult, 1, 'off');
        friedmanR((k - 1) * functionsNumber + i, :) = stats.meanranks;
    end
end

fridmanRes = zeros(experimentNumber, functionsNumber + 1, algorithmsNumber);
for i = 1 : functionsNumber
    for j = 1 : experimentNumber
        fridmanRes(j, i, :) = friedmanR(functionsNumber * (j - 1) + i, :);
    end
end
for i = 1 : algorithmsNumber
    for j = 1 : experimentNumber
        fridmanRes(j, functionsNumber + 1, i) = mean(fridmanRes(j, 1 : functionsNumber, i));
    end
end

fridmanResult = []; wStatistic = []; fStatistic = []; wStatisticS = {};
for i = 1 : experimentNumber
    fridmanResult = [fridmanResult, squeeze(fridmanRes(i,:,:))];
    wStatistic = [wStatistic, sumWilcon((i-1)*algorithmsNumber+1:i*algorithmsNumber,:)];
    fStatistic = [fStatistic; fridmanRes(i,functionsNumber+1,:)];
    for j = 1 : algorithmsNumber
        wStatisticS(j, i) = cellstr(strcat(int2str(wStatistic(j, 3 * (i - 1) + 1)), '/', int2str(wStatistic(j, 3 * (i - 1) + 2)), '/', int2str(wStatistic(j, 3 * i))));
    end
end

fStatistic = squeeze(fStatistic);
fSortIndex = 1;
if(experimentNumber == 1)
    fStatistic = fStatistic';
else
    fSortIndex = fSortIndex + experimentNumber;
    fStatistic = [fStatistic; mean(fStatistic)];
end
algoritmNames = upper(algorithms);
[~,fI] = sort(fStatistic(fSortIndex,:));
sortedAlgorithm = algoritmNames(fI);

xlswrite(filename, result, 'Wilcoxon');
xlswrite(filename, fridmanResult, 'Friedman');
xlswrite(filename, algoritmNames', 'Wstatistic', 'A2');
xlswrite(filename, experimentsName(1:experimentNumber), 'Wstatistic', 'B1');
xlswrite(filename, wStatisticS, 'Wstatistic', 'B2');
xlswrite(filename, algoritmNames, 'Fstatistic', 'B1');
xlswrite(filename, experimentsName', 'Fstatistic', 'A2');
xlswrite(filename, sortedAlgorithm, 'Fstatistic', 'B7');
xlswrite(filename, fStatistic, 'Fstatistic', 'B2');
disp('Bitti :)');