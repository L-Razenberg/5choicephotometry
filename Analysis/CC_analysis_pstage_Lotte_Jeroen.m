function [data stageUpIndex s] = CC_analysis_pstage_Lotte_Jeroen(fileName)

% Gives output struct of MED-PC output per stage to be further analysed.
%
% 1. Place this script and MED-PC output script in same folder.
% 2. Add folder to Matlab path.
% 3. Change MED-PC output 'filename.MPC' to 'filename.txt'
% 4. run 'CC_analysis_pstage_Lotte_Jeroen(fileName)' in Matlab
% Command Window. Where 'fileName' = 'filename.txt'
%

% Find where x array starts and imports the whole X array.
[data,delimiterOut,headerlinesOut]=importdata(fileName,'',60);
idx = (find(ismember(data(:,1), 'X:')))+60-headerlinesOut;
clear data delimiterOut headerlinesOut
data=importdata(fileName,' ',idx);

% Cut the X array off and the last trial.
endFind = [1: 5: length(data.data)];

for basta = endFind
    if data.data(basta,1) == 0
        final = basta;
        data.data = data.data(1:final,:);
        break
    end
end

data = data.data;

clearvars -except data stage

% Divide in matrix per stage.
stageData = [];

indexing = 1:5:length(data)-1;
stageUpIndex = 1;
n = 2;
for i = indexing
    if data(i+5,3) - data(i,3) > 0
        stageUpIndex(n) = i+5;
        n = n+1;
    end
    
end
stages = size(stageUpIndex,2);

% create a struct
for x = 1:size(stageUpIndex,2)
    nr = num2str(x);
    s(x).stage = nr;
end 

%s = struct('name',{'1','2','3','4','5','6','7','8','9','10','11','12','13'},'data',[])

n=1;
for i = 1:stages
stageUpIndex(i);
    if (i<stages) 
output = data(stageUpIndex(i):stageUpIndex(i+1)-1,:);
else
    output = data(stageUpIndex(i):end,1:end);
    end
s(n).data=output;
n=n+1;
end

