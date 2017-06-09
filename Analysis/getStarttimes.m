function [c times] = getStarttimes
% This function gets the starttimes from each new stage 

fileName = uigetfile % choose a text file

[a, b, struct]= CC_analysis_pstage_Lotte_Jeroen(fileName);
% save the struct of stages from the medpc output file

stages = 1:size(struct,2);

times = zeros(size(stages,2),1);

for i = stages
    times(i,1) = struct(i).data(1);
    
end


c = zeros(size(stages,2),1);
for n = 1:(size(stages,2)-1)
    if n+1 ~ 0;
    c(n,1) = times(n+1)-times(n);
    else
        c(n,1)=0;
    end 
end 
end 