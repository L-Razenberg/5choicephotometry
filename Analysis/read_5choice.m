function [performance respT magLat]=read_5choice(filename)
%function read_5choice(filename) Will import the data from filename and
%perform simple analysis on the data stored on row 41 downwards

M=importdata(filename,' ',40); %import the data, consider line 1-39 as text
realdata=M.data;               %realdata now stores the numbers
endFind=mean(realdata,2);

yesTrial=find(endFind>0);
lastTrial=ceil(yesTrial(end)/12);  %The total number of trials

clear endFind yesTrial

startT=realdata(1:12:lastTrial*12,1); %starttime
respT=zeros(lastTrial,1); %response times
laser=zeros(lastTrial,1); %laser on =1
corct=zeros(lastTrial,1); %Correct response =1
incr=zeros(lastTrial,1); %Incorrect response =1
prem=zeros(lastTrial,1); %premature response >0
omis=zeros(lastTrial,1); %Omission trial =1
magLat=zeros(lastTrial,1); %magazine latency
cueLoc=zeros(lastTrial,1); %Cue location
compCor=zeros(lastTrial,1); %Compulsive responses correct location
compIncor=zeros(lastTrial,1);%Compulsive responses incorrect location

latPostTO=zeros(lastTrial,1); %not sure what TO is



for i=1:lastTrial
    cueLoc(i)=realdata((i-1)*12+1,3);
    latPostTO(i)=realdata((i-1)*12+1,4);
    
    
%     if realdata(i*12,1)~0
%         laser(i)=1;     %Was the laser on?
%     end
    if sum(realdata((i-1)*12+3,:))>0 %Correct trial
        corct(i)=1;

        respT(i)=sum(realdata((i-1)*12+3,:))-sum(realdata((i-1)*12+8,:));
        magLat(i)=realdata((i-1)*12+1,5);

    end
    if sum(realdata((i-1)*12+4,:))>0 %Incorrect trial
        incr(i)=1;
%         respT(i)=sum(realdata((i-1)*12+4),:)-sum(realdata((i-1)*12+6),:);
        magLat(i)=realdata((i-1)*12+1,5);
    end
    if sum(realdata((i-1)*12+9,:))>0
        prem(i)=sum(realdata((i-1)*12+9,:)); %will store number of premature responses
        if prem(i)>0
            prem(i)=1;
        end
%         respT(i)=sum((realdata((i-1)*12+2,:))/sum(realdata((i-1)*12+9,:)); %avg premature time
        
    end
    
    om=prem(i)+corct(i)+incr(i); %Will be 0 if no premature,correct or incorrect (==omission)
    if om==0
        omis(i)=1;
    end
    
%     if sum(realdata((i-1)*12+10,:))>0
%         compCor(i)=sum(realdata((i-1)*12+10,:));
%     end
%     if sum(realdata((i-1)*12+11,:))>1
%         compIncor=sum(realdata((i-1)*12+11,:));
%     end
%     
end


performance=corct+incr*2+prem*3+omis*4;


end

