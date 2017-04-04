function [corCalc, incCalc, omiCalc, preCalc, deltaFoverFcor, deltaFoverFinc, deltaFoverFomi, deltaFoverFpre, artefactTimes]=controlincluded_61hz_4(filename)



close all
load(filename); %edited 170306#
trigger=data(:,3); %Edited 170306
signal=data(:,2); %Edited 170306
% signal_lockin=data(:,1);

daqframe = 1/0.0004995; %Hz

lockframe = 61; %Hz


%% Start of Calcium signal analysis
pointsPerLED = ceil(daqframe/lockframe); %daqframes per lockframe
% plot(signal(1:10000))
% pause
% floor(numel(daqframe)/lockframe)


i=1; %will loop through signal
inVect=1;   %When not perfectly synced i will loop too fast,
%sometimes skipping a peak therefore inVect to store peak indices


while i<floor(numel(signal)/(2*pointsPerLED))-3 %Look for excitation peak per frame
    % in the photodetector output
    
    [m1 indx1]=max(signal((i-1)*(pointsPerLED*2)+1:i*(pointsPerLED*2))); %find index of peak
    [m2 indx2]=max(signal((i-1)*(pointsPerLED*2)+pointsPerLED:i*(pointsPerLED*2)+pointsPerLED));
    
    if indx1>pointsPerLED & indx2<pointsPerLED %if in second half of frame
        
        if m1>=m2 %take the maximum peak
            indxDef(inVect)=indx1+(i-1)*(pointsPerLED*2);
        else
            indxDef(inVect)=indx2+(i-1)*(pointsPerLED*2)+pointsPerLED;
        end
        
    elseif indx1<=pointsPerLED %If peak in first half of frame select index
        indxDef(inVect)=indx1+(i-1)*pointsPerLED*2;
    else
        indxDef(inVect)=indx2+(i-1)*pointsPerLED*2+pointsPerLED;
    end
    if indxDef>1
        if signal(indxDef(inVect))<signal(indxDef(inVect)-1) %if a peak got skipped over
            indxDef(inVect)=indxDef(inVect)-1;
        end
    end
    
    if inVect>1 & indxDef(inVect)>indxDef(inVect-1)+2*pointsPerLED %if the script misses a peak take it from indx1
        tempInd=indxDef(inVect); %temporary store this value
        indxDef(inVect)=indx1+(i-1)*(pointsPerLED*2);
        indxDef(inVect+1)=tempInd; %put in as last value
        inVect=inVect+2; %vector increases in size with 2
    else
        inVect=inVect+1; %vector increases in size with 1
    end
    i=i+1; %while counter
end

indDif=diff(indxDef);
is_same=find(indDif<2); %take out double detected peaks
indxDef(is_same+1)=[];
hold on
plot(indxDef(1:100),signal(indxDef(1:100)),'*r')
% pause
% for i=1:floor(numel(signal)/window)
%     [m ind]=max(signal((i-1)*window+1:i*window));
%     localpeaks(i)=ind+(i-1)*window;
% end
hold off
h470=[]; %Large peak currently this is calcium signal
l405=[]; %Isosbestic wavelength used for normalization/internal control



for i=1:numel(indxDef) % photodetector output
    h470(i)=mean(signal(indxDef(i)+5:indxDef(i)+floor(pointsPerLED/2)));
    l405(i)=mean(signal(indxDef(i)+pointsPerLED+5:indxDef(i)+pointsPerLED+floor(pointsPerLED/2)));
end


calAct=h470./l405;
size(calAct)

afname=[num2str(filename(1:end-4)) '_af.mat'];
allFiles=dir('*.mat');

for i=1:numel(allFiles)
    if strcmp(allFiles(i).name, afname)
        load(allFiles(i).name);
        artefactTimes=[startEnd(1:2:end) startEnd(2:2:end)];
        plot(calAct)
        axis([1 numel(calAct) 0 10])
        hold on
        plot(startEnd(1:end),calAct(startEnd(1:end)),'*r')
        break
    elseif i==numel(allFiles)
        plot(calAct)
        axis([1 numel(calAct) 0 10])
        artefactN=input('how many artefacts are there? ');
        [startEnd unused]=ginput(artefactN*2);
        startEnd=round(startEnd);
        artefactTimes=[startEnd(1:2:end) startEnd(2:2:end)];
        hold on
        plot(startEnd(1:end),calAct(startEnd(1:end)),'*r')
        startEnd
        artefacts=[num2str(filename(1:end-4)) '_af.mat'];
        save(artefacts, 'artefactN', 'startEnd');
        break
    end
end


if startEnd(end)>numel(calAct)
    startEnd(end)=numel(calAct); %make sure the calcium signal stays same length
end

for i=1:artefactN
    linLength=startEnd(i*2)-startEnd(i*2-1)+1;
    %     numel( calAct_noArtefact(startEnd(i*2-1):startEnd(i*2)))
    %     pause
    
    calAct(startEnd(i*2-1):startEnd(i*2))=linspace(calAct(startEnd(i*2-1)),calAct(startEnd(i*2)),linLength);
end
% plot(calAct)
% hold on
% plot(calAct_noArtefact,'r')
% 
% axis([1 numel(calAct) 0 10])
% pause

% x=[1:numel(calAct)]; %Linear fit to remove bleaching/slow drift
% P=polyfit(x,calAct,1);
% yfit=P(1)*x+P2;
% calAct=calAct-yfit;


% X = 1:100;
% Y = 3*exp(-0.5.*(X))+10;
% f = @(c,x) c(1)*exp(c(2)*x)+c(3);
% D = lsqcurvefit(f,zeros(3,1),X,Y)
% Z = f(D,X);
% plot(X,Z)
% grid on
% hold on
% plot(X,Y,['g','o'])

% X = 1:100;
% Y = 3*exp(-0.5.*(X))+10;
% f = @(c,x) c(1)*exp(c(2)*x)+c(3);
% D = lsqcurvefit(f,zeros(3,1),x,calAct)
% Z = f(D,x);

% grid on
% hold on
% plot(x,calAct,'r')
% plot(x,Z)
%
%
% pause




% cLPF=zeros(1,numel(calAct_noArtefact));

% cLPF=lowpassF(calAct,30.5,0.1);
% cLPF_lockin=lowpassF(calAct_lockin,30.5,0.1);

plot(calAct)
hold on


for i=1:artefactN %create lowpass signal to subtract 
    
%     if artefactN == 1
%         cLPF(1:startEnd(i*2-1))=lowpassF(calAct_noArtefact(1:startEnd(i*2-1)),30.5,0.05); %before first artefact
%         cLPF(startEnd(i*2):end)=lowpassF(calAct_noArtefact(startEnd(i*2):end),30.5,0.05); %after last artefact
%     
%     elseif artefactN>1
%         if i==1
%             cLPF(1:startEnd(i*2-1))=lowpassF(calAct_noArtefact(1:startEnd(i*2-1)),30.5,0.05)%before first artefact
%             cLPF(startEnd(i*2):startEnd(i*2+1))=lowpassF(calAct_noArtefact(startEnd(i*2):startEnd(i*2+1)),30.5,0.05)
%         elseif i==artefactN
%             cLPF(startEnd(i*2):end)=lowpassF(calAct_noArtefact(startEnd(i*2):end),30.5,0.05) %after last artefact
%         else
%             cLPF(startEnd(i*2):startEnd(i*2+1))=lowpassF(calAct_noArtefact(startEnd(i*2):startEnd(i*2+1)),30.5,0.05);%between artefacts
% 
%         end
% 
%     end
   
    linLength=startEnd(i*2)-startEnd(i*2-1)+1;
    %     numel( calAct_noArtefact(startEnd(i*2-1):startEnd(i*2)))
    %     pause
    calAct(startEnd(i*2-1):startEnd(i*2))=linspace(calAct(startEnd(i*2-1)),calAct(startEnd(i*2)),linLength);
    %Interpolate during artefacts

end

%% matching adjusted trace to MPC epochs
% finds exact timepoints of MPC epochs in the fluorescence acquisition
trig=get_triggerTimes(trigger);


cd ..
cd 'bhv'

lst=dir('*.txt');
for i=1:numel(lst)
    if filename(1:19)==lst(i).name(1:19)
        newFilename=lst(i).name;
    end
end

[perf]=read_5choice(newFilename);

cd ..
cd 'calcium'

i=1;

if strcmp(filename,'!2017-03-15_15h34m.SubjectGMDM1-LEFT.mat')
    perf(48:end)=[];
end


%MED does not give new trigger after premature response. Following loop
%deletes the premature trials from the behavioral data
while i<numel(perf)
    
    
    for j=1:artefactN
        if trig(i)/pointsPerLED-15>artefactTimes(j,1)/30.5 & trig(i)/pointsPerLED+15<artefactTimes(j,2)/30.5
            trig(i)=[];
            perf(i)=[];
        end
    end
    i=i+1;
end

errors=0; %in case of no problems with number of started trials and triggers

%Sometimes animal starts trial just before being taken out of 5choice. This
%leads to last trial being ommission. Delete this trial

% perf(1)=[];
% trig(1)=[];
if perf(end)==4
    trig(end)=[];
    perf(end)=[];
    errors=2; %error at the end
end

%In rare cases that the calcium acquisition is started after the rat has
%started performig the 5choice, delete the first trial, since it does not
%have a trigger or calcium signal (definition might need improvement)
if numel(trig)<numel(perf)
    perf(1)=[];
    if errors==2
        errors=3; %error both at start and end
    else
        errors=1; %error at start of protocol
    end
end
%%Until here behavior

trigstart=trig/daqframe;
% creates fluorescence traces for every trial
traceWindow = 5;
totalWindow = traceWindow*pointsPerLED*2;

for i=1:length(trigstart)-1 %Fluorescent signal during trials
    
    %     capitalF=mean(calAct(trigstart(i):trigstart(i)-100));%what is the fluorescence in 1 s before trial start
    %     i
    
    %     size(calAct)
    % %     trigstart(i)
    % %     trigstart(i)+1000
    %    trigstart(i)-61
    %    trigstart(i) + 10*61
    %    numel(calAct)
    %    trigstart(end-10:end)
    %
    adjustedFluorescence=calAct(round((trigstart(i)-traceWindow)*30.5):round((trigstart(i)-1)*30.5)+336)'; %Delta F
%     pause
    CalAct_trigstart(i,:)= adjustedFluorescence;
%     pause
    %dF over F
    Fzero = mean(CalAct_trigstart(i,(traceWindow*pointsPerLED-2*pointsPerLED:traceWindow*pointsPerLED)));
    F = CalAct_trigstart(i,:);
    deltaFoverF(i,:)=100*((F-Fzero)/Fzero);
    %
    
    %     adjustedFluorescence_lockin=calAct_lockin(round((trigstart(i)-1)*30.5):round((trigstart(i)-1)*30.5)+336)'; %Delta F
    %     CalAct_trigstart_lockin(i,:)= adjustedFluorescence_lockin;
    
end

for i=1:length(trigstart)-1 %Assigns traces to trial outcomes
    cor=find(perf==1);
    inc=find(perf==2);
    omi=find(perf==4);
    pre=find(perf==3);
    if numel(cor)>0
        corCalc=CalAct_trigstart(cor,:);
        deltaFoverFcor=deltaFoverF(cor,:);
    else
        corCalc=nan(1,ceil(traceWindow*30.5));
        deltaFoverFcor=nan(1,ceil(traceWindow*30.5));
    end
    
    if numel(inc)>0
        incCalc=CalAct_trigstart(inc,:);
        deltaFoverFinc=deltaFoverF(inc,:);
    else
        incCalc=nan(1,459);
        deltaFoverFinc=nan(1,459);
    end
    
    if numel(omi)>0
        omiCalc=CalAct_trigstart(omi,:);
        deltaFoverFomi=deltaFoverF(omi,:);
    else
        omiCalc=nan(1,459);
        deltaFoverFomi=nan(1,459);
    end
    
    if numel(pre)>0
        preCalc=CalAct_trigstart(pre,:);
        deltaFoverFpre=deltaFoverF(pre,:);
    else
        preCalc=nan(1,459);
        deltaFoverFpre=nan(1,459);
    end
end

%% plot outcomes
output(1) = figure; % plot adjusted fluorescence readout
plot([-traceWindow-1/30.5:1/30.5:10],mean(corCalc),'b')
hold on
plot([-traceWindow-1/30.5:1/30.5:10],mean(omiCalc),'k')
if numel(inc)>0
    plot([-traceWindow-1/30.5:1/30.5:10],mean(incCalc),'r')
end 
yL = get(gca,'YLim');
line([0 0],yL,'Color','m');
line([5 5],yL,'Color','g');
title('GCaMP6m fluorescence during 5-choice trials')%title
xlabel('time(s)') % x-axis label
ylabel('dF470/F405') % y-axis label
legend('Correct responses', 'Omissions', 'Incorrect responses')

output (2) = figure; % plot dF/F for adjusted fluorescence
plot([-traceWindow-1/30.5:1/30.5:10],mean(deltaFoverFcor))
hold on

title('GCaMP6m fluorescence during 5-choice trials')%title
xlabel('time(s)') % x-axis label
ylabel('dF/F') % y-axis label
plot([-traceWindow-1/30.5:1/30.5:10],mean(deltaFoverFomi),'k')
if numel(inc)>-1
    plot([-traceWindow-1/30.5:1/30.5:10],mean(deltaFoverFinc),'r')
end
legend('Correct responses', 'Omissions', 'Incorrect responses')
yL = get(gca,'YLim'); 
line([0 0],yL,'Color','m');
line([5 5],yL,'Color','g');


NewFolder = ([num2str(filename(20:end-4)) '_output']); %Creates new folder to save figures in
mkdir (NewFolder)
cd (NewFolder)
savefig (output,'two_plots.fig') 
cd .. 


end