function [cor, incor, omiss, pre, dfcor, dfincor, dfomiss, dfpre, allFiles]=loop_through_folder

%% file selection
% go into the folder where you stored the files
folder = uigetdir
cd (folder)
list = dir()

%select the desired matfiles, use shift or ctrl to select multiple files 
[filename, pathname, filterindex] = uigetfile( ...
{  '*.mat','MAT-files (*.mat)'; ...
   '*.slx;*.mdl','Models (*.slx, *.mdl)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'on');

sessions=numel(filename);

%% 
cor=[];
incor=[];
omiss=[];
pre=[];
dfcor=[];
dfincor=[];
dfomiss=[];
dfpre=[];

% errors=[];
% pause
for i=1:numel(filename)
    filename{i}

    [cc, ic, om, pr, dcc, dic, dom, dpr]=controlincluded_61hz_4(filename{i});
    %     ncc=nan(size(cc,1),459); % create tables
    %     nic=nan(size(ic,1),459);
    %     nom=nan(size(om,1),459);F
    %     npr=nan(size(pr,1),459);
    
          
    cor = vertcat(cor, cc);
    if size(ic,1)>1
        incor = vertcat(incor, ic);
    end
    omiss = vertcat(omiss, om);
    if size(pr,1)>1
        pre = vertcat(pre, pr);
    end
    %     ndcc=nan(size(dcc,1),459);
    %     ndic=nan(size(dic,1),459);
    %     ndom=nan(size(dom,1),459);
    %     ndpr=nan(size(dpr,1),459);
    
    dfcor = vertcat(dfcor, dcc);
    if size(dfincor,1)>1
        dfincor = vertcat(dfincor, dic);
    end
    dfomiss = vertcat(dfomiss, dom);
    if size(dfpre,1)>1
        dfpre = vertcat(dfpre, dpr);
    end
end
% 
% plot([1/1000:1/1000:10.001],mean(cor))
% hold on
% plot([1/1000:1/1000:10.001],nanmean(incor),'r')
% plot([1/1000:1/1000:10.001],mean(omiss),'k')

figure(2) % plot dF/F for adjusted fluorescence
plot([-5-1/30.5:1/30.5:10],mean(cor))
hold on
yL = get(gca,'YLim');
line([0 0],yL,'Color','m');
line([5 5],yL,'Color','g');
title('GCaMP6m fluorescence during 5-choice trials')%title
xlabel('time(s)') % x-axis label
ylabel('dF/F') % y-axis label
plot([-5-1/30.5:1/30.5:10],mean(omiss),'k')
if numel(incor)>0
    plot([-5-1/30.5:1/30.5:10],mean(incor),'r')
end
if numel(pre)>0
    plot([-5-1/30.5:1/30.5:10],mean(pre),'m')
end


figure(3) % plot dF/F for adjusted fluorescence
plot([-5-1/30.5:1/30.5:10],mean(dfcor))
hold on
yL = get(gca,'YLim');
line([0 0],yL,'Color','m');
line([5 5],yL,'Color','g');
title('GCaMP6m fluorescence during 5-choice trials')%title
xlabel('time(s)') % x-axis label
ylabel('dF/F') % y-axis label
plot([-5-1/30.5:1/30.5:10],mean(dfomiss),'k')
if numel(dfincor)>0
    plot([-5-1/30.5:1/30.5:10],mean(dfincor),'r')
end
if numel(dfpre)>0
    plot([-5-1/30.5:1/30.5:10],mean(dfpre),'m')
end



allFiles=filename;

% [h,p,ci,stats]=ttest2(cor,omiss)