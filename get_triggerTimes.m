function [trig]=get_triggerTimes(ADCchannel)
% function [TTLts]=get_triggerTimes(ADCchannel) will return timestamps on
% which the TTL was up (from 200ms piezo protocol)

ADChigh=find(ADCchannel<2);ADCd=diff(ADChigh);
firstHigh=find(ADCd>300);
trig=[ADChigh(1); ADChigh(firstHigh+1)]; %timestamp in samples

