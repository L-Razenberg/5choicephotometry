function fSig=highpassF(sign,fs,fc)

orderN=3;
Wn=fc/(fs/2);
[B,A] = butter(orderN,Wn,'high');

fSig=filtfilt(B,A,sign);