function [fEst AEst]=DFTEstimate(clip, dX0)


%Using DFT to estimate the frequency 
T=dX0;
Fs=1/T;
L=size(clip,2);
NFFT=2^nextpow2(L);
fEst=0;
AEst=0;
for i=1:size(clip,1)
y=clip(i,:);    
Y=fft(y,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2+1);
coef=2*abs(Y(1:NFFT/2+1));
% plot(f,coef); 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

coef1=coef(2:end);
[row,col]=find(coef1==max(coef1));
fEst=fEst+f(row,col);
AEst=AEst+coef(row,col);
end
fEst=fEst/size(clip,1);
AEst=AEst/size(clip,1);