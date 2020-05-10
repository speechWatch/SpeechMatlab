clear all
clc

ploton = 1;
% This Script is used to generate the kasier window
pi=3.14159265;
kwin = kaiser(512, 1.9*pi);  % Kaiser Window
kwsigma = sqrt(sum(kwin.^2)/512*4);
kwin = kwin/kwsigma;
% analyze the performance of window
if ploton
    Wf      = fft(kwin);
    Wd      = Wf(1);
    W       = kwin;
    x       = zeros(1,1024);
    x(1:512)= kwin.';
    x((1+128):(512+128))=x((1+128):(512+128))+kwin.';
    x((1+128*2):(512+128*2))=x((1+128*2):(512+128*2))+kwin.';
    x((1+128*3):(512+128*3))=x((1+128*3):(512+128*3))+kwin.';
    x((1+128*4):(512+128*4))=x((1+128*4):(512+128*4))+kwin.';
    x=x.*128/Wd;
    plot(x)
end
save window kwin