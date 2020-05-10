% this script is used to test the overlap¡¡save method
addpath(genpath('./testdata'));

fs = 16000;
h_ms = 200 / 1000;
decay_time = 80 / 1000;
decay_len  = round(decay_time * fs);
h_len = round(h_ms * fs);
h = randn(1,h_len);
h_amp = 10.^(-(0:h_len-1)/decay_len);
h = h .* h_amp;
h = h/sqrt(sum(h.^2));
t = 0 : 1/fs : (h_len - 1)/fs;
%%
[data fin] = audioread('CleanSpeech.wav');
spk = resample(data,fs, fin);

%% filter
mic1 = filter(h, 1, spk);
% audiowrite('output.wav',mic1,fs);

%% OLS
blockLen = h_len;
frmLen  = floor(length(data) / blockLen);
dataVec = zeros(2 * blockLen,1);
h_      = zeros(2 * blockLen,1);
h_(1:blockLen) = h;
mic2    = zeros(size(mic1));
validLen = frmLen * blockLen;
for i = 1 : frmLen
    index = (i - 1) * blockLen + 1 : i * blockLen;
    dataVec(1 : blockLen) = dataVec(blockLen + 1: end);
    dataVec(blockLen + 1: end) = spk(index);
    Ytmp = fft(dataVec) .* fft(h_);
    Ytmp = ifft(Ytmp);
    mic2(index) = Ytmp(blockLen + 1: end);
end
diff1 = mic2(1:validLen) - mic1(1:validLen);
fprintf('OLSÎó²îºÍÎª£º%f\n',sum(abs(diff1)));