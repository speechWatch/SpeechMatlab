% This script is used to test the STFT and Synthsis for speech
clear 

addpath(genpath('./testdata'))

%% initialize the paramters
msg=audioread('CleanSpeech.wav');
filter=load('window.mat');
coef=filter.kwin;
fft_num=512;
offset_per=1/4;
overlap_per=1-offset_per;
frame_shift=offset_per*fft_num;
ch_num=size(msg,2);
msg_len=size(msg,1);
block_num=ceil(msg_len/frame_shift);
extra_len = block_num * frame_shift - msg_len;
xout=zeros(fft_num/2+1,block_num,ch_num,'like',msg);

%% STFT analysis
for ch_idx=1:ch_num
    x_tmp=msg(:,ch_idx).';
    x_tmp=[zeros(1,fft_num*overlap_per),x_tmp];
    for block_idx=1:block_num-1
        pos=(block_idx-1)*frame_shift;
        block_msg=x_tmp(pos+1:pos+fft_num);
        block_filter=block_msg.*coef.';
        block_fft=fft(block_filter);
        xout(:,block_idx,ch_idx)=block_fft(1,1:fft_num/2+1);
    end 
end
xout = [xout, zeros(fft_num/2+1,overlap_per/offset_per)];
%% Synthesis filter bank overlap add
msg_len=frame_shift*size(xout,2);
msg_out = zeros(1,msg_len,'like',xout);
tmp_out = zeros(1,msg_len+fft_num*overlap_per,'like',xout);
ch_idx=1;
j=sqrt(-1);
X=xout(:,:,ch_idx);
tdl = zeros(1,fft_num,'like',xout);
coef_fft=fft(coef);
for msg_idx = 1 : size(xout,2)
    X_fft=[X(:,msg_idx);conj(X(end-1:-1:2,msg_idx))];
    x_n=real(ifft(X_fft));
    y_n=x_n.*coef;%smooth the point between two different frames
    k=frame_shift*(msg_idx-1)+1;
    tmp_out(k:k+fft_num-1)=tmp_out(k:k+fft_num-1)+y_n.';
    msg_out(k:k+frame_shift-1)=tmp_out(k:k+frame_shift-1);%.*frame_shift/coef_fft(1);
end

% remove the extra zeros
msg_out = msg_out(1 : end - extra_len);

% remove the delay
delay = fft_num * overlap_per;
msg_out = msg_out(delay + 1 :end);

diff = mean(abs(real(msg_out) - msg'));
fprintf('The average synthesis error is: %d\n',diff);
