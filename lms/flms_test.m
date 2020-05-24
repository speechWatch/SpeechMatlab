% fnlms test
% load data
clc
clear all

aecdata = load('../testdata/aecdata.mat');
recData = aecdata.mic;
refData = aecdata.spk;

% set parameters
cfg.block_num  = 20;
cfg.block_len  = 128;                  % 8ms for 16khz sampling rate
cfg.filter_len = 2 * cfg.block_len;    % filter length  
% this filter length can supress the echo within (block_len * block_num / 16) ms
cfg.sig_min    = 1e-4;
cfg.threshold  = 0.1;                 % an important parameter
% are these parameters presented below can be formulated ?
cfg.gama       = 0.98;
cfg.alpha      = 0.01;
cfg.len        = length(recData);

% frequency NLMS
[res, filter_coeff] = FLMS(cfg, recData, refData);

% show the estimated room impulse response
figure
t = 1 : length(recData);
plot(t,recData, '-r', t, res,'-g')
legend('mic data','aec output data')
figure
plot(filter_coeff)
grid on