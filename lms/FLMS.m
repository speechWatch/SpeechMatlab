function [res, filter_coeff] = FLMS(cfg, recData, refData)
% This Scripts is used to implement the AEC based on the frequency normalized LMS
% input 
       % cfg : parameterts
       % recData : micrphone received data
       % refData : reference data for aec
% output
       % res : aec output
       % filter_coef : estimated room impulse response

% 1) FLMS filter initialization
dataLen       = cfg.len;
alpha         = cfg.alpha;
threshold     = cfg.threshold;
block_len     = cfg.block_len;
block_num     = cfg.block_num;
filter_len    = cfg.filter_len;
sig_min       = cfg.sig_min .* eye(filter_len);
Vec           = zeros (filter_len, block_num);
curVec        = zeros (filter_len, 1);
W             = zeros (filter_len, block_num);
e             = zeros (filter_len, 1);
phi           = zeros (filter_len, block_num);
PHI           = zeros (filter_len, block_num);
res           = zeros(dataLen,1);
block_number  = floor(dataLen / block_len);
Pi            = 0 .* eye(filter_len);

% 2£© start FLMS simulation
for i = 1 : block_number    
    gama                    = cfg.gama;
    avePow                  = zeros (filter_len, 1);
    Y                       = zeros (filter_len, 1);
    index                   = (i - 1) * block_len + 1;
    d                       = recData(index : index + block_len -1);
    curVec(1 : block_len)   = curVec(block_len+1 : end);
    curVec(block_len+1 : end)  = refData(index : index + block_len - 1); 
    Fvec                    = fft(curVec);
    % construct the block matrix
    Vec(:,2:end)            = Vec(:,1:end-1);
    Vec(:,1)                = Fvec;
    % calculate the average power
    for j = 1 : block_num
        avePow = avePow + conj(Vec(:,j)) .* Vec(:,j);
    end
    U                       = diag(avePow);
    Pi                      = gama * Pi + (1 - gama) .* U;
    avePi                   = Pi ./ block_num;
    D                       = inv(avePi + sig_min);
    
    for j = 1 : block_num
        Y = Y + Vec(:, j) .* W(: , j);
    end
    y                       = ifft(Y);
    e(block_len+1 : end)    = d - y(block_len+1 : end);
    E                       = fft(e); 
    % error normalization
    E1                      = D * E;
    % threshold for error
    absEf                   = max(abs(E1), threshold);
    absEf                   = ones(filter_len,1)*threshold./absEf;
    E1                      = E1.*absEf;
    % step size factor
    E1                      = alpha .* E1;
    %
    for j = 1 : block_num
        phi(:, j) = ifft(conj(Vec(: , j)) .* E1);
    end
    PHI(1 : block_len,:)    = phi(1 : block_len,:);
    % weight update
    W                       = W + fft(PHI); 

    % refiltering
    Y                       = zeros (filter_len, 1);
    for j = 1 : block_num
        Y = Y + Vec(:, j) .* W(: , j);
    end
    y                       = ifft(Y);
    e(block_len+1 : end)    = d - y(block_len+1 : end);
    
    % output
    res(index : index + block_len -1)      = real(e(block_len+1 : end));
    
end

% 3) calculate the weight values in time domain
filter_coeff = ifft(W);
filter_coeff = real(filter_coeff(1:block_len,:));
filter_coeff = reshape(filter_coeff, [block_num * block_len, 1]);
end


