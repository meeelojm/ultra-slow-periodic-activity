function tai_ang_fra = cal_tai_ang_fra(tai_ang_uni, tau)
    fra_rat.man = 120;% fps
    ifi.man = 1/fra_rat.man;
    
    tbin = ifi.man;

    kernellength = round(2*tau/tbin);
    kernel1 = exp(-(0:kernellength-1)/(tau/tbin));  % exponential kernel
    
    % Full convolution
    y_full = conv(tai_ang_uni, kernel1, 'full');
    
    % Trim so the output starts at the same time as the input (causal)
    tai_ang_fra = y_full(1:length(tai_ang_uni));
end