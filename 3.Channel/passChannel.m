function [recframe] = passChannel(frame,SNR)
%passChannel: add some NBI noise/ AWGN/ delay/ CFO or something else.
%   
    global N;
    global delay I0_dB int_CFO frc_CFO Rs fcw_k;
%% delay
    deFrame = [zeros(delay,1);frame;zeros(N,1)];
    %deFrame = [frame(end-delay+1:end);frame;zeros(N,1)];    % [last frame; this frame; next frame]

%% CFO
    CFO = int_CFO+frc_CFO;
    cfFrame = deFrame.*exp(2*pi*1i*CFO*(0:length(deFrame)-1)'/ N);      % additive in time domain
%% NBI
    nbiFrame = cfFrame + 10^(I0_dB/20)*exp(2*pi*1i*fcw_k*(0:length(cfFrame)-1)'/N);
%% AWGN
    recframe = awgn(nbiFrame,SNR,'measured');

    %recframe = nbiFrame;
end

