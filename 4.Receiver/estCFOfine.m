function [finCFO] = estCFOfine(recFrame,toa,coaCFO)
%estCFOfine: To estimate the accuracy CFO of certain OFDM frame within
%        preamble structure 
%   
    global N;
    global N1 k1 N2 k2 N3 k3;
    global delay I0_dB int_CFO frc_CFO Rs fcw_k;
    % modifying recFrame with estimated parameter 'coaCFO'
    modiFrame = recFrame.*exp(-2*pi*1i*coaCFO*(0:length(recFrame)-1)'/ N);
    % estimate the 'finCFO' by modified frame.
    len = N/k1;     L = N / k1 * 2;
    corr = sum(modiFrame(toa+N1*len:toa+(N1+1)*len-1).*conj(modiFrame(toa+(N1+1)*len:toa+(N1+2)*len-1)));
    finCFO = - angle(corr) * N / L / pi;
end

