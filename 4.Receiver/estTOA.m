function [auCorr,toa] = estTOA(recFrame)
%estTOA: To estimate the time of arrival of certain OFDM frame within
%        preamble structure 
%input:
%   recFrame: the consecutive discrete input into receiver.
%output:
%   auCorr: the SAC curve;
%   toa: the estimated time of arrival for the first 'short' symbol 
%% Acquisition
%SAC implementation
    global Fuc Fus beta Ndf Nhd Ngi Fsc N;
    global N1 k1 N2 k2 N3 k3;
    global coefficient;
    load '.\4.Receiver\looup_table.mat' 'looup_table';
    len = length(recFrame);
    auCorr = zeros(len-2.0*N/k1+1,1);      % auto-correlation of the received signal
    auCorr(1) = sum(recFrame(1:N/k1).*conj(recFrame(N/k1+1:2*N/k1)));
    for index = 2:len-2.0*N/k1+1
        auCorr(index) = auCorr(index-1) + recFrame(index-1+N/k1).*conj(recFrame(index-1+2*N/k1))...
            - recFrame(index-1).*conj(recFrame(index-1+N/k1));
    end
    plot(abs(auCorr));
%pick the location within first 'short' symbol using threshold
    MAXV = 10000;
    maxval = max(auCorr);   minval = min(auCorr);
    threshold = maxval - coefficient * (maxval - minval);
    f = abs(auCorr-threshold);   loc = zeros(1,4);
    for index = 1:4
        [~,loc(index)] = min(f);
        f(loc(index)) = MAXV;
    end
    [pickone,~] = min(loc);
%% Time tracking
    start = pickone + N/k1;     % within the second 'short' symbol
    candi = fft(recFrame(start:start+N/k1*4-1));
    temp = repmat(candi,[1,N/k1]);
    % normalization of the lookup-table using energy
    aveT = sum(abs(looup_table(:,1)).^2);    aveP = sum(abs(candi).^2);   
    temp = temp / sqrt(aveP / aveT);
    % calculating and comparing
    MSE = sum(abs((temp-looup_table)).^2)/(4*N1/k1);
    [~,shift] = min(MSE);       % shift index in second 'short' symol
    toa = pickone - (shift - 1);
%{
    load 'TEMP.mat' 'TEMP';
    load 'preamble.mat' 'preamble';
    load '.\4.Receiver\looup_table.mat' 'looup_table';
    load 'factor.mat' 'factor';
    load 'frame.mat' 'frame';
    %}
    end

