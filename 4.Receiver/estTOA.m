function [auCorr,toa] = estTOA(recFrame,label)
    if label == 1
        [auCorr,toa] = estTOA_A(recFrame);
    elseif label == 2
        [auCorr,toa] = estTOA_B(recFrame);
    elseif label == 3
        [auCorr,toa] = estTOA_C(recFrame);
    elseif label == 4
        [auCorr,toa] = estTOA_D(recFrame);
    end
end

function [auCorr,toa] = estTOA_A(recFrame)
%estTOA: To estimate the time of arrival of certain OFDM frame within
%        preamble structure 
%input:
%   recFrame: the consecutive discrete input into receiver.
%output:
%   auCorr: the SAC curve;
%   toa: the estimated time of arrival for the first 'short' symbol 
%% Acquisition method: A
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

function [auCorr,toa_seg] = estTOA_B(recFrame)
%% Acquisition method: B
% window length: five 'short' symbols
    global Fuc Fus beta Ndf Nhd Ngi Fsc N;
    global N1 k1 N2 k2 N3 k3;
    global coefficient slide;
    len = length(recFrame);
    segCorr = ones(len-(N1-2)*N/k1+1,1);      % auto-correlation of the received signal
    auCorr = ones(len-(N1-2)*N/k1+1,1);
    %auCorr(1) = sum(recFrame(1:N/k1).*conj(recFrame(N/k1+1:2*N/k1)));
    for index = 1:len-(N1-2)*N/k1+1
        %auCorr(index) = auCorr(index-1) + recFrame(index-1+N/k1).*conj(recFrame(index-1+2*N/k1))...
        %    - recFrame(index-1).*conj(recFrame(index-1+N/k1));
        temp = recFrame(index:index+(N1-2)*N/k1-1);
        for k = 1:N1-3
            temp1 = temp(1+(k-1)*N/k1:k*N/k1);    temp2 = temp(1+k*N/k1:(k+1)*N/k1);
            auCorr(index) = auCorr(index) * abs(sum(temp1.*conj(temp2)));
            segCorr(index) = segCorr(index) * Segment(temp1,temp2,slide);
        end
    end
    % display
    figure;     hold on;
    subplot(1,2,1);    plot(abs(segCorr));      title('segment');
    subplot(1,2,2);    plot(abs(auCorr));       title('normal Correlation');
    % pick the beginning
    [~,toa_seg] = max(segCorr);     [~,toa_nor] = max(abs(auCorr));
    toa_seg = toa_seg - beta/2;     toa_nor = toa_nor - beta/2;
end

function [auCorr,toa_seg] = estTOA_C(recFrame)
%% Acquisition method: C
% window length: seven 'short' symbols
    global Fuc Fus beta Ndf Nhd Ngi Fsc N;
    global N1 k1 N2 k2 N3 k3;
    global coefficient slide;
    len = length(recFrame);
    segCorr = ones(len-(N1)*N/k1+1,1);      % auto-correlation of the received signal
    auCorr = ones(len-(N1)*N/k1+1,1);
    %auCorr(1) = sum(recFrame(1:N/k1).*conj(recFrame(N/k1+1:2*N/k1)));
    for index = 1:len-(N1)*N/k1+1
        %auCorr(index) = auCorr(index-1) + recFrame(index-1+N/k1).*conj(recFrame(index-1+2*N/k1))...
        %    - recFrame(index-1).*conj(recFrame(index-1+N/k1));
        temp = recFrame(index:index+(N1)*N/k1-1);
        for k = 1:N1-1
            temp1 = temp(1+(k-1)*N/k1:k*N/k1);    temp2 = temp(1+k*N/k1:(k+1)*N/k1);
            auCorr(index) = auCorr(index) * abs(sum(temp1.*conj(temp2)));
            segCorr(index) = segCorr(index) * Segment(temp1,temp2,slide);
        end
    end
    % display
    figure;     hold on;
    subplot(1,2,1);    plot(abs(segCorr));      title('segment');
    subplot(1,2,2);    plot(abs(auCorr));       title('normal Correlation');
    % pick the beginning
    [~,toa_seg] = max(segCorr);     [~,toa_nor] = max(abs(auCorr));
    toa_seg = toa_seg;     toa_nor = toa_nor;
end

function [auCorr,toa_seg] = estTOA_D(recFrame)
%% Acquisition method: D
% window length: seven 'short' symbols and two 'long' symbols
    global Fuc Fus beta Ndf Nhd Ngi Fsc N;
    global N1 k1 N2 k2 N3 k3;
    global coefficient slide;
    len = length(recFrame);
    segCorr = ones(len-(N1+N2)*N/k1+1,1);      % auto-correlation of the received signal
    auCorr = ones(len-(N1+N2)*N/k1+1,1);
    %auCorr(1) = sum(recFrame(1:N/k1).*conj(recFrame(N/k1+1:2*N/k1)));
    for index = 1:len-(N1+N2)*N/k1+1
        %auCorr(index) = auCorr(index-1) + recFrame(index-1+N/k1).*conj(recFrame(index-1+2*N/k1))...
        %    - recFrame(index-1).*conj(recFrame(index-1+N/k1));
        % Section one of the preamble
        temp = recFrame(index:index+(N1+N2)*N/k1-1);
        for k = 1:N1-1
            temp1 = temp(1+(k-1)*N/k1:k*N/k1);    temp2 = temp(1+k*N/k1:(k+1)*N/k1);
            auCorr(index) = auCorr(index) * abs(sum(temp1.*conj(temp2)));
            segCorr(index) = segCorr(index) * Segment(temp1,temp2,slide);
        end
        % Section two of the preamble
        for k = 1:N2-1
            temp1 = temp(1+N1*N/k1+(k-1)*N/k2:k*N/k2+N1*N/k1);    temp2 = temp(1+N1*N/k1+k*N/k2:(k+1)*N/k2+N1*N/k1);
            auCorr(index) = auCorr(index) * abs(sum(temp1.*conj(temp2)));
            segCorr(index) = segCorr(index) * Segment(temp1,temp2,slide);
        end
    end
    % display
    figure;     hold on;
    subplot(1,2,1);    plot(abs(segCorr));      title('segment');
    subplot(1,2,2);    plot(abs(auCorr));       title('normal Correlation');
    % pick the beginning
    [~,toa_seg] = max(segCorr);     [~,toa_nor] = max(abs(auCorr));
    toa_seg = toa_seg;     toa_nor = toa_nor;
end

