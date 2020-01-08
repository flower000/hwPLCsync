function [header] = hea_gener()
%hea_gener: To generator the header of the certain frame
%              for power lines
%   Reference: ITU-T G.9960 p89
%input: 
    global Fuc Fus beta Ndf Nhd Ngi Fsc N;
%% do ifft
    HEAD = (2*randi([0,1],[1,N])-1) + 1i*(2*randi([0,1],[1,N])-1);
    head = ifft(HEAD,N);  head(1) = 0;  % time-domain payload
    %pay = reshape(pay,1,N*(l-1));
%% adding the CP to OFDM block
    Ncp = Nhd + beta;
    CP_head = [head(N-Ncp+1:end),head];
%% Windowing
    % the first two payload
    len = N+Ncp;
    win = ones(len,1);
    win(1:beta) = 1.0 / (beta+1) * [1:beta];
    win(len-beta+1:end) = 1 - 1.0 / (beta+1) * [1:beta];
    winOFDM = win .* CP_head';
    % return value
    header = winOFDM;
%% display header
    %{
    figure;    hold on;
    plot(abs(winOFDM));
    set(gca,'xlim',[1,2850],'ylim',[0,0.05]);
    xlabel('discrete time');
    ylabel('amplitude');
    title('Header structure in time domain');
    legend('header');
%}
end

