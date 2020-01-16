% segment correlation
%{
function [coaCFO] = estCFOcoar(recFrame,auCorr,toa)
%estCFOcoar: To estimate the coarse CFO of certain OFDM frame within
%        preamble structure 
%   
    global N1 k1 N2 k2 N3 k3 N;
    Sum = 0;
    L = N / k1 * 2;
    for index = 1:N1-1
        %Sum = Sum - angle(auCorr(toa+(index-1)*L)) * N / L / pi;    % angle: the degree in rad mode 
        %temp = sum((toa+N1*len:toa+(N1+1)*len-1)...
        %    .*conj(recFrame(toa+(N1+1)*len:toa+(N1+2)*len-1)));
        
        Sum = Sum - angle(temp) * N / L / pi;
    end
    coaCFO = Sum/(N1-1);
end
%}



function [coaCFO] = estCFOcoar(recFrame,auCorr,toa)
%estCFOcoar: To estimate the coarse CFO of certain OFDM frame within
%        preamble structure 
%   
    global N1 k1 N2 k2 N3 k3 N;
    Sum = 0;
    L = N / k1 * 2;
    for index = 1:N1-1
        Sum = Sum - angle(auCorr(toa+(index-1)*L)) * N / L / pi;    % angle: the degree in rad mode 
    end
    coaCFO = Sum/(N1-1);
end




