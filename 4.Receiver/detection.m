function [est_time,est_fc] = detection(recFrame)
%detection: Receiver detect the signal and implement synchronization
%   
    global N1 k1 N2 k2 N3 k3;
%% TOA(time of arrival) estimationL: SAC implementation
    [auCorr,toa] = estTOA(recFrame);
%% coarse estimation of the frequency
    [coaCFO] = estCFOcoar(auCorr,toa);
%% fine estimation of the frequency
    [finCFO] = estCFOfine(recFrame,toa,coaCFO);      % use coaCFO to modify the data in 'long' symbols
%% return value
    est_time =toa;
    est_fc = coaCFO + finCFO;
end

