clear all, close all, clc;
% Note: all the parameters, methods, structures and practical constants are from ITU-T G.9960
% author: Shiyuan Sun
% Finish date: None
% last modifying date: 
%% Initialize of the global parameters
% OFDM control parameters for power lines-----p102-103
global Fuc Fus beta Ndf Nhd Ngi Fsc N;
Fus = 25e6; Fuc = 0; Fsc = 24.4140625e3;
N = 2048; Ngi = N/32; Nhd = N/4; Ndf = N/4;
beta = N/8;
% preamble structure for power lines------p103
global N1 k1 N2 k2 N3 k3;
N1 = 7; k1 = 8;
N2 = 2; k2 = 8;
N3 = 0; k3 = 0;
% frame structure
global l;
l = 10;
% Channel effect
global delay I0_dB int_CFO frc_CFO Ts Rs fcw_k;
delay = floor(1.5 * N);
int_CFO = 0; frc_CFO = 0;
I0_dB = -2;
Ts = 1.0/Fsc/N;
Rs = Fsc * N;
fcw_k = N/100;
SNR = 10;
% detection
global coefficient slide;
coefficient = 1/9;
slide = N/k1/2;

%% generate header and payload with PLC standards in ITU-T G.9960
header = hea_gener();       % header of a certain frame
payload = pay_gener();      % payload of a certain frame
%% generate preamble for power lines with 50MHz
pream = pream_gener(2);      % preamble of a certain frame
frame = frame_gener(pream,header,payload);
%% go through the certain power line with 50MHz
recFrame = passChannel(frame,SNR);
%% Receiver: do detectind, timing and frequency evaluation
[est_time,est_fc] = detection(recFrame);    % the delay estimator and CFO estimator


