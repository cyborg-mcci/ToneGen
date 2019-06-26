%% Master.m
%%  Microelectronic Circuits Centre Ireland (www.mcci.ie)
% 
%% 
% *Filename: *    Master.m
%%                    
% *Written by: *  Brendan O'Callaghan
%% 
% *Created on:*  20th June 2019
% 
% *Revised on:*   -
% 
% 
% 
% *File Description:*
% 
%  _Script creating ideal Sinewave Lookup values & Testing them in an FFT_
% 
% 
% 
% _* Copyright 2019 Anthony Wall, Brendan O'Callaghan, MCCI, Tyndall, UCC*_

%% Initialisation Section

clearvars -EXCEPT nch pch
global ON OFF s FigureCounter f;

ON = 1;
OFF = 0;

s = tf('s');

FigureCounter = 0;

maxNumCompThreads(12);

%% Parameter Declaration
<<<<<<< HEAD
sample_cycle_ratio = 739;
=======
MC=41;
N=2048;
sample_cycle_ratio = MC/N;
>>>>>>> parent of 7d8c2e9... Including basic filter test
Full_Scale = 2;
num_bits = 12;
f_in = 1e6;

[Dig_Out,Time_Out] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
Normalised_time = Time_Out.*(1/f_in); % normalising time series

[Stitched_DArray,Stitched_TArray] = StitchedArray(Dig_Out,Normalised_time);
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
plot(Stitched_TArray,Stitched_DArray);
grid on


OSR=1;
L = length(Stitched_TArray);
y = log2(L);
N = floor(y);

<<<<<<< HEAD
M = 59;
cycinsample = Stitched_TArray(end)*f_in;
numdatapts = length(Stitched_TArray);
fs = f_in*2^N/41; %%%%%%%% fi/fs = #cyclesperwindow/#dataptsperwindow

SDA = Stitched_DArray(1:2^N)';

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(SDA, 2^N, fs, fi, OSR);
=======
fs=f_in/sample_cycle_ratio;
OSR=1;

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(Dig_Out', N, fs, f_in, OSR);
>>>>>>> parent of 7d8c2e9... Including basic filter test

% num = [5];
% den = [1];
% t = Stitched_TArray(1:length(Stitched_DArray));

% Response = Filtering(num,den,Stitched_DArray,t);

