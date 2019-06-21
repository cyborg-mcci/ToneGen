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
global ON OFF s FigureCounter T kb fmin fmax FreqPoints wmax wmin w_test f;

ON = 1;
OFF = 0;

s = tf('s');

FigureCounter = 0;

maxNumCompThreads(12);

sample_cycle_ratio = 1021;
Full_Scale = 2;
num_bits = 16;
f_in = 10e3;

[Dig_Out,Time_Out] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
Normalised_time = Time_Out.*(1/f_in); % normalising time series

[Stitched_DArray,Stitched_TArray] = StitchedArray(Dig_Out,Normalised_time);
plot(Stitched_TArray,Stitched_DArray);

fs = sample_cycle_ratio * f_in;
OSR=1;
L = length(Stitched_TArray);
y = log2(L);
N = floor(y);

Stitched_TArray = Stitched_TArray(1:2^N)';

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(Stitched_TArray, 2^N, fs, fi, OSR);

% [Stitched_Array, time_series] = StitchedArray(Dig_Array,sample_cycle_ratio,Full_Scale,num_bits);