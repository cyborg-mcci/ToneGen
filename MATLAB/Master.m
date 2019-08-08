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

NumBits = 16;
MC = 41;
NumSamples = 8192;
FullScale = (1e-3)*100;
f_in = 1e4;
fs = f_in*NumSamples/MC;

%%%% DAC
FcornerDAC = 1e4;
DAC_ThermalN = 50e-12;
[DAC_Output,Normalised_Time,DAC_Noise] = DAC(FcornerDAC,FullScale,MC,NumBits,NumSamples,DAC_ThermalN);
Time = Normalised_Time/f_in;

%%%%Attenuator
Divider = 2;
FCorner_Att =1e3;
Att_ThermalN = 1e-12;
[AttOutput,Att_Noise] = Attenuator(DAC_Output,Divider,FCorner_Att,Att_ThermalN);

%%%%Filter
Filter_TNoise = 1e-12;
FCornerFilter = 1e8;
TFnum = 1;
TFden = [1e-4 1]; %fc = 1.6kHz
Filtered_Att = Filtering(TFnum,TFden,AttOutput,Time,FCornerFilter,Filter_TNoise);

% p1 = subplot(3,1,1);
% plot(Time,DAC_Output)
% title('DAC')
% p2 = subplot(3,1,2);
% plot(Time,AttOutput)
% title('Attenuator')
% p3 = subplot(3,1,3);
% plot(Time,Filtered_Att)
% title('Filtered')
% linkaxes([p1,p2,p3],'xy');

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(DAC_Noise', length(DAC_Noise), fs, f_in,1);




