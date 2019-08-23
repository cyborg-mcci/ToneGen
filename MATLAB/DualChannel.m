%% DualChannel.m
%%  Microelectronic Circuits Centre Ireland (www.mcci.ie)
% 
%% 
% *Filename: *    Master.m
%%                    
% *Written by: *  Brendan O'Callaghan
%% 
% *Created on:*  15th August 2019
% 
% *Revised on:*   -
% 
% 
% 
% *File Description:*
% 
% _Script modelling a low noise current source with dual channel DACs and 
%   multiplication for uncorrelated noise cancelling_
% 
% 
% 
% _* Copyright 2019 Anthony Wall, Brendan O'Callaghan, MCCI, Tyndall, UCC*_

%% Initialisation Section

clearvars -EXCEPT nch pch
global ON OFF s FigureCounter k;

ON = 1;
OFF = 0;

s = tf('s');

FigureCounter = 0;

maxNumCompThreads(12);

%% Parameter declaration

%_i indicates a current quantity
%_v indicates a voltage quantity

k = physconst('Boltzmann'); %Boltzmann's constant
Temperature = 300; %Temperature in Kelvin

% DAC parameters (ADI LTC1668)
FCornerDAC = 1e4; %estimate; not in datasheet
TNoiseDAC_i = 50e-12; %per sqrt(Hz) - Value from Datasheet
FullScale_i = 20e-3; % Value from Datasheet - this is a pkpk value
num_bits = 16;
fs = 50e6; % DAC Sample Frequency

% Input signal parameters and coherent sampling condition
NumSamples = 2^22;
f_in_desired= 50e3;
NearestPrime = max(primes(f_in_desired*(NumSamples)/fs));
f_in = fs*NearestPrime/(NumSamples);

%Howland op-amp (AD797)
e_n = 0.9e-9; %%op-amp input referred voltage noise density
AOL = 10000; %%open-loop gain
ep = 0.1; %%mismatch error
GBWP = 110E6; %%gain-bandwidth product

%Current DAC Output to Voltage Conversion
Converter_R = 1e3; % Resistor at output of DAC to convert Current to voltage

% Sample RLC filter at the DAC output: parameters
f_f1 = f_in; % resonant frequency
BW_f1 = 5000; %bandwidth
L_f1 = 1e-3; %inductance 
Rpar_f1 = 50e-3; %Parasitic resistance of the inductor
C_f1 = 1/(4*(pi^2)*(f_f1^2)*L_f1); %Capacitance
num_f1 = [1/(L_f1*C_f1)]; %Transfer function numerator coefficients, in descending powers of s
den_f1 = [1 Rpar_f1/L_f1 1/(L_f1*C_f1)]; %Transfer function denominator coefficients, in descending powers of s
filter_1 = tf(num_f1,den_f1); %Filter transfer function
FCornerFilter_1 = 0; %Filter 1/f noise corner frequency
Filter_TNoise_1_v = sqrt(4*k*Temperature*Rpar_f1); %Filter thermal noise

% Analog multiplier (Mixer) parameters: (HA-2556)
FCorner_Mixer = 1e4; %1/f noise corner frequency - Estimated Value not in Datasheet?
TNoise_Mixer_v = 0.15e-6; %specified at f = 1kHz (worst case) - Output Referred Noise

%Sample LC current-mode filter at Howland output
f_f2 = 2*f_f1;
L_f2 = 1e-3;
C_f2 = 1/(4*(pi^2)*(f_f2^2)*L_f2); %% MAKE SURE THIS IS CORRECT
Rpar_f2 = 50e-3;
ZL = 1e4;
num_f2 = [L_f2 Rpar_f2];
den_f2 = [L_f2*C_f2*ZL (Rpar_f2*ZL*C_f2 + L_f2) (ZL + Rpar_f2)];
filter_2= tf(num_f2,den_f2);
FCornerFilter_2 = 0;
Filter_TNoise_2_i = sqrt(4*k*Temperature/Rpar_f2);

%% Model
Phase_Shift = 0;
% DAC Output, includes Converter_R noise 
[DAC_Output_1_i,DAC_NormalisedTime_1] = DAC(FCornerDAC,FullScale_i,NearestPrime,num_bits,NumSamples,TNoiseDAC_i,fs,Phase_Shift);
[DAC_Output_2_i,DAC_NormalisedTime_2] = DAC(FCornerDAC,FullScale_i,NearestPrime,num_bits,NumSamples,TNoiseDAC_i,fs,Phase_Shift);

% Scaling time vector depending on the input frequency
Time = DAC_NormalisedTime_1/f_in;

% PSD of the DAC output
[DAC_Output_1_i_Spectrum, DACOutput_1_i_f_TS, PSD_DAC_Output_1_i, DAC_Output_1_f_OS, DAC_Output_1_Window] = wall_fresp(DAC_Output_1_i, Time, @blackman, 0);
%[DAC_Output_2_i_Spectrum, DAC_Output_2_i_f_TS, PSD_DAC_Output_2_i, DAC_Output_2_f_OS, DAC_Output_2_Window] = wall_fresp(DAC_Output_1_i, Time, @hann, 0);

% Plotting DAC Output in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,DAC_Output_1_i)
xlabel('Time (s)')
ylabel('DAC output 1 (Current)')
title('DAC output 1 (A)')
subplot(2,1,2)
semilogx(DAC_Output_1_f_OS,PSD_DAC_Output_1_i)
grid on
xlabel('Frequency (Hz)')
ylabel('PSD')

String = sprintf('SNR of Current Coming from DAC: %f dB', snr(DAC_Output_1_i));
disp(String)
% FigureCounter = FigureCounter + 1;
% figure(FigureCounter)
% clf
% subplot(4,1,3)
% plot(Time,DAC_Output_2_i)
% xlabel('Time (s)')
% ylabel('DAC output 2 (Current)')
% title('DAC output 2 (A)')
% subplot(4,1,4)
% semilogx(DAC_Output_2_f_OS,PSD_DAC_Output_2_i)
% xlabel('Frequency (Hz)')
% ylabel('PSD')


% TIA to create two voltage channels
Channel_1_v = TIA_Converter(DAC_Output_1_i,Converter_R,Temperature,e_n);
Channel_2_v = TIA_Converter(DAC_Output_2_i,Converter_R,Temperature,e_n);

String = sprintf('SNR of Voltage Coming from TIA: %f dB', snr(Channel_1_v));
disp(String)

% Voltage signal
[Channel_1_v_Spectrum, Channel_1_v_f_TS, PSD_Channel_1_v, Channel_1_v_f_OS, Channel_1_v_Window] = wall_fresp(Channel_1_v, Time, @hann, 0);
% [Channel_2_v_Spectrum, Channel_2_v_f _TS, PSD_Channel_2_v, Channel_2_v_f _OS, Channel_2_v_Window] = wall_fresp(Channel_2_v, Time, @hann, 0);


% Plotting Channel 1 (V) in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Channel_1_v)
grid on
xlabel('Time (s)')
ylabel('TIA output (Voltage)')
title('TIA output (V)')
subplot(2,1,2)
semilogx(Channel_1_v_f_OS,PSD_Channel_1_v)
grid on
xlabel('Frequency (Hz)')
ylabel('PSD')

% Channel 2
% subplot(4,1,3)
% plot(Time,Channel_2_v)
% xlabel('Time (s)')
% ylabel('DAC output (Voltage)')
% title('DAC output (V)')
% subplot(4,1,4)
% semilogx(Channel_2_f _OS,PSD_Channel_2_v)
% xlabel('Frequency (Hz)')
% ylabel('PSD')


% Filter output
Filter_Output_1_v = Filtering(num_f1,den_f1,Channel_1_v,Time,FCornerFilter_1,Filter_TNoise_1_v,fs,0);
Filter_Output_2_v = Filtering(num_f1,den_f1,Channel_2_v,Time,FCornerFilter_1,Filter_TNoise_1_v,fs,0);

[Filter_Output_1_v_Spectrum, Filter_Output_1_v_f_TS, PSD_Filter_Output_1_v, Filter_Output_1_v_f_OS, Filter_Output_1_v_Window] = wall_fresp(Filter_Output_1_v, Time, @hann, 0);
% [Filter_Output_2_v_Spectrum, Filter_Output_2_v_f_TS, PSD_Filter_Output_2_v, Filter_Output_2_v_f_OS, Filter_Output_2_v_Window] = wall_fresp(Filter_Output_2_v, Time, @hann, 0);

String = sprintf('SNR of Voltage Coming from the Bandpass Filter: %f dB', snr(Filter_Output_1_v));
disp(String)

% SettlingTime = stepinfo(sys,'SettlingTimeThreshold',0.005);

% Plotting Filter output in time and frequency domain (One-sided
% FFT),removing transient from time domain output by only plotting the
% second half of the data

FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Filter_Output_1_v)
grid on
xlabel('Time')
ylabel('Filter Output (c1)')
title('Filter output 1')
subplot(2,1,2)
semilogx(Filter_Output_1_v_f_OS,PSD_Filter_Output_1_v)
grid on

% FigureCounter = FigureCounter + 1;
% figure(FigureCounter)
% clf
% subplot(2,1,1)
% plot(Time(length(Time)/2 +1:end),Filter_Output_2_v(length(Time)/2 +1:end))
% xlabel('Time')
% ylabel('Filter Output (c2)')
% title('Filter output 2')
% subplot(2,1,2)
% semilogx(Filter_Output_2_v_f_OS,PSD_Filter_Output_2_v)


% Bode plot of DAC filter's transfer function, axis in Hz
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
options = bodeoptions;
options.FreqUnits = 'Hz';
bode(filter_1,options)
title('Filter for the DAC output')

Mixer_Output_v = Mixer(Filter_Output_1_v,Filter_Output_2_v,TNoise_Mixer_v,FCorner_Mixer,fs);
[Mixer_Spectrum, Mixer_f_TS, PSD_Mixer, Mixer_f_OS, Mixer_Window] = wall_fresp(Mixer_Output_v, Time, @hann, 0);

% Plotting Mixer output in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Mixer_Output_v)
xlabel('Time')
ylabel('Mixer Output')
title('Mixer output')
subplot(2,1,2)
semilogx(Mixer_f_OS,PSD_Mixer)

%Offsetting the squared sine wave to be centred at zero (need to fix this
%to only subtract average of steady-state data)
% Mixer_Output_offset_v = detrend(Mixer_Output_v, 'constant');


% Plotting Mixer output(Howland input) in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time(length(Time)/2 +1:end),Mixer_Output_v(length(Time)/2 +1:end))
xlabel('Time')
ylabel('Final Output (c1)')
title('Final output')
subplot(2,1,2)
semilogx(Mixer_f_OS,PSD_Mixer)
xlabel('Frequency(Hz)')
ylabel('PSD')


% Doubled output frequency due to squaring sine wave
f_out = 2*f_in;

%Howland current source: V -> I: calculations based on 'Biocompatible,high
%precision, wideband, improved Howland current source with lead-lag
%compensation' IEEE paper

R2A = 150e6;
R2B = 150e6;
R2 = R2A + R2B;
R3 = 1e6;
R4 = 30e6;
G = R4/R3;
R1 =(R2)/G; %%balanced resistor bridge,calculating R1 as a function of other Rs

Rx = ((R2A + R2B)^2)/(R1*R2B^2) + R4^2/(R3*R2B^2) + (R2A + R2B + R4)/(R2B^2);
Co = (R3 + R4)*(R1 + R2A + R2B)/((2*pi*GBWP*R3)*(R1 + R2A)*R2B);

w = 2*pi*f_out;

RO1 = (R2B*(R1 + R2A))/(ep*(R2A + R2B)); %%accounting for resistor bridge mismatch
RO2 = (((R1 + R2A)*(R2B))/((R1 + R2A + R2B)))*((AOL*R3)/(R3+R4)+1); %%accounting for finite open-loop gain
XCo = 1/(w*Co); %%accounting for finite bandwidth

Zo = (1/RO1 + 1/RO2 + 1/XCo)^(-1); %%output impedance

ThermalNoise = sqrt(4*k*Temperature*Rx);
Opampnoise = e_n*(R1 + R2)/(R1*R2B);
TotalNoise = sqrt(ThermalNoise^2 + Opampnoise^2);

%Vin applied to non-inverting input
I_out_i = Mixer_Output_v*(R2A + R2B)/(R1*R2B);

%Vin applied to inverting input
%Output_current_i = Howland_Input_v*(-R4)/(R3*R2B);

Load_current_w_noise_i = I_out_i*(Zo/(ZL + Zo));

% Final_Output_i =









