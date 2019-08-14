%Dual channel DAC w/ mixing

%DAC output
FCornerDAC = 1e6;
FullScale = 1;
k = physconst('Boltzmann');
num_bits = 16;
%f_in = 1e3/(2*pi);
f_in = 50e3;
fs = 50e6;
NumSamples = 2^18;
NearestPrime = max(primes(f_in*(NumSamples/2)/fs));
f_in = fs*NearestPrime/NumSamples;

TNoiseDAC = 50e-12;
[DAC_Output,DAC_NormalisedTime,~] = DAC(FCornerDAC,FullScale,NearestPrime,num_bits,NumSamples,TNoiseDAC);
Time = DAC_NormalisedTime/f_in;

[snr1, enob1, pot_signal_B1, f1, PSD1] = gs_fresp(DAC_Output', NumSamples, fs, f_in, 1);

figure(1)
clf
subplot(2,1,1)
plot(Time,DAC_Output)
xlabel('Time (s)')
ylabel('DAC output (Current)')
title('DAC output (A)')
subplot(2,1,2)
semilogx(f1,PSD1)
xlabel('Frequency (Hz)')
ylabel('PSD')

% -> Voltage 
R_sig = 1e3;
%Two voltage channels
c1 = DAC_Output*R_sig;
c2 = DAC_Output*R_sig;

[snr2, enob2, pot_signal_B2, f2, PSD2] = gs_fresp(c1', NumSamples, fs, f_in, 1);
[snr3, enob3, pot_signal_B3, f3, PSD3] = gs_fresp(c2', NumSamples, fs, f_in, 1);

figure(2)
clf
subplot(2,1,1)
plot(Time,c1)
xlabel('Time (s)')
ylabel('DAC output (Voltage)')
title('DAC output (V)')
subplot(2,1,2)
semilogx(f2,PSD2)
xlabel('Frequency (Hz)')
ylabel('PSD')
%channel 2 is the exact same
% subplot(4,1,3)
% plot(Time,c2)
% xlabel('Time (s)')
% ylabel('DAC output (Voltage)')
% title('DAC output (V)')
% subplot(4,1,4)
% semilogx(f3,PSD3)
% xlabel('Frequency (Hz)')
% ylabel('PSD')


%Filter parameters+transfer function (resonant frequency = 1000/2pi Hz,BW = 50Hz)
L = 1e-3;
R = 50*L;
C = 0.04e-6;
num1 = [R/L 0];
den1 = [1 R/L 1/(L*C)];
num2 = [R/L 0];
den2 = [1 R/L 1/(L*C)];
FCornerFilter = 0;
Filter_TNoise = sqrt(4*k*300*R);

Filter_Output_1 = Filtering(num1,den1,c1,Time,FCornerFilter,Filter_TNoise);
Filter_Output_2 = Filtering(num2,den2,c2,Time,FCornerFilter,Filter_TNoise);

[snr4, enob4, pot_signal_B4, f4, PSD4] = gs_fresp(Filter_Output_1(length(Time)/2 +1:end)', NumSamples/2, fs, f_in, 1);
[snr5, enob5, pot_signal_B5, f5, PSD5] = gs_fresp(Filter_Output_2(length(Time)/2 +1:end)', NumSamples/2, fs, f_in, 1);

figure(3) 
clf
subplot(2,1,1)
plot(Time(length(Time)/2:end),Filter_Output_1(length(Time)/2:end))
xlabel('Time')
ylabel('Filter Output (c1)')
title('Filter output')
subplot(2,1,2)
semilogx(f4,PSD4)



FCorner_Mixer = 1e3;
TNoise_Mixer = 1e-12;
kf_Mixer = FCorner_Mixer*TNoise_Mixer^2;

[Flicker_Mixer,~] = f_alpha(NumSamples,kf_Mixer,0.5,1);
Thermal_Mixer = randn(1,NumSamples)*TNoise_Mixer;
kf_Mixer = FCorner_Mixer*TNoise_Mixer.^2;

Mixer_Noise = Flicker_Mixer' + Thermal_Mixer;
Mixer_Output = Filter_Output_1.*Filter_Output_2 + Mixer_Noise;
num3 = [R/L 0];
den3 = [1 R/L 1/(L*C)];
Final_Output = Filtering(num3,den3,Mixer_Output,Time,FCornerFilter,Filter_TNoise);


f_out = 2*f_in;


