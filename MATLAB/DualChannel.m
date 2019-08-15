%Dual channel DAC w/ mixing

global FigureCounter;
FigureCounter = 0;
k = physconst('Boltzmann');

%DAC parameters
FCornerDAC = 1e6;
TNoiseDAC = 50e-12;
FullScale = 1;
num_bits = 16;
NumSamples = 2^22;

%Input signal parameters
f_in_desired= 50e3;
fs = 50e6;
NearestPrime = max(primes(f_in_desired*(NumSamples/2)/fs));
f_in = fs*NearestPrime/(NumSamples/2);

Converter_R = 1e3; %Resistor at output of DAC to convert to voltage signal

%DAC Output (in Amps), includes Converter_R noise 
[DAC_Output,DAC_NormalisedTime,~] = DAC(FCornerDAC,FullScale,NearestPrime,num_bits,NumSamples,TNoiseDAC,Converter_R);

%Scaling time vector depending on the input frequency
Time = DAC_NormalisedTime/f_in;


%FFT of the DAC output (A)
[DAC_Spectrum, DAC_f_TS, PSD_DAC, DAC_f_OS, DAC_Window] = wall_fresp(DAC_Output, Time, @rectwin, 0);

%[snr1, enob1, pot_signal_B1, f1, PSD1] = gs_fresp(DAC_Output', NumSamples, fs, f_in, 1);
[Spectrum, f_TS, PSD_OSdB, f_OS, Window] = wall_fresp(DAC_Output, Time, @hann, 1);

%Plotting DAC Output (A) in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,DAC_Output)
xlabel('Time (s)')
ylabel('DAC output (Current)')
title('DAC output (A)')
subplot(2,1,2)
semilogx(DAC_f_OS,PSD_DAC)
xlabel('Frequency (Hz)')
ylabel('PSD')

%Two voltage channels
c1 = DAC_Output*Converter_R;
c2 = DAC_Output*Converter_R;

%Voltage signal
[C1_Spectrum, C1_f_TS, PSD_C1, C1_f_OS, C1_Window] = wall_fresp(c1, Time, @rectwin, 0);
% [C2_Spectrum, C2_f_TS, PSD_C2, C2_f_OS, C2_Window] = wall_fresp(c2, Time, @rectwin, 0);

%Plotting DAC Output (V) in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,c1)
xlabel('Time (s)')
ylabel('DAC output (Voltage)')
title('DAC output (V)')
subplot(2,1,2)
semilogx(C1_f_OS,PSD_C1)
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

%Sample RLC Bandpass Filter parameters (resonant frequency = f0, BW = bandwidth)
%Fixed L,calculate C based on f0 and L
f0 = f_in;
BW = 50;
L = 1e-3;
R = 2*pi*BW*L;
C = 1/(4*(pi^2)*(f0^2)*L);
num = [R/L 0];
den = [1 R/L 1/(L*C)];
filter_1 = tf(num,den);
FCornerFilter = 0;
T = 300;
Filter_TNoise = sqrt(4*k*T*R);

%Filter output
Filter_Output_1 = Filtering(num,den,c1,Time,FCornerFilter,Filter_TNoise);
Filter_Output_2 = Filtering(num,den,c2,Time,FCornerFilter,Filter_TNoise);

[Filter_1_Spectrum, Filter_1_f_TS, PSD_Filter_1, Filter_1_f_OS, Filter_1_Window] = wall_fresp(Filter_Output_1, Time, @rectwin, 0);
% [snr5, enob5, pot_signal_B5, f5, PSD5] = gs_fresp(Filter_Output_2(length(Time)/2 +1:end)', NumSamples/2, fs, f_in, 1);

%Plotting Filter output in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Filter_Output_1)
xlabel('Time')
ylabel('Filter Output (c1)')
title('Filter output')
subplot(2,1,2)
semilogx(Filter_1_f_OS,PSD_Filter_1)

%Bode plot of DAC filter's transfer function
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
options = bodeoptions;
options.FreqUnits = 'Hz';
bode(filter_1,options)
title('DAC output Filter')

%Mixer parameters
FCorner_Mixer = 1e3;
TNoise_Mixer = 1e-12;

Mixer_Output = Mixer(Filter_Output_1,Filter_Output_2,TNoise_Mixer,FCorner_Mixer);
[M_Spectrum, M_f_TS, PSD_M, M_f_OS, M_Window] = wall_fresp(Mixer_Output, Time, @rectwin, 0);
%Mixer filter,resonates at f1 = 2*f0 (as squared sine wave has double the input
%frequency)
f1 = 2*f0;
BW1 = 50;
L1 = 1e-3;
R1 = 2*pi*BW*L1;
C1 = 1/(4*(pi^2)*(f1^2)*L1);
num1 = [R1/L1 0];
den1 = [1 R1/L1 1/(L1*C1)];
filter_2 = tf(num1,den1);

%Mixer output after filtering
Final_Output = Filtering(num1,den1,Mixer_Output,Time,FCornerFilter,Filter_TNoise);

[Final_Spectrum, Final_f_TS, PSD_Final, Final_f_OS, Final_Window] = wall_fresp(Final_Output, Time, @rectwin, 0);

Mixer_Output = detrend(Mixer_Output, 'constant');

%Plotting Mixer output in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Mixer_Output)
xlabel('Time')
ylabel('Mixer Output')
title('Mixer output')
subplot(2,1,2)
semilogx(M_f_OS,PSD_M)

%Plotting Final output in time and frequency domain (One-sided FFT)
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
subplot(2,1,1)
plot(Time,Final_Output)
xlabel('Time')
ylabel('Final Output (c1)')
title('Final output')
subplot(2,1,2)
semilogx(Final_f_OS,PSD_Final)

%Bode plot of mixer filter's transfer function
FigureCounter = FigureCounter + 1;
figure(FigureCounter)
clf
options = bodeoptions;
options.FreqUnits = 'Hz';
bode(filter_1,options)
title('Mixer output filter')
 

%Doubled output frequency due to squaring sine wave
f_out = 2*f_in;

%Chopper op-amp

