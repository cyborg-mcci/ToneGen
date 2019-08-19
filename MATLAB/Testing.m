%Testing space for ToneGen

% % NumBits = 12;
% % MC = 41;
% % NumSamples = 2^16;
% % FullScale = 2;
% % f_in = 1e4;
% % fs = f_in*NumSamples/MC;
% % sample_cycle_ratio = MC/NumSamples;

%ADC
% [DO,T] = ADC(sample_cycle_ratio,FullScale,NumBits, MC);
% plot(T,DO);
% title('Digital Output vs Normalised time');
% xlabel('Time (s)');
% ylabel('Digital Output');

% [snrADC, enobADC, pot_signal_B_ADC, f, PSDADC] = gs_fresp(DO', length(DO), fs, f_in,1);


% 
% %DACNoise
%  FCornerDAC = 1e4;
% TNoiseDAC = 1e-11;
% [DAC_Out,DAC_T,DACNoise] = DAC(FCornerDAC,FullScale,MC,NumBits,NumSamples,TNoiseDAC);
% Time = DAC_T/f_in;
% DAC Noise PSD
%[snrDAC, enobDAC, pot_signal_B_DAC, f, PSDDAC] = gs_fresp(DACNoise', length(DACNoise), fs, f_in,1);


%Attenuator Noise PSDDivider = 2;
% FCorner_Att = 1e5;
% TNoise_Att = 1e-9;
% [Attenuator_Output,Att_Noise]  = Attenuator(DAC_Out,Divider,FCorner_Att,TNoise_Att);
% [snrAtt, enobAtt, pot_signal_B_Att, f, PSDAtt] = gs_fresp(Att_Noise', length(Attenuator_Output), fs, f_in,1);


%Filtering Test
% Filter_TNoise = 1e-12;
% FCornerFilter = 1e8;
% TFnum = 1;
% TFden = [1e-4 1]; %fc = 1.6kHz, Low-pass
% Filtered_Att = Filtering(TFnum,TFden,Attenuator_Output,Time,FCornerFilter,Filter_TNoise);
% plot(Time,Filtered_Att)



% fs = 10e3;
% [Pxx, f] = periodogram(RThermal*sqrt(fs*2/pi), [], [], fs);
% 

% % fs = 10e3;

% t = 0:1/fs:(NumSamples-1)/fs;

% % TNoiseDAC = 1e-11;
% % FCornerDAC = 1000;
% % 
% % kf_DAC = FCornerDAC*(TNoiseDAC^2);
% sample_cycle_ratio = MC/NumSamples;


%[FlickerDAC,~] = f_alpha(NumSamples,kf_DAC,1,1);



% f_test = logspace(log10(0.1), log10(fs/2), 1000);
% Sn_fln_test = kf_DAC./f_test;
% in_fln_test = sqrt(Sn_fln_test);


% Sn_fln = kf_DAC * randn(1,NumSamples);
% s = tf('s');
% shape = s^(-0.5);
% Sn_fln = abs(lsim(shape, Sn_fln, t));
% in_fln = sqrt(Sn_fln);


% % in_fln = f_alpha(NumSamples,kf_DAC, 0.5, 1);
% % 
% % 
% % RThermal = randn(1,NumSamples)*TNoiseDAC;
%DAC_Noise = RThermal + FlickerDAC;

%DAC_Output = Dig_Out + DAC_Noise;
% % DAC_Output = 0;
% % DAC_Noise = 0;


% % [Pxx_th, f] = periodogram(RThermal*sqrt(fs*2/pi), [], [], fs);
% % [Pxx_fl, f_fl] = periodogram(in_fln*sqrt(fs*2/pi), [], [], fs);



% figure(1)
% loglog(f, sqrt(Pxx_th), f_test, in_fln_test, f_fl, sqrt(Pxx_fl))
%semilogx(f_test, 10*log10(Sn_fln_test), f_fl, 10*log10(Pxx_fl))
% grid on

% N = 2^10;
% t = linspace(0,2*pi,N);
% fi = 1e1;
% y = sin(2*pi*fi*t);
% fs = fi*N/41;
% 
% [snr, enob, pot_signal_B, f, PSD] = gs_fresp(y', N, fs, fi, 1);


%% Testing the DAC Flicker Noise

% k = physconst('Boltzmann'); %Boltzmann's constant
% Temperature = 300; %Temperature in Kelvin
% 
% % DAC parameters (ADI LTC1668)
% FCornerDAC = 1e4; %estimate; not in datasheet
% TNoiseDAC_i = 50e-12; %per sqrt(Hz) - Value from Datasheet
% 
% fs = 50e6; % DAC Sample Frequency
% 
% NumSamples = 2^16;
% 
% f_in_desired= 50e3;
% NearestPrime = max(primes(f_in_desired*(NumSamples)/fs));
% f_in = fs*NearestPrime/(NumSamples);
% 
% 
% 
% kf_DAC = FCornerDAC*(TNoiseDAC^2);
% sample_cycle_ratio = MC/NumSamples;
% 
% [Dig_Out,DAC_NormalisedTime] = ADC(sample_cycle_ratio,FullScale,num_bits,MC);
% 
% [FlickerDAC,~] = f_alpha(NumSamples,kf_DAC,0.5,1);
% Thermal = randn(1,NumSamples)*TNoiseDAC;
% RThermal = randn(1,NumSamples)*sqrt(4*physconst('Boltzmann')*Temperature/Converter_R);
% DAC_Noise = Thermal + FlickerDAC' + RThermal;
% DAC_Output = DAC_Noise;
% 
% Time = DAC_NormalisedTime_1/f_in;


NoSamples = 2^22;
thn = 50e-12;
fc_fln = 1000;
k_fln = fc_fln * (thn^2);

ts = 1e-9;

[Flicker, ~] = f_alpha(NoSamples, k_fln, 0.5, 1);
Thermal = randn(1,NoSamples)*thn;

%Noise = Flicker' + Thermal;
Noise = Flicker';

t = 0:ts:(NoSamples-1)*ts;

[DAC_Output_1_i_Spectrum, DACOutput_1_i_f_TS, PSD_OS, f_OS, DAC_Output_1_Window] = wall_fresp(Noise, t, @rectwin, 0);
loglog(f_OS, db2mag(PSD_OS), f_OS, sqrt(k_fln./f_OS))
