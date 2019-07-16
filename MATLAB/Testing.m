%Testing space for ToneGen

NumBits = 12;
MC = 41;
NumSamples = 8192;
FullScale = 2;
f_in = 1e4;
fs = f_in*NumSamples/MC;
sample_cycle_ratio = MC/NumSamples;

%ADC
% [DO,T] = ADC(sample_cycle_ratio,FullScale,NumBits, MC);
% plot(T,DO);
% title('Digital Output vs Normalised time');
% xlabel('Time (s)');
% ylabel('Digital Output');

% [snrADC, enobADC, pot_signal_B_ADC, f, PSDADC] = gs_fresp(DO', length(DO), fs, f_in,1);



%DACNoise
% FCornerDAC = 1e4;
% TNoiseDAC = 1e-11;
% [DAC_Out,DAC_T,DACNoise] = DAC(FCornerDAC,FullScale,MC,NumBits,NumSamples,TNoiseDAC);
% Time = DAC_T/f_in;
% % DAC Noise PSD
% [snrDAC, enobDAC, pot_signal_B_DAC, f, PSDDAC] = gs_fresp(DACNoise', length(DACNoise), fs, f_in,1);


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






