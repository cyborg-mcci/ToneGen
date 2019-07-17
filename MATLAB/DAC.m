function [DAC_Output,DAC_NormalisedTime,DAC_Noise] = DAC(FCornerDAC,FullScale,MC,num_bits,NumSamples,TNoiseDAC)
%%%%Returns output of DAC
%%%%TNoiseDAC in A/sqrt(Hz)

kf_DAC = FCornerDAC*(TNoiseDAC^2);
sample_cycle_ratio = MC/NumSamples;

[Dig_Out,DAC_NormalisedTime] = ADC(sample_cycle_ratio,FullScale,num_bits,MC);

[FlickerDAC,~] = f_alpha(NumSamples,sqrt(kf_DAC),1,1);
RThermal = randn(1,NumSamples)*TNoiseDAC;
%DAC_Noise = RThermal + FlickerDAC;

%DAC_Output = Dig_Out + DAC_Noise;
DAC_Output = 0;
DAC_Noise = 0;

fs = NumSamples/max(DAC_NormalisedTime);
[Pxx_th, f] = periodogram(RThermal*sqrt(fs*2/pi), [], [], fs);
[Pxx_fl, ~] = periodogram(FlickerDAC*sqrt(fs*2/pi), [], [], fs);



figure(1)
loglog(f, sqrt(Pxx_th), f, sqrt(Pxx_fl))
grid on


end

