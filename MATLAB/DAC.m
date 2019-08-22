function [DAC_Output,DAC_NormalisedTime] = DAC(FCornerDAC,FullScale,NearestPrime,num_bits,NumSamples,TNoiseDAC,fs,Phase_Shift)
%%%%Returns output of DAC (in Amps)
%%%%TNoiseDAC in A/sqrt(Hz)
%global k;

sample_cycle_ratio = NearestPrime/NumSamples;

[Dig_Out,DAC_NormalisedTime] = ADC(sample_cycle_ratio,FullScale,num_bits,NearestPrime,Phase_Shift);

sigma_Thermal = TNoiseDAC*sqrt(fs/2);
Thermal = randn(1,NumSamples)*sigma_Thermal;
var_Flicker = FCornerDAC*(TNoiseDAC^2)*log(NumSamples/2);

[FlickerDAC,~] = f_alpha(NumSamples,var_Flicker,1,1);


DAC_Noise = Thermal + FlickerDAC';
% DAC_Noise = Thermal + FlickerDAC';
DAC_Output = Dig_Out + DAC_Noise;

end

