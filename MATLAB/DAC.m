function [DAC_Output,DAC_NormalisedTime] = DAC(FCornerDAC,FullScale,NearestPrime,num_bits,NumSamples,TNoiseDAC,fs,Phase_Shift)
%%%% A function which models the sine wave output of the DAC

%%%% INPUT PARAMETERS
%%% FCornerDAC - the corner frequency of the 1/f noise (unit = Hz)
%%% FullScale - peak-to-peak amplitude of the sine wave (A)
%%% NearestPrime - The nearest prime number of cycles governed by the
% coherent sampling condition
%%% num_bits - the bit precision of the DAC
%%% NumSamples - the total number of sample points 
%%% TNoiseDAC - the DAC Thermal noise (unit = A/sqrt(Hz))
%%% fs = DAC sampling frequency (unit = Hz) (LTC1668: 50 MSPS)
%%% Phase shift of the sine wave (units = rads)

%%%%OUTPUT PARAMETERS
%%% DAC_Output - DAC Sine wave (including 1/f and thermal noise) (unit = A)
%%% DAC_NormalisedTime = Normalised Time Vector

sample_cycle_ratio = NearestPrime/NumSamples;

[Dig_Out,DAC_NormalisedTime] = Digitiser(sample_cycle_ratio,FullScale,num_bits,NearestPrime,Phase_Shift);

sigma_Thermal = TNoiseDAC*sqrt(fs/2);
Thermal = randn(1,NumSamples)*sigma_Thermal;
var_Flicker = FCornerDAC*(TNoiseDAC^2)*log(NumSamples/2);

[FlickerDAC,~] = f_alpha(NumSamples,var_Flicker,1,1);


DAC_Noise = Thermal + FlickerDAC';
DAC_Output = Dig_Out + DAC_Noise;

end

