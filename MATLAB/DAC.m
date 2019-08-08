function [DAC_Output,DAC_NormalisedTime,DAC_Noise] = DAC(FCornerDAC,FullScale,MC,num_bits,NumSamples,TNoiseDAC)
%%%%Returns output of DAC
%%%%TNoiseDAC in A/sqrt(Hz)

kf_DAC = FCornerDAC*(TNoiseDAC^2);
sample_cycle_ratio = MC/NumSamples;

[Dig_Out,DAC_NormalisedTime] = ADC(sample_cycle_ratio,FullScale,num_bits,MC);

[FlickerDAC,~] = f_alpha(NumSamples,kf_DAC,0.5,1);
Thermal = randn(1,NumSamples)*TNoiseDAC;
RThermal = randn(1,NumSamples)*sqrt(4*physconst('Boltzmann')*293*100);
DAC_Noise = Thermal + FlickerDAC' + RThermal;

DAC_Output = Dig_Out + DAC_Noise;


end

