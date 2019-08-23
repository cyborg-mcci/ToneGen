function Mixer_Output = Mixer(Channel_1,Channel_2,TNoiseMixer,FCornerMixer,fs)
%%%%Returns output of DAC (in Amps)
%%%%TNoiseDAC in A/sqrt(Hz)

Product = Channel_1.*Channel_2;

sigma_Thermal = TNoiseMixer*sqrt(fs/2);
Thermal = randn(1,length(Product))*sigma_Thermal;
var_Flicker = FCornerMixer*(TNoiseMixer^2)*log(length(Product)/2);

[FlickerMixer,~] = f_alpha(length(Product),var_Flicker,1,1);

Mixer_Noise = Thermal + FlickerMixer';
%Mixer_Output = Product + Mixer_Noise;
Mixer_Output = Product;
end