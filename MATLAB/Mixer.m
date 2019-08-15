function Mixer_Output = Mixer(Channel_1,Channel_2,TNoise_Mixer,FCorner_Mixer)
%%%%Returns output of DAC (in Amps)
%%%%TNoiseDAC in A/sqrt(Hz)

Product = Channel_1.*Channel_2;

kf_Mixer = FCorner_Mixer*(TNoise_Mixer^2);
[Flicker_Mixer,~] = f_alpha(length(Channel_1),kf_Mixer,0.5,1);
Thermal = randn(1,length(Channel_1))*TNoise_Mixer;
Mixer_Noise = Thermal + Flicker_Mixer';
Mixer_Output = Product + Mixer_Noise;
end