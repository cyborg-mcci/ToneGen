function Mixer_Output = Mixer(Channel_1,Channel_2,TNoiseMixer,FCornerMixer,fs)
%%%% A function which multiplies two signals together and returns the
%%%% product 

%%%% INPUT PARAMETERS
%%% Channel_1 - the first input channel to the mixer
%%% Channel_2 - the second input channel to the mixer
%%% TNoiseMixer - the mixer thermal noise (output-referred by default)
%%% FCornerMixer - the mixer 1/f noise corner frequency (unit = Hz)
%%% fs - the DAC sampling frequency

%%%%OUTPUT PARAMETERS
%%% Mixer_Output - the multplied signal

Product = Channel_1.*Channel_2;

sigma_Thermal = TNoiseMixer*sqrt(fs/2);
Thermal = randn(1,length(Product))*sigma_Thermal;
var_Flicker = FCornerMixer*(TNoiseMixer^2)*log(length(Product)/2);

[FlickerMixer,~] = f_alpha(length(Product),var_Flicker,1,1);

Mixer_Noise = Thermal + FlickerMixer';
Mixer_Output = Product + Mixer_Noise;

end