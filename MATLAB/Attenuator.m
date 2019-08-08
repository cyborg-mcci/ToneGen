function [Attenuator_Output,Att_Noise]  = Attenuator(DAC_Output,Divider,FCorner_Att,TNoise_Att)
%%Returns output of Attenuator
%%%%TNoise_Att in A/sqrt(Hz)

Div_Output = DAC_Output/Divider;

kf_Att = FCorner_Att*TNoise_Att^2;
[Att_Flicker,~] = f_alpha(length(DAC_Output),kf_Att,0.5,1);
Att_Thermal = randn(1,length(DAC_Output))*TNoise_Att;

Att_Noise = Att_Flicker' + Att_Thermal;
Attenuator_Output = Div_Output + Att_Noise;

end