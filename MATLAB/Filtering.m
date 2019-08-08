function [Filter_Output] = Filtering(TFnum,TFden,AttenuatorOutput,Time,FCornerFilter,Filter_TNoise)
%%Returns filtered output of Attenuator
% Filter_TNoise in A/sqrt(Hz)

kf_Filter = FCornerFilter*Filter_TNoise^2;

TransferFunction = tf(TFnum,TFden);
Filtered = lsim(TransferFunction,AttenuatorOutput,Time)';

Filter_Thermal = randn(1,length(AttenuatorOutput))*Filter_TNoise;
[FilterFlicker,~] = f_alpha(length(AttenuatorOutput),kf_Filter,0.5,1);
Filter_Noise = Filter_Thermal + FilterFlicker';

Filter_Output = Filtered + Filter_Noise;

end
