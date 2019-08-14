function [Filtered] = Filtering(TFnum,TFden,Signal,Time,FCornerFilter,Filter_TNoise)
%%Returns filtered output of Attenuator
% Filter_TNoise in A/sqrt(Hz)

kf_Filter = FCornerFilter*(Filter_TNoise^2);

TransferFunction = tf(TFnum,TFden);

Filter_Thermal = randn(1,length(Signal))*Filter_TNoise;
[FilterFlicker,~] = f_alpha(length(Signal),kf_Filter,0.5,1);
Filter_Noise = Filter_Thermal + FilterFlicker';
Input = Signal + Filter_Noise;
Input = Signal;

Filtered = lsim(TransferFunction,Input,Time)';
Filtered = Filtered + Filter_Noise;

end
