function [Filtered] = Filtering(TFnum,TFden,Signal,Time,FCornerFilter,Filter_TNoise)
%%Returns filtered output of Attenuator
% Filter_TNoise in A/sqrt(Hz)

kf_Filter = FCornerFilter*(Filter_TNoise^2);

TransferFunction = tf(TFnum,TFden);
TransferFunction_ss = ss(TransferFunction);

Filter_Thermal = randn(1,length(Signal))*Filter_TNoise;
[FilterFlicker,~] = f_alpha(length(Signal),kf_Filter,0.5,1);
Filter_Noise = Filter_Thermal + FilterFlicker';
Input = Signal + Filter_Noise;
%Input = Signal;

%s = tf('s')
%TransferFunction = 1/(314.59265e3 + s);

Filtered = lsim(TransferFunction_ss,Input,Time)';
end
