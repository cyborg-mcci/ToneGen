function [Filtered] = Filtering(TFnum,TFden,Signal,Time,FCornerFilter,TNoiseFilter,fs)
%%Returns filtered output of Attenuator
% Filter_TNoise in A/sqrt(Hz)

TransferFunction = tf(TFnum,TFden);


sigma_Thermal = TNoiseFilter*sqrt(fs/2);
Thermal = randn(1,length(Signal))*sigma_Thermal;
var_Flicker = FCornerFilter*(TNoiseFilter^2)*log(length(Signal)/2);

[FlickerFilter,~] = f_alpha(length(Signal),var_Flicker,1,1);

Filter_Noise = Thermal + FlickerFilter';
Input = Signal + Filter_Noise;
%Input = Signal;

%s = tf('s')
%TransferFunction = 1/(314.59265e3 + s);

Filtered = lsim(TransferFunction,Input,Time)';
end
