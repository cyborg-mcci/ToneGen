function [Filtered_Signal] = Filtering(TFnum,TFden,Signal,Time,FCornerFilter,TNoiseFilter,fs,IO_Refer)
%%Returns filtered output 
% TNoiseFilter in A/sqrt(Hz)
%IO_Refer indicates whether the noise is input or output referred
%(0=Input,1=Output), input is default

TransferFunction = tf(TFnum,TFden);


sigma_Thermal = TNoiseFilter*sqrt(fs/2);
Thermal = randn(1,length(Signal))*sigma_Thermal;
var_Flicker = FCornerFilter*(TNoiseFilter^2)*log(length(Signal)/2);

[FlickerFilter,~] = f_alpha(length(Signal),var_Flicker,1,1);

Filter_Noise = Thermal + FlickerFilter';

if IO_Refer == 0
    Input = Signal + Filter_Noise;
    Filtered_Signal = lsim(TransferFunction,Input,Time)';
elseif IO_Refer==1
    Input = Signal;
    Filtered_Signal = lsim(TransferFunction,Input,Time)' + Filter_Noise;
else
    error('Must specify if noise is referred to the input or output');

end
