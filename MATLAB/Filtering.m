function [Filtered_Signal] = Filtering(TFnum,TFden,Signal,Time,FCornerFilter,TNoiseFilter,fs,IO_Refer)
%%%% A function which filters signal with a given transfer function and
%%%% returns the filtered signal

%%%% INPUT PARAMETERS
%%% TFnum - the transfer function numerator coefficients (specified in
% descending powers of s) 
%%% TFden - the transfer function denominator coefficients (specified in
% descending powers of s) 
%%% Signal - the input signal to be filtered 
%%% Time - the signal's time vector
%%% FCornerFIlter - the corner frequency of the filter 1/f noise (unit =
% Hz)
%%% TNoiseFilter - the filter's thermal noise (specified per sqrt(Hz))
%%% fs - DAC sampling frequency
%%% IO_Refer - an indicator to refer the noise to the input/output:
% set IO_Refer = 0: Input-referred noise
% set IO_Refer = 1: Output-referred noise

%%%%OUTPUT PARAMETERS
%%% Filtered_Signal - the filtered output signal 

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
