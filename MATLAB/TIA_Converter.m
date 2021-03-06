function [DAC_Output_v] = TIA_Converter(Signal,Converter_R,Temperature,Op_Amp_Noise)
%%%% A function which converts the current output of the DAC to a voltage
%%%% signal

%%%% INPUT PARAMETERS
%%% Signal - the current output of the DAC (unit = A)
%%% Converter_R - the TIA resistor (unit = Ohms)
%%% Temperature - the temperature (Unit = K)
%%% Op_Amp_Noise - input-referred voltage noise density (V/sqrt(Hz))

%%%%OUTPUT PARAMETERS
%%% DAC_Output_v - the DAC output converted to a voltage 

global k;
Resistor_Noise = randn(1,length(Signal))*sqrt(4*k*Temperature*Converter_R);
DAC_Output_v = -1*Signal*Converter_R + sqrt(Resistor_Noise.^2 + (Op_Amp_Noise*randn(1,length(Signal))).^2);

end