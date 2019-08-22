function [DAC_Output_v] = TIA_Converter(Signal,Converter_R,Temperature,Op_Amp_Noise)

global k;
Resistor_Noise = randn(1,length(Signal))*sqrt(4*k*Temperature*Converter_R);
DAC_Output_v = -1*Signal*Converter_R + sqrt(Resistor_Noise.^2 + (Op_Amp_Noise*randn(1,length(Signal))).^2);

end