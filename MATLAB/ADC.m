function [DigitalOutput,Normalised_Time] = ADC(cycle_sample_ratio,Full_Scale,num_bits, Nearest_Prime,Phase_Shift)

samples = cycle_sample_ratio; 
k = 0:samples:Nearest_Prime-samples; % sample points
FS = 2^num_bits; % original full-scale to give precise quantisation for digitised sine wave

%Storing sampled values in Digital array
D = round((FS/2).*sin(2*pi*k + Phase_Shift)); %1 Hz sine wave 
DigitalOutput = (D.*Full_Scale)/(FS); %scaling to full scale
Normalised_Time = k;

end