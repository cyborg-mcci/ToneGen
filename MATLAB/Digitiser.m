function [DigitalOutput,Normalised_Time] = Digitiser(cycle_sample_ratio,Full_Scale,num_bits, NearestPrime,Phase_Shift)
%%%% A function which outputs a digitised sine wave and a normalised time
%%%% vector

%%%% INPUT PARAMETERS
%%% cycle_sample_ratio - ratio of number of cycles to number of sample
% points
%%% Full_Scale - - peak to peak amplitude of the sine wave (A)
%%% num_bits - the DAC bit precision
%%% NearestPrime - The nearest prime number of cycles governed by the
% coherent sampling condition
%%% Phase shift of the sine wave (units = rads)

%%%%OUTPUT PARAMETERS
%%% DigitalOutput - A digitised sine wave (unit = A)
%%% Normalised_Time = Normalised Time Vector


samples = cycle_sample_ratio; 
k = 0:samples:NearestPrime-samples; % sample points
FS = 2^num_bits; % original full-scale to give precise quantisation for digitised sine wave

%Storing sampled values in Digital array
D = round((FS/2).*sin(2*pi*k + Phase_Shift)); %1 Hz sine wave 
DigitalOutput = (D.*Full_Scale)/(FS); %scaling to full scale
Normalised_Time = k;

end