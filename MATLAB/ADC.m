function [DigitalOutput,sample_times] = ADC(sample_cycle_ratio,Full_Scale,num_bits)

MC=41;
samples = sample_cycle_ratio; 
k = 0:samples:MC-samples; % sample points
FS = 2^num_bits;

D = zeros(1,length(k));

%Storing sampled values in I array

D = round((FS/2).*sin(2*pi*k));
DigitalOutput = (D.*Full_Scale)/(2^num_bits); %scaling to full scale
sample_times = k;

plot(k,DigitalOutput);
end