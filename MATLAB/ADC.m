function [DigitalOutput,sample_times] = ADC(sample_cycle_ratio,Full_Scale,num_bits)

samples = sample_cycle_ratio; 
k = linspace(0,1,samples); % sample points
FS = 2^num_bits;

D = zeros(1,samples);

%Storing sampled values in I array

D = round((FS/2).*sin(2*pi*k));
DigitalOutput = (D.*Full_Scale)/(2^num_bits); %scaling to full scale
sample_times = k;

plot(k,DigitalOutput);
end