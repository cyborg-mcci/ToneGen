function [DigitalOutput] = ADC(sample_cycle_ratio,Full_Scale,num_bits)

Ts = 1;
samples = sample_cycle_ratio*Ts; %MIGHT SCALE BY NUMBER OF CYCLES
k = linspace(0,Ts,samples); % sample points
FS = 2^num_bits;

D = zeros(1,samples);

%Storing sampled values in I array

D = round((FS/2).*sin(2*pi*k));
DigitalOutput = (D.*Full_Scale)/(2^num_bits); %scaling to full scale

plot(k,DigitalOutput);
end