function [DigitalOutput] = ADC(sample_cycle_ratio,Full_Scale,num_bits)

Ts = 1;
samples = sample_cycle_ratio*Ts; %MIGHT SCALE BY NUMBER OF CYCLES
k = linspace(0,Ts,samples); % sample points
Full_Scale = 2^num_bits; 
f = 1;
w = 2*pi*f;

DigitalOutput = zeros(1,samples);

y = zeros(1,samples);
%Storing sampled values in I array
for j = 1:samples
    DigitalOutput(j) = round((Full_Scale/2)*sin(w*k(j)));
    y(j) = (Full_Scale/2)*(sin(w*k(j)));
end

plot(k,y);
hold on
plot(k,DigitalOutput);
end