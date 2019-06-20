function [Stitched_Array, time_series] = StitchedArray(sample_cycle_ratio,sample_period,Full_Scale,num_bits)
[D,t1] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
n = 10; %number of replications of array
Stitched_Array = repmat(D,1,n); %stitching digital output
terminator = n*t1(end);
time_series = linspace(0,terminator,n*length(t1));
plot(time_series,Stitched_Array);

end
