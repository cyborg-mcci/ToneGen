function [Stitched_Array, sample_times] = StitchedArray(sample_cycle_ratio,Full_Scale,num_bits)
[D,t1] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
n = 10; %number of replications of array
Stitched_D_Array = repmat(D,1,n);
sample_times = 

end
