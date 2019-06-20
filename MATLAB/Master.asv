sample_cycle_ratio = 1000;
Full_Scale = 2;
num_bits = 5;

[D,t1] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
sample_period = t1(2) - t1(1);
f_in = 1/sample_cycle_ratio;




[Stitched_Array, time_series] = StitchedArray(sample_cycle_ratio,sample_period,Full_Scale,num_bits);