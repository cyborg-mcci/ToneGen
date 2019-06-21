sample_cycle_ratio = 1000;
Full_Scale = 2;
num_bits = 5;
f_in = 0.5;

[Dig_Out,Time_Out] = ADC(sample_cycle_ratio,Full_Scale,num_bits);
Normalised_time = Time_Out.*(1/f_in); % normalising time series

[Stitched_DArray,Stitched_TArray] = StitchedArray(Dig_Out,Normalised_time);
plot(Stitched_TArray,Stitched_DArray);

fs = sample_cycle_ratio * f_in;
OSR=1;
L = length(Stitched_TArray);
y = log2(L);
N = floor(y);

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(Stitched_TArray, 2^N, fs, fi, OSR);

% [Stitched_Array, time_series] = StitchedArray(Dig_Array,sample_cycle_ratio,Full_Scale,num_bits);