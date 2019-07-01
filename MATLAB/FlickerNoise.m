function FlickerNoise = FlickerNoise(kf,Frequency_Range,Vector_Length)

S = kf*randn(1,Vector_Length)./Frequency_Range;
FlickerNoise = sqrt(abs(S));

end