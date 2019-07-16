function Flicker = FlickerNoise(kf,Frequency_Vector,Vector_Length)

S = kf./Frequency_Vector;
Flicker = randn(1,Vector_Length).*sqrt(abs(S)); 

end