function RThermal_Noise = RThermal_Noise(Resistance,Temperature,vector_length)

k = physconst('Boltzmann');
noise_ampl = sqrt(4*k*Temperature/Resistance);
RThermal_Noise = noise_ampl*randn(1,vector_length);

end