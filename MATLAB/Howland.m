function O_R_Current_Noise_Density = Howland(R1,R2A,R2B,R3,R4,Op_Amp_Noise_Input_Referred)
%Returns output-referred current noise density 
% e_n = input-referred voltage noise density (V/sqrt(Hz))

if R4/R3 ~= (R2A + R2B)/R1
    warning('Must balance the resistor bridge for Howland Current Source');

else
    e_n = Op_Amp_Noise_Input_Referred;
    Rx = ((R2A + R2B)^2)/(R1*R2B^2) + R4^2/(R3*R2B^2) + (R2A + R2B + R4)/(R2B^2);

    ThermalNoise = sqrt(4*k*Temperature*Rx);
    Opampnoise = e_n*(R1 + R2)/(R1*R2B);
    O_R_Current_Noise_Density = sqrt(ThermalNoise^2 + Opampnoise^2);

  
end