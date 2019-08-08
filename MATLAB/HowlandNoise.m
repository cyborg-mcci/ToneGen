%%%%Howland Current Source


k = physconst('Boltzmann');
T = 293.15; %%room temp
en = 0.9e-9; %%op-amp input referred voltage noise density
AOL = 10000; %%open-loop gain
ep = 0.05; %%mismatch error
GBWP = 110E6; %%gain-bandwidth product

f = 10e6; %%input signal frequency (worst case)
w = 2*pi*f;

%%Resistor values

R2A = 1e3;
R2B = 1e3;
R2 = R2A + R2B;
R3 = 1e3;
R4 = 1e3;
G = R4./R3;

R1 =(R2)./G; %%balanced resistor bridge

Rx = ((R2A + R2B).^2)./(R1.*R2B.^2) + R4.^2./(R3.*R2B.^2) + (R2A + R2B + R4)./(R2B.^2);
Co = (R3 + R4).*(R1 + R2A + R2B)./((2*pi*GBWP*R3).*(R1 + R2A).*R2B);

RO1 = (R2B.*(R1 + R2A))./(ep*(R2A + R2B)); %%accounting for resistor bridge mismatch
RO2 = (((R1 + R2A).*(R2B))./((R1 + R2A + R2B))).*((AOL*R3 + R3 + R4)./(R3 + R4)); %%accounting for finite open-loop gain
XCo = 1./(w*Co); %%accounting for finite bandwidth

Zo = (1./RO1 + 1./RO2 + 1./XCo).^(-1); %%output impedance

ThermalNoise = sqrt(4*k*T*Rx);
Opampnoise = en*(R1 + R2)./(R1.*R2B);
TotalNoise = sqrt(ThermalNoise.^2 + Opampnoise.^2);

subplot(2,1,1)
semilogx(R2A,TotalNoise)
title('Increasing R2A and R2B together')
xlabel('R2A');
ylabel('i_{noise} (A/\surdHz)');
subplot(2,1,2)
semilogx(R2B,TotalNoise)
xlabel('R2B');
ylabel('i_{noise} (A/\surdHz)');








