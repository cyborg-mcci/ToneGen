MC = 41;
M = 2048;

fs = 1000;
fin = fs * (MC/M);

Ts = 1/fs;
Tin = 1/fin;

t = 0:Ts:(MC*Tin)-Ts;

V = sin(2*pi*fin*t);

figure
plot(t, V)

[snr, enob, pot_signal_B, f, PSD] = gs_fresp(V', 2048, fs, fin, 1);
