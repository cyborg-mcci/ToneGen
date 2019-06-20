%%%% Coder:   Gerardo Salgado.
%%%% Company: MCCI/Tyndall/UCC.
%%%% Date:    Jun 9th, 2017, 16h50hrs.
%%%% Update:  Aug 12th: Added Periodic to Hann, explicit use of PSD.
%%%%          Aug 14th: Corrected FFT normalization Factor, it is amplitude not RMS.  
%%%% Description of function: It calculates frequency response of periodically
%%%% sampled signal, it plots Power-Spectral-Density and computes SNDR/ENOB
%%%% x = time domain vector in the form [N,1], i.e. [N x 1] matrix.
%%%% N = Number of data point in the analized vector x. Must be power of
%%%% two, i.e, 1024, 2048, 4096...etc.
%%%% fs = the sampling rate of the analized data
%%%% fi = the frequency of the input data x
%%%% OSR= Over-Sampling-Ration, set it to 1 if you are using nyquist rate. 
%%%% Comments: the function is fixed to use hann window.
%%%% Copyright UCC-Tyndall-MCCI, 2017.
function [snr, enob, pot_signal_B, f, PSD] = gs_fresp(x, N, fs, fi, OSR)

%%%%Frequency response of a real signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FigureCounter = 0;
ventana=hann(N,'periodic');                  %Window for the data.   %%%blackmanharris
win_data=x.*ventana;                         %Windowed data.
spectrum=2*fft(win_data)/(sum(ventana));     %FFT Single side.
psd=10*log10(abs(spectrum.^2));              %PowerSpectralDensity.
FigureCounter = FigureCounter + 1;
figure(FigureCounter);
plot([0:N/2-1].*(fs/N), psd(1:N/2)); %Plotting options
xlabel('Frequency, Hz','FontSize',14)        %Plotting options  
ylabel('PSD, dB','FontSize',14)              %Plotting options 
set(gca,'FontSize',14)                       %Plotting options
grid on    
f = [0:N/2-1].*(fs/N);
PSD = psd(1:N/2);
%Plotting options
% figure                                       %Plotting options
% semilogx([0:N/2-1].*(fs/N), psd(1:N/2));     %Plotting options
% xlabel('Frequency, Hz','FontSize',14)        %Plotting options
% ylabel('PSD, dB','FontSize',14)              %Plotting options
% set(gca,'FontSize',14)                       %Plotting options
% grid on                                      %Plotting options

%%%Signal-to-Distortion Ratio (no harmonics suppresion) and ENOB%%%%% 

BW=fs/(2*OSR);                               %Bandwidth determination
bin_signal=round(fi*N/fs);                   %Bin of signal
bin_band=floor(N*BW/fs);                     %Latest Bin in bandwidth   


bin_signal=bin_signal+1; %%To account for DC as the first bin.
aux=1;   %this has to change as different windows are used 3 for blackman and harris
%pot_signal=10*log10(sum(abs(spectrum(bin_signal:bin_signal+2)).^2)); %Signal Power
pot_signal=10*log10(sum(abs(spectrum(bin_signal-aux:bin_signal+aux)).^2)); %Signal Power
pot_signal_B=sum(abs(spectrum(bin_signal-aux:bin_signal+aux)).^2); %Signal Power
pot_noise=10*log10(sum(abs(spectrum(3:bin_band)).^2)-sum(abs(spectrum(bin_signal-aux:bin_signal+aux).^2))); %%Noise power
pot_noise_lab=10*log10(sum(abs(spectrum(3:bin_band)).^2));

snr=pot_signal-pot_noise;                    %SNR calculation
enob=(snr-1.76)/6.02;                        %ENOB calculation   