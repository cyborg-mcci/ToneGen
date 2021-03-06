function [Spectrum, f_TS, PSD_OSdB, f_OS, Window] = wall_fresp(Signal, t, windowType, PlotFlag)

global FigureCounter

if(~exist('FigureCounter', 'var'))
    FigureCounter = 0;
end
% Calculating the sample frequency
Ts = t(2)-t(1);
fs = 1/Ts;
% Calculating the Sample Length
SampleLength = length(Signal);

% Creating the Window
if(isequal(windowType, @rectwin))
    Window = window(@rectwin, SampleLength);
else
    Window = window(windowType, SampleLength, 'periodic');
end


% Ensuring Vector orientation of window and signal are the same
if size(Signal,1) ~= size(Window,1)
    Window = Window';
end
% Windowing the input signal
WindowedSignal = Signal .* Window;

% Performing the fft to get the raw spectrum
Spectrum = fft(WindowedSignal);

% Accounting for Window gain
Spectrum  = Spectrum ./ sqrt(sum(Window));

% Calculating the 2-sided Frequency Vector
f_TS = (0:SampleLength) * fs/SampleLength;

if(PlotFlag)
    FigureCounter = FigureCounter + 1;
    figure(FigureCounter) 
    clf
    loglog(f_TS, abs(Spectrum))
    grid on
    title('Two Sided FFT')
end
% Creating the One-Sided PSD
OneSidedFFT = Spectrum(1:end/2+1);
PSD_OS = (1/(fs)) * abs(OneSidedFFT).^2;
PSD_OS(2:end-1) = 2 * PSD_OS(2:end-1);
PSD_OSdB = 10*log10(PSD_OS);
f_OS = (0:SampleLength/2) * fs/SampleLength;


% Plotting the One-Sided PSD
if(PlotFlag)
    FigureCounter = FigureCounter + 1;
    figure(FigureCounter)
    clf
    semilogx(f_OS, PSD_OSdB)
    grid on
    title('One-Sided PSD in dB')
end
end
