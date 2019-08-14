function [Signal_RC, t_RC] = wall_invfresp(Spectrum, f_TS, Window, PlotFlag)

global FigureCounter

if(~exist('FigureCounter', 'var'))
    FigureCounter = 0;
end

fs = f_TS(2)-f_TS(1);

if(length(Spectrum)~= length(Window))
   error('Length of the Spectrum not equal to Length of the Window. Are you giving one-sided spectrum by accident?')  
end
SampleLength = length(Spectrum);

% Un-Normalising the Data
Spectrum = Spectrum .* sum(Window);

% Inverse Fourier Transform
Signal_RC = ifft(Spectrum);
% Un-Windowing the Data
Signal_RC = Signal_RC ./ Window;

Ts = 1/fs;

t_RC = 0:Ts:(SampleLength-1)*Ts;

% Plotting the Recreated Signal
if(PlotFlag)
    FigureCounter = FigureCounter + 1;
    figure(FigureCounter)
    clf
    plot(t_RC, Signal_RC);
    title('Reconstructed Signal')
end

end
