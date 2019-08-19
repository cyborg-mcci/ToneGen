Z = [0 2^18 0 2^22 0 2^26];

for i = 2:2:6

    k = physconst('Boltzmann'); %Boltzmann's constant
    Temperature = 300; %Temperature in Kelvin
    Converter_R = 1e3;

    % DAC parameters (ADI LTC1668)
    FCornerDAC = 1e4; %estimate; not in datasheet
    TNoiseDAC_i = 50e-12; %per sqrt(Hz)
    FullScale_i = 1e-3; 
    num_bits = 16;
    fs = 50e6;

    % Input signal parameters and coherent sampling condition
    f_in_desired= 50e3;
    N = Z(i);
    NearestPrime = max(primes(f_in_desired*(N)/fs));
    f_in = fs*NearestPrime/N;

    % Sample RLC filter at the DAC output: parameters
    f_f1 = f_in; % resonant frequency
    BW_f1 = 50; %bandwidth
    L_f1 = 1e-3; %inductance 
    R_f1 = 2*pi*BW_f1*L_f1; %Resistance
    C_f1 = 1/(4*(pi*pi)*(f_f1*f_f1)*L_f1); %Capacitance
    num_f1 = [R_f1/L_f1 0]; %Transfer function numerator coefficients, in descending powers of s
    den_f1 = [1 R_f1/L_f1 1/(L_f1*C_f1)]; %Transfer function denominator coefficients, in descending powers of s
    filter_1 = tf(num_f1,den_f1); %Filter transfer function
    FCornerFilter_1 = 0; %Filter 1/f noise corner frequency
    Filter_TNoise_1_v = sqrt(4*k*Temperature*R_f1); %Filter thermal noise

    %Model

    % DAC Output, includes Converter_R noise 
    [DAC_Output_1_i,DAC_NormalisedTime_1] = DAC(FCornerDAC,FullScale_i,NearestPrime,num_bits,N,TNoiseDAC_i,Converter_R,Temperature);
    [DAC_Output_2_i,DAC_NormalisedTime_2] = DAC(FCornerDAC,FullScale_i,NearestPrime,num_bits,N,TNoiseDAC_i,Converter_R,Temperature);

    % Scaling time vector depending on the input frequency
    Time = DAC_NormalisedTime_1/f_in;

    % Creating two voltage channels
    Channel_1_v = DAC_Output_1_i*Converter_R;
    Channel_2_v = DAC_Output_2_i*Converter_R;

    % Filter output
    Filter_Output_1_v = Filtering(num_f1,den_f1,Channel_1_v,Time,FCornerFilter_1,Filter_TNoise_1_v);
    Filter_Output_2_v = Filtering(num_f1,den_f1,Channel_2_v,Time,FCornerFilter_1,Filter_TNoise_1_v);

    [Filter_Output_1_v_Spectrum, Filter_Output_1_v_f_TS, PSD_Filter_Output_1_v, Filter_Output_1_v_f_OS, Filter_Output_1_v_Window] = wall_fresp(Filter_Output_1_v, Time, @rectwin, 0);
    % [Filter_Output_2_v_Spectrum, Filter_Output_2_v_f_TS, PSD_Filter_Output_2_v, Filter_Output_2_v_f_OS, Filter_Output_2_v_Window] = wall_fresp(Filter_Output_2_v, Time, @rectwin, 0);

    subplot(3,2,i-1)
    plot(Time,Filter_Output_1_v)
    str2 = sprintf('Time response for N =%d',N);
    title(str2)
    xlabel('Time')
    ylabel('Filter Output')
    subplot(3,2,i)
    semilogx(Filter_Output_1_v_f_OS,PSD_Filter_Output_1_v)
    str2 = sprintf('Time response for N =%d',N);
    title(str2)
end