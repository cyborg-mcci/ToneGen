%Pure sine wave PSD

% t = linspace(0,2*pi,1200);
% f_in = 1e6;
% i =sin(2*pi*f_in.*t);
% n = floor(log2(length(t)));
% SCR = 1000;
% i = i(1:2^n)';
% fs = 2^n*f_in/53;
% [snr, enob, pot_signal_B, f, PSD] = gs_fresp(i, 2^n, fs,f_in , 1);



%ADC and Stitched Array%

% [DO,ST] = ADC(1000,2,12);
% [SA, STA] = StitchedArray(DO,ST);



%LSIM%

% TF = tf([1],[1 1]);
% t = linspace(-10,10,1001);
% x = 5*ones(1,length(t));
% y = 5*(1-exp(-t));
% 
% y = lsim(TF,x,t);
% plot(t,y);



%Filtering%

% f_in = 1e3;
% [i,t] = gensig('sin',1/f_in);
% 
% R = 1e3;
% C = 1e-7; %fc = 1.6kHz
% TFnum = [1];
% TFden = [R*C 1];
% TF = tf(TFnum, TFden);
% y = Filtering(TFnum, TFden, i, t);
% plot(t,i);
% hold on
% plot(t,y);
% grid on




%Noise%

% t = linspace(0,2*pi,1200);
% f_in = 1e6;
% i =sin(2*pi*f_in.*t);

% r_th = RThermal_Noise(R,T,length(i));

% t = linspace(0,2*pi,16384);
% f_in = 1e2;
% i =sin(2*pi*f_in.*t);
% R = 1e3;
% T = 293.15;
% 
% TF = tf(1,[1 0]);
% 
% kf = 1; 
% wn = RThermal_Noise(R,T,length(i));
% 
% fs = 4e8;
% f = linspace(0,fs/2,length(i)/2);
% 
% fnoise = 2*pi*kf.*randn(1,length(i));
% S = lsim(TF,fnoise,t);
% i_fnoise = sqrt(abs(S));
% 
% fs = f_in*length(i_fnoise)/41;
% 
% [Pxx, w] = periodogram(i_fnoise,[],[],fs,'onesided');
% f = w./(2*pi);
% ixx = sqrt(Pxx);
% Line = sqrt(kf./f);
% grid on
% N = floor(log2(length(ixx)));
% [snr, enob, pot_signal_B, f, PSD] = gs_fresp(ixx(1:2^N), 2^N, fs, f_in,1);


% t = linspace(0,2*pi,16384);         
% f_in = 1e5;
% i =sin(2*pi*f_in.*t);
% 
% fs = f_in*(length(i))/41;
% y = FlickerNoise(1,fs,t);
% % totalnoise = sqrt(th_n.^2 + y.^2);
% N = floor(log2(length(y)));
% 
% [snr, enob, pot_signal_B, f, PSD] = gs_fresp(y(1:2^N), 2^N, fs, f_in,1);

kf = 10;
f = linspace(0,10^6,10000);
S = kf*randn(1,10000)./f;
fnoise = sqrt(abs(S));
fnoise = 20*log10(fnoise);












