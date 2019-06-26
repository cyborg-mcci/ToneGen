%Pure sine wave PSD%%
t = linspace(0,2*pi,1200);
f_in = 1e6;
i =sin(2*pi*f_in.*t);
n = floor(log2(length(t)));
SCR = 1000;
i = i(1:2^n)';
fs = 2^n*f_in/53;
[snr, enob, pot_signal_B, f, PSD] = gs_fresp(i, 2^n, fs,f_in , 1);

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

TFnum = 1;
TFden = [1 1];
t = linspace(0,10,1001);
x = 5*ones(1,length(t));

y = 5*(1 - exp(-t));

R = Filtering(TFnum,TFden,x,t);
plot(t,y)
plot(t,R)






















    














 


















