% s_n = 1000; %number of samples
% Ts = 1; %cycle
% k = linspace(0,Ts,s_n); % sample points
% N = 3; %number of bits
% I_pp = 2^N; 
% f = 1;
% w = 2*pi*f;
% 
% D_I = zeros(1,s_n);
% LSB = (I_pp)/(2^N);
% 
% y = zeros(1,s_n);
% %Storing sampled values in I array
% for j = 1:s_n
%     D_I(j) = round((I_pp/2)*sin(w*k(j)));
%    
% end
% 
% plot(k,D_I);

%ADC
% [D,t] = ADC(1000,2,5);

n = 3;
A = [1 2 3];
B = repmat(A,1,n)












 


















