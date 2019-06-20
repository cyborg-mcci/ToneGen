%%%%%%Fucntion taken from a Boris Murman presentation%%%%%%%%%%

function [inl, dnl]=gs_inldnl(ADC_codes, Nbit, alpha)
FigureCounter = 0;
alpha=alpha;

first_bin = 0;%min(ADC_codes);
last_bin = 2^Nbit-1;%max(ADC_codes);
h=hist(ADC_codes, first_bin:last_bin);

for i=1:1:2^Nbit-1
    if i==1
    dnl(i)=h(i)/(alpha-1) -1;      %%%Accounts for the fact that the ramp
                                   %%%Vector starts at 1 rather than 0.
    else
    dnl(i)=h(i)/alpha -1;          %%%Main DNL formula
    end
    inl(i)=dnl(i)+sum(dnl(1:i));   %%%Main INL formula
end

FigureCounter = FigureCounter + 1;
figure(FigureCounter);
plot(dnl)
grid on
xlabel('Output code','FontSize',14)        %Plotting options
ylabel('DNL, LSB','FontSize',14)              %Plotting options
set(gca,'FontSize',14)  
FigureCounter = FigureCounter + 1;
figure(FigureCounter);
plot(inl)
grid on
xlabel('Output code','FontSize',14)        %Plotting options
ylabel('INL, LSB','FontSize',14)              %Plotting options
set(gca,'FontSize',12)


