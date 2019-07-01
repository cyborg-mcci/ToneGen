function [Filtered] = Filtering(num,den,OutputArray,TimeArray)

TransferFunction = tf(num,den);
Filtered = lsim(TransferFunction,OutputArray,TimeArray);

end
