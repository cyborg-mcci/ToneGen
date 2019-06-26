function [Filtered] = Filtering(num,den,Stitched_DArray,Stitched_TArray)

%might try to take TF as input to function instead if num and den
TransferFunction = tf(num,den);
Filtered = lsim(TransferFunction,Stitched_DArray,Stitched_TArray);

end
