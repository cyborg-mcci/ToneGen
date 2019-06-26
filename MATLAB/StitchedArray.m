function [Stitched_DArray, Stitched_TArray] = StitchedArray(Dig_Out,Time_Out)
n = 41; %number of replications of array
Stitched_DArray = repmat(Dig_Out,1,n); %stitching digital output
terminator = n*Time_Out(end);
Stitched_TArray = linspace(0,terminator,n*length(Time_Out));

% L = length(Stitched_TArray);
% N = 2^L;
plot(Stitched_TArray,Stitched_DArray);
end
