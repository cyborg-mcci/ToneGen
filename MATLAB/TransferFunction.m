function [Response] = TransferFunction(TF,Stitched_DArray,Stitched_TArray)

Response = lsim(TF,Stitched_DArray,Stitched_TArray);
plot(Stitched_TArray,Response);

end

