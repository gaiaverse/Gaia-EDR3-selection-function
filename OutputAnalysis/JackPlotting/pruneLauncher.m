for i = 100:139
target = "../../../CodeOutput/Diagnostic29_PostProcessing/PostProcessing/" + num2str(i) + ".dat";
data = "../../../Data/MainData/" + num2str(i) + ".csv";
output = "../../../Data/PrunedData/" + num2str(i) + ".csv";
pruneData(target,data,output,0.5)
end