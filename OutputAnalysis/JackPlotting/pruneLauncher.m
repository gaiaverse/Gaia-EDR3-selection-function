matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
for i = 100:139
target = "../../../Output/Diagnostic29_PostProcessing/PostProcessing/" + num2str(i) + ".dat";
data = "../../../Data/MainData/" + num2str(i) + ".csv";
output = "../../../Data/PrunedData/" + num2str(i) + ".csv";
pruneData(target,data,output,0.1)
end

for i = 100:100
    shuffleCommand = "shuf ../../../Data/PrunedData/" + num2str(i) + ".csv.temp > ../../../Data/PrunedData/" + num2str(i) + ".csv";
    system(shuffleCommand);
end

system("rm ../../../Data/PrunedData/*.temp");