function [outputArg1,outputArg2] = pruneData(pathToPrune,pathToData,pathToOutput,anomalyCriteria)
    formats = {'%d','%d'};
    names = ["k","n"];
    for i = 1:500
        formats{end+1} = '%d';
        names(end+1) = "Obs" + num2str(i);
    end
    tf = tabularTextDatastore(pathToPrune,"ReadVariableNames",true);
    td = tabularTextDatastore(pathToData,"ReadVariableNames",false,"VariableNames",names,"TextscanFormats",formats);

    f = tall(tf);
    d = tall(td);
   
   mod = abs(f.FlattenedGap);
   orig = abs(f.OriginalContribution);
   diff = (mod - orig)./orig;

   anomaly = gather(diff > anomalyCriteria);
   
   a = sum(anomaly);
   n = size(anomaly,1);
   fprintf("Loaded in file with %d lines, found %d anomalies " + pathToOutput + " should have %d stars in it\n",n,a,n-a);
   ids = [1:n];
   badIDs = ids(anomaly);
   
   createFileCommand = "cp " + pathToData + " " +pathToOutput + ".temp";
   system(createFileCommand);
   
   deleteLinesCommand = "sed -i -e '";
   for i = 1:a
       deleteLinesCommand = deleteLinesCommand + num2str(badIDs(i)) + "d";
       if i < a
           deleteLinesCommand = deleteLinesCommand + ";";
       end
   end
   deleteLinesCommand = deleteLinesCommand + "' " + pathToOutput + ".temp";
   system(deleteLinesCommand);
   
   shuffleCommand = "shuf " + pathToOutput + ".temp >" + pathToOutput;
   system(shuffleCommand);
   
   removeTemp = "rm " + pathToOutput + ".temp";
%    system(removeTemp);

end

