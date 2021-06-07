% clear;
folder = "../../../CodeOutput/Diagnostic29_PostProcessing/PostProcessing/";
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
n = 100:100;
files = string(n)+".dat";
datas = "../../../Data/MainData/" + string(n) + ".csv";
T = table;
ns = [];
ks = [];
times = [];
clf;
Lorig = [];
Lmodif = [];
for i = 1:length(files)
    file = files(i);
    data = datas(i);
    disp("Opening " + file);
    
    [tt,nn,kk,tLine] = analyseFile(file,data,folder,0.01);
    times = [times,tt];
    ns = [ns,nn];
    ks = [ks,kk];
    T = [T;tLine];
%     plotHists(times,ns,ks);

%     f = readtable(folder+file,"ReadVariableNames",true);
%     Lorig = [Lorig;f.OriginalContribution(:)];
%     Lmodif = [Lmodif;f.FlattenedGap(:)];
end
% 
% cut = -5;
% Lorig = Lorig(Lorig > cut);
% Lmodif = Lmodif(Lmodif > cut);
% nBins = 1000;
% subplot(2,1,1);
% hist(Lorig,nBins);
% xlabel("Original L");
% ylabel("Counts");
% subplot(2,1,2);
% hist(Lmodif,nBins);
% xlabel("Smoothed L");
% ylabel("Counts");
% T.Properties.VariableNames = {'File','nStars','nAnomalies','AnomalyPercentage'};
disp(T)
plotHists(times,ns,ks)

function [times,ns,ks,tLine] = analyseFile(file,data,folder,anomalyCriteria)
    times = [];
    ns = [];
    ks = [];
    formats = {'%d','%d'};
    names = ["n","k"];
    for i = 1:500
        formats{end+1} = '%d';
        names(end+1) = "Obs" + num2str(i);
    end
    tf = tabularTextDatastore(folder+file,"ReadVariableNames",true);
    td = tabularTextDatastore(data,"ReadVariableNames",false,"VariableNames",names,"TextscanFormats",formats);
    
%     td.VariableNames=names;
%     td.TextscanFormats = formats;
    
    
%     td.Formats = formats;
    f = tall(tf);
    d = tall(td);
   
   mod = abs(f.FlattenedGap);
   orig = abs(f.OriginalContribution);
   diff = (mod - orig)./orig;
   
   
%    plot(diff);
%    xlim([0,100])
   anomaly = gather(diff > anomalyCriteria);
   
   n = height(anomaly);
   a = sum(anomaly);
   fr = a/n*100;   
   
   tLine = {file,n,a,fr};
   
   
   faultyStars = gather(d(anomaly,:));

   
   step = 5;
   threshold = step;
   for k = 1:a
       ks(end+1) = faultyStars{k,1};
       n = faultyStars{k,2};
       ns(end+1) = n;
       for z = 1:n
           j = faultyStars{k,2+z}; 
           
          times(end+1) = j;
       end
       
       progress= round((k/a)*100);
       if progress >= threshold
          while threshold <= progress
              threshold = threshold + step;
          end
          fprintf("\t"+num2str(progress) + "%%\n");
       end
   end

end

function plotHists(times,ns,ks)
    smallCut = 1e9;
    
    f = readmatrix('../../../Data/MainData/100.csv');
    g = f(:,3:end);
    h = reshape(g,1,[]);
    Nfull = height(f);
    nReduced = length(ns);
    
    minTime = 0e6;
    maxTime = 10e6;
    selector = (times < maxTime) & (times > minTime);
    selectorBase = (h < maxTime) & (h > minTime);
    h = h(selectorBase);
    smallStars = times(selector);
    
    span = (max(smallStars) - min(smallStars));
    nBins = min(span,300);
    baseLine = linspace(min(smallStars),max(smallStars),nBins);
    subplot(2,2,[1,3]);
    
    [normal,edges] = histcounts(h,nBins);
    sBins = histcounts(smallStars,edges);
    
    normalPredicted = normal * nReduced/Nfull;
    
    bar(baseLine,sBins,1);
    
    xlabel("Gaia Timesteps (/10s)");
    ylabel("Counts");
    title("Observations of stars classed as anomalous")
    % set(gca,'yscale','log');
    subplot(2,2,2);
    bar(baseLine,normalPredicted,1);
    title("Predicted Observations Assuming No Temporal Bias")
    xlabel("Gaia Timesteps (/10s)");
    ylabel("Predicted Counts");
    subplot(2,2,4);
    bar(baseLine,(sBins - normalPredicted)./normalPredicted,1)
    title("Bias in Anomalous Stars")
    xlabel("Gaia Timesteps (/10s)");
    ylabel("Fractional Bias");
    drawnow;
end