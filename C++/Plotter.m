set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

files = ["TestBaseNew","TestInversion"];
folder = "TestBaseNew";
% getData(10)


temporalPlot(folder);
progressPlot(files)

function getData(timeGap)
    f = load("SyncTime.mat");
    SyncCurrentTime = datetime('now');
    timeSince = seconds(SyncCurrentTime - f.SyncCurrentTime);

    if timeSince > timeGap
        system(' rsync -avr "jackfraser@hydra.physics.ox.ac.uk:/mnt/extraspace/GaiaSelectionFunction/Code/C++/Output/" Output/');
        save("SyncTime.mat","SyncCurrentTime");
    end
end
function temporalPlot(folder)
    gaps = readtable("Output/edr3_gaps.csv");
    properties = readtable("Output/" + folder + "/Optimiser_Properties.dat");
    z= readmatrix("Output/" + folder + "/TempPositions/TempPosition_TransformedParameters.dat");
    t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
    Nt = properties.Nt(1);
    
    figure(1);
    T = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    nexttile(T);
    f = z(1:Nt);
    m = z(Nt+1:end);
    x = linspace(t(1),t(2),length(f));
    q = 1./(1 + exp(-f));
    z = q;
    plot(x,z,'k');
    hold on;

    for i = 1:height(gaps)
        t1 = floor(gaps.tbeg(i))-1;
        t2 = ceil(gaps.tend(i))+1;
        fill([t1,t1,t2,t2],[0,1,1,0],'b','LineStyle','None','FaceAlpha',0.3);
    end

    hold off;
    legend("pt","Known Gaps");
end
function progressPlot(files)

    figure(2);

    patterns = ["-","-","-","--"];
    cols = colororder;
    cols2 = min(1,cols*1.5);
    minLim = 1;
    clf;
    hold on;
    for i = 1:length(files)
        file = "Output/" + files(i) + "/OptimiserProgress.txt";
        f = readtable(file);
    
        fullEpoch = f(f.Batch == -1,:);
        miniBatches = f(f.Batch > -1,:);
		
		cutx = false(1,height(fullEpoch));
		for j = 2:height(fullEpoch)
			up = fullEpoch.nBatches(j);
			down = fullEpoch.nBatches(j-1);
			if up ~= down
				cutx(j-1) = true;
			end
		end
		shrinkLines = fullEpoch(cutx,:);
		
        subplot(2,1,1);
        hold on;

        miniX = miniBatches.Epoch-1 + miniBatches.Batch ./ miniBatches.nBatches+1e-3;
        L0 = miniBatches.F(1);
        xB = miniX; %miniBatches.Elapsed;
        xF = fullEpoch.Epoch; %fullEpoch.Elapsed;

%         plot(fullEpoch.Elapsed,xF);
        plot(miniBatches.Elapsed,xB,'Color',cols2(i,:),'LineWidth',0.5);   
        plot(fullEpoch.Elapsed,xF,'LineWidth',2,'Color',cols(i,:),'HandleVisibility','Off');
        
		cx = [];
		cy = [];
		cz = [];
		for j = 1:height(shrinkLines)
			cx(end+1) = shrinkLines.Elapsed(j);
			cy(end+1) = shrinkLines.Epoch(j);
			cz(end+1) = shrinkLines.F(j)/L0;
		end
		scatter(cx,cy,40,cols(i,:),'Filled','HandleVisibility','Off');
%         xlim([400,max(xF)])

        ylabel("Complete Epochs");
        xlabel("Time Elapsed (s)");
        grid on;

        set(gca,'yscale','log')
        set(gca,'xscale','log')
        
        subplot(2,1,2);
        hold on;
        plot(miniBatches.Elapsed,miniBatches.F/L0,'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
        
        plot(fullEpoch.Elapsed,fullEpoch.F/L0,'LineWidth',1.4,'Color',cols(i,:));
		scatter(cx,cz,40,cols(i,:),'Filled','HandleVisibility','Off');
        set(gca,'yscale','log')
        set(gca,'xscale','log')
          xlabel("Elapsed Time (s)");
        ylabel("$L/L_0$");
        hold off;
        grid on;
%         hold on;
%         dL = fullEpoch.dF;
%         dL(1) = NaN;
%         set(gca,'yscale','log')
% 
% 
%         dL = dL./fullEpoch.F;
%         dLm = dL;
%         dLp = dL;
%         dLm(dLm > 0) = 1e-8;
%         dLp(dLp < 0) = 1e-8;
%         plot(xF,abs(dLm),'b','LineStyle',patterns(i));
% 
%         plot(xF,dLp,'r','LineStyle',patterns(i));
% 
%         legend('Negative Steps','Positive Steps');
%         set(gca,'yscale','log');
%         set(gca,'xscale','log')
%         xlabel("Elapsed Time");
%         ylabel("$\frac{\Delta L}{L}$","FontSize",15);
%         grid on;

    end
    subplot(2,1,1);
    legend('Baseline Model',"Inverted Transform Order");
    % figure(3);
    % cla;
    % M = 213;
    % for i = [1:M:height(m)-M]
    %    msub = m(i:i+M);
    %    psub = 1./(1 + exp(-msub));
    %    hold on;
    %    plot(psub);
    %    hold off;
    % end

end

