set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

files = ["temptest_tinysigma"];
folder = files(1);
% getData(60);

N1 =0;
N2 = 2780;
gap = 50;
progressPlot(files,0)
% gifPlot(folder,N1,N2,gap,"mum1_evolution.gif",false);
temporalPlot(folder,N2);


function gifPlot(folder,startN,maxN,gap,fileName,includeFinal)
    for i = startN:gap:maxN
        temporalPlot(folder,i);
%         subplot(2,1,1)
       
        frame = getframe(gcf); 
          im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
        if i == startN
          imwrite(imind,cm,fileName,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,fileName,'gif','WriteMode','append'); 
      end 
    end
    
    if includeFinal == true
        temporalPlot(folder,-1);
%         subplot(2,1,1)
       
        frame = getframe(gcf); 
          im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
        if i == startN
          imwrite(imind,cm,fileName,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,fileName,'gif','WriteMode','append'); 
      end 
    end
end

function getData(timeGap)
    f = load("SyncTime.mat");
    SyncCurrentTime = datetime('now');
    timeSince = seconds(SyncCurrentTime - f.SyncCurrentTime);

    if timeSince > timeGap
        system(' rsync -avr "jackfraser@hydra.physics.ox.ac.uk:/mnt/extraspace/GaiaSelectionFunction/Code/C++/Output/" Output/');
        
        SyncCurrentTime = datetime('now');
        save("SyncTime.mat","SyncCurrentTime");
    end
end
function temporalPlot(folder,number)
    figure(1);
    clf;
    ny = 3;
    nx = 1;
    t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
    xmin = t(1);%2310;
    xmax = 2250;%2340;
    ymin = -10;
    ymax = 11.5;
    gaps = readtable("Output/edr3_gaps.csv");
    properties = readtable("Output/" + folder + "/Optimiser_Properties.dat");
    
    name = "Output/" + folder + "/TempPositions/TempPosition";
    if number > -1
        name = name + num2str(number);
    end
    name = name + "_TransformedParameters.dat";
    if number == -1
        name = "Output/" + folder + "/FinalPosition_TransformedParameters.dat";
    end
    
    
    z= readmatrix(name);
    
    
    
    Nt = properties.Nt(1);
    Nl = properties.Nl(1);
    Nm = properties.Nm(1);
%     figure(1);
%     T = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
%     nexttile(T);
    f = z(1:Nt);
    m = z(Nt+1:end);
    x = linspace(t(1),t(2),length(f));
    cutT = (x > xmin) & (x < xmax);
    q = 1./(1 + exp(-f));
    z = q;
    
%     [sx,sz] = bottomOut(x,z,1);
    
    subplot(ny,nx,1);
    plot(x(cutT),z(cutT),'k');
    hold on;
    subplot(ny,nx,2);
    plot(x(cutT),f(cutT),'k');
    hold on;

    for i = 1:height(gaps)
        t1 = (gaps.tbeg(i));
        t2 = (gaps.tend(i));
        subplot(ny,nx,1);
        fill([t1,t1,t2,t2],[0,2,2,0],'b','LineStyle','None','FaceAlpha',0.3);
        subplot(ny,nx,2);
        fill([t1,t1,t2,t2],[ymin,ymax,ymax,ymin],'b','LineStyle','None','FaceAlpha',0.3);
    end
    
    hold off;
    
    subplot(ny,nx,1);
    frameTitle = "Frame " + num2str(number);
    if number == -1
        frameTitle = "Final Position";
    end
     title("$P_t$ " + frameTitle );
     legend("pt","Known Gaps");
     xlim([xmin,xmax])
     ylim([0,1.01])
    xlabel("OBMT (Revolutions)");
    ylabel("Detection Efficiency,$P_t$")
        
     subplot(ny,nx,2)
     title("$x_t$ " + frameTitle);
    xlabel("OBMT (Revolutions)");
    ylabel("Detection Parameter,$x_t$")
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    grid on;
    
    subplot(ny,nx,3);
    
    
    q = zeros(1,Nm);
    
    hold on;
    ms = reshape(m,Nm,Nl);
    for j = 1:Nl
       
        q = q + ms(:,j);
    end
    if Nm > 1
        plot(q/Nl)
    else
        scatter(1,q/Nl);
    end
    hold off;
%     plot(m);
    title("Spatial Components " + frameTitle);
    xlabel("$i$");
    ylabel("Magnitude Bin")
%     xlim([xmin,xmax])
%     ylim([8,11])
    grid on;
    
end
function progressPlot(files,minLim)

    figure(2);

    patterns = ["-","-","-","--"];
    cols = colororder;
    cols2 = min(1,cols*1.5);
    clf;
    hold on;
    for i = 1:length(files)
        file = "Output/" + files(i) + "/OptimiserProgress.txt";
        f = readtable(file);
    
        fullEpoch = f(f.Batch == -1,:);
        miniBatches = f(f.Batch > -1,:);
        fE = 0;
        mE = 0;
         L0 = 0;
        if height(fullEpoch) > 0
            fE = fullEpoch.Elapsed(end);
             L0 = fullEpoch.F(1);
        end
        if height(miniBatches) > 1
            mE = miniBatches.Elapsed(end);
             L0 = miniBatches.F(1);
        end
        ender = max(fE,mE);
		cutx = false(1,height(fullEpoch));
		for j = 2:height(fullEpoch)
			up = fullEpoch.nBatches(j);
			down = fullEpoch.nBatches(j-1);
			if up ~= down
				cutx(j-1) = true;
			end
		end
		shrinkLines = fullEpoch(cutx,:);
		
        subplot(2,2,1);
        hold on;

        miniX = miniBatches.Epoch + (miniBatches.Batch +1)./ miniBatches.nBatches-1;
       
        xB = miniX; %miniBatches.Elapsed;
        xF = fullEpoch.Epoch-1; %fullEpoch.Elapsed;

%         plot(fullEpoch.Elapsed,xF);
        plot(miniBatches.Elapsed,xB,'Color',cols2(i,:),'LineWidth',0.5);   
        plot(fullEpoch.Elapsed,xF,'LineWidth',2,'Color',cols(i,:),'HandleVisibility','Off');
        
		cx = [];
		cy = [];
		cz = [];
        cz1 = [];
        cz2 = [];
		for j = 1:height(shrinkLines)
			cx(end+1) = shrinkLines.Elapsed(j);
			cy(end+1) = shrinkLines.Epoch(j)-1;
			cz(end+1) = shrinkLines.F(j)/L0;
            cz1(end+1) = abs(shrinkLines.dF(j)/shrinkLines.F(j));
            cz2(end+1) = shrinkLines.GradNorm(j);
		end
		scatter(cx,cy,40,cols(i,:),'Filled','HandleVisibility','Off');
        xlim([minLim,ender])
       
        ylabel("Complete Epochs");
        xlabel("Time Elapsed (s)");
        grid on;

%         set(gca,'yscale','log')
%         set(gca,'xscale','log')
        
        subplot(2,2,2);
        hold on;
        plot(miniBatches.Elapsed,miniBatches.F/L0,'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
        
        plot(fullEpoch.Elapsed,fullEpoch.F/L0,'LineWidth',1.4,'Color',cols(i,:));
		scatter(cx,cz,40,cols(i,:),'Filled','HandleVisibility','Off');
        set(gca,'yscale','log')
%         set(gca,'xscale','log')
          xlabel("Elapsed Time (s)");
        ylabel("$L/L_0$");
        hold off;
        grid on;
        xlim([minLim,ender])
        
        
        subplot(2,2,3);
        map = [cols2(i,:); 1-cols2(i,:); cols(i,:); 1-cols(i,:)];
        colormap(map);
        hold on;
        
        if height(miniBatches) > 1
            miniBatches.dF(1) = 0;
            x1 = miniBatches.Elapsed';
            y1 = (miniBatches.dF./miniBatches.F)';
            zq = zeros(size(x1));
            s1 = abs(y1(abs(y1) > 0 & miniBatches.Elapsed' > minLim));
            minner = min(s1);
            maxer = max(s1);
            col = (y1>0)*1.0;


            surface([x1;x1],[abs(y1);abs(y1)],[zq;zq],[col;col],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',0.5);
            hold on;
        end
        if height(fullEpoch) > 1

            fullEpoch.dF(1) = fullEpoch.F(1) - L0;
            x2 = fullEpoch.Elapsed';
            y2 = (fullEpoch.dF./fullEpoch.F)';
            zq = zeros(size(x2));
            col = (y2>0)*1.0+2;
            surface([x2;x2],[abs(y2);abs(y2)],[zq;zq],[col;col],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);
        
             s2 = abs(y2(abs(y2) > 0 & fullEpoch.Elapsed' > minLim));

             if height(miniBatches) > 1
      
                minner = min(min(s1),min(s2));
                maxer = max(max(s1),max(s2));
             else
                 minner = min(s2);
                 maxer = max(s2);
             end
        end
        caxis([0,3]);
		scatter(cx,cz1,40,cols(i,:),'Filled','HandleVisibility','Off');
        hold off;
                xlabel("Elapsed Time (s)");
        ylabel("$\Delta L / L $");
        set(gca,'yscale','log')
%         set(gca,'xscale','log')

        hold off;
        grid on;
        xlim([minLim,ender])
   
        ylim([minner,maxer]);
       
       
        grid on;
         xlim([minLim,ender])
%         set(gca,'yscale','log')
        
        subplot(2,2,4);
        hold on;
        plot(miniBatches.Elapsed,miniBatches.GradNorm,'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
        
        plot(fullEpoch.Elapsed,fullEpoch.GradNorm,'LineWidth',1.4,'Color',cols(i,:));
		scatter(cx,cz2,40,cols(i,:),'Filled','HandleVisibility','Off');
        set(gca,'yscale','log')
%         set(gca,'xscale','log')
          xlabel("Elapsed Time (s)");
        ylabel("$|\nabla L|$");
        hold off;
        grid on;
        xlim([minLim,ender])

    end
end

function [sx,sy] = bottomOut(x,y,factor)
    sx = [];
    sy = [];
    
    i = 1;
    
    while i < length(x)
        top = min(length(x),i+factor);
        sx(end+1) = mean(x(i:top));
        sy(end+1) = min(y(i:top));
        i = i + factor;
    end
end
