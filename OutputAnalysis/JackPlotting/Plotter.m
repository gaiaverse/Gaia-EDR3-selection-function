set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

files = "Diagnostic76_gapPrior_" + ["alpha4_beta8","alpha4_beta10", "alpha2_beta10", "alpha2_beta8"];
% files = "hometest";
getData(60);
files = files([4]);
N1 =2;
N2 = 96;
gap = 2;
progressPlot(files,0)
% gifPlot(files,N1,N2,gap,"evolution4.gif",false);
temporalPlot(files,N2);



function getData(timeGap)
f = load("SyncTime.mat");
SyncCurrentTime = datetime('now');
timeSince = seconds(SyncCurrentTime - f.SyncCurrentTime);

if timeSince > timeGap
	system(' rsync -avr --exclude "*PostProcessing" "jackfraser@hydra.physics.ox.ac.uk:/mnt/extraspace/GaiaSelectionFunction/Output/" ../../../CodeOutput/');
	
	SyncCurrentTime = datetime('now');
	save("SyncTime.mat","SyncCurrentTime");
end
end
function temporalPlot(folders,number)
figure(1);
clf;
ny = 2;
nx = 2;
t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
xmin = t(1);
xmax = t(2);
% xmin = 2390;
xmax = 1500;
ymin = -16;
ymax = 16;
gaps = readtable("edr3_gaps.csv");
map = colororder;
nc =  [0.83 0.14 0.14;
             1.00 0.54 0.00;
             0.47 0.25 0.80;
             0.25 0.80 0.54];
map = [ map; nc];

while height(map) < length(folders)
	map = [map; [0,0,0]];
end
subplot(ny,nx,3);
hold on;
for j = 1:length(folders)
	plot([-1000,-1000],[0,1],'Color',map(j,:));
end
hold off;
for i = 1:length(folders)
	folder = folders(i);
	properties = readtable("../../../CodeOutput/" + folder + "/OptimiserProperties.dat","ReadRowNames",true,"Delimiter","=");
	pData = table2array(properties)';
    vnames = properties.Properties.RowNames;
    properties = array2table(pData,"VariableNames",vnames);
	name = "../../../CodeOutput/" + folder + "/TempPositions/TempPosition";
		
	if number > -1
		name = name + num2str(number);
	end
	name = name + "_TransformedParameters.dat";
	if number == -1
		name = "../../../CodeOutput/" + folder + "/FinalPosition_TransformedParameters.dat";
	end
	
	
	z= readmatrix(name);
	
	
	
	Nt = properties.Nt(1);
	Nl = properties.Nl(1);
	Nm = properties.Nm(1);
    
    varianceSegment = z(Nt+Nl*Nm+1:end);
    
    pop = properties.NVariancePopulations(1);
    pow = properties.hyperOrder(1);
    fprintf("\nVariance output for " + name+"\n");
    for k = 1:pop
       ps = [];
        for j = 0:pow
           ps(end+1) = varianceSegment(j*pop+k);  
        end
        frac = varianceSegment((1+pow)*pop + k);
        fprintf("\tPop %d has fraction %.8f and variance model ",k,frac)
        for j = 0:pow
            s = "%f ";

            if j > 0
                s = " + " + s + "n^%d";
                powFac = j;
                fprintf(s,varianceSegment(j*pop+k)^powFac,j);  
            else
                 fprintf(s,varianceSegment(j*pop+k));  
            end
           
        end
        fprintf("\n");
    end
    
	%     figure(1);
	%     T = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
	%     nexttile(T);
	f = z(1:Nt);
	

	
	m = z(Nt+1:Nt+Nm*Nl);
	magT = z(Nt+Nm*Nl+1:end);
	x = linspace(t(1),t(2),length(f));
	cutT = (x >= xmin) & (x <= xmax);
	q = 1./(1 + exp(-f));
	z = q;
	
	%     [sx,sz] = bottomOut(x,z,1);
	
	subplot(ny,nx,[1,2]);
	
	hold on;
	plot(x(cutT),z(cutT),'Color',map(i,:),"HandleVisibility","Off");
	hold off;
	subplot(ny,nx,3);
	
	hold on;
	plot(x(cutT),f(cutT),'Color',map(i,:),"HandleVisibility","Off");
	
	subplot(ny,nx,[1,2]);
	frameTitle = "Frame " + num2str(number);
	if number == -1
		frameTitle = "Final Position";
	end
	title("$P_t$ " + frameTitle );
	%          legend("pt","Known Gaps");
	xlim([xmin,xmax])
	ylim([0,1.01])
	xlabel("OBMT (Revolutions)");
	ylabel("Detection Efficiency,$P_t$")
	
	subplot(ny,nx,3)
	title("$x_t$ " + frameTitle);
	xlabel("OBMT (Revolutions)");
	ylabel("Detection Parameter,$x_t$")
	xlim([xmin,xmax])
	ylim([ymin,ymax])
	grid on;
	
	subplot(ny,nx,4);
	
	
	q = zeros(1,Nm);
	
	hold on;
	ms = reshape(m,Nm,Nl);
	
	for j = 1:Nl
		q = q + ms(:,j)';
    end
    q = q';
    maxes = max(ms,[],2);
    
    mins = min(ms,[],2);
	zs = [0:length(mins)-1];
    xRow = [zs, fliplr(zs)]';
    yRow = [maxes; flipud(mins)];

    expMode = true;
    if Nm > 1
        alpha = 0.5*log(2);
        if expMode == true
            p = @(x)  exp(-alpha*2*exp(-x));
           yl = [0,1];
        else
             p = @(x) x;
              yl = [-12,18];
             
        end
		plot(zs,p(q/Nl),'Color',map(i,:),"HandleVisibility","Off")
        
        fill(xRow,p(yRow),map(i,:),'FaceAlpha',0.1,'EdgeColor','None');
%         plot(maxes,'Color',map(i,:),"HandleVisibility","Off")
%         plot(mins,'Color',map(i,:),"HandleVisibility","Off")
	else
		scatter(0,q/Nl,'MarkerEdgeColor',map(i,:),"HandleVisibility","Off");
    end
  
	hold off;
	%     plot(m);
	title("Spatial Components " + frameTitle);
	xlabel("Source file, $i$.csv");
	ylabel("Mean $x_{ml}$ on sky")
	xlim([0,Nm-1])
	ylim(yl)
	grid on;
	
	
end


for i = 1:height(gaps)
	t1 = (gaps.tbeg(i));
	t2 = (gaps.tend(i));
	subplot(ny,nx,[1,2]);
	hold on;
	fill([t1,t1,t2,t2],[0,2,2,0],'b','LineStyle','None','FaceAlpha',0.1,"HandleVisibility","Off");
	hold off;
	subplot(ny,nx,3);
	hold on
	fill([t1,t1,t2,t2],[ymin,ymax,ymax,ymin],'b','LineStyle','None','FaceAlpha',0.1,"HandleVisibility","Off");
	hold off;
	
	%         subplot(ny,nx,4);
	%         hold on
	%         fill([t1,t1,t2,t2],[ymin,ymax,ymax,ymin],'b','LineStyle','None','FaceAlpha',0.3,"HandleVisibility","Off");
	%         hold off;
end

subplot(ny,nx,3);

legend(folders,"Interpreter","None")
end
function progressPlot(files,minLim)

figure(2);

patterns = ["-","-","-","--"];
cols = colororder;

while height(cols) < length(files)
	cols = [cols; [rand,rand,rand]];
end

cols2 = min(1,cols*1.5);
clf;
hold on;
miniSmooth = 1;

map = [];
for i = 1:length(files)
	map = [map; cols2(i,:); 1-cols2(i,:); cols(i,:); 1-cols(i,:)];
end

for i = 1:length(files)
	file = "../../../CodeOutput/" + files(i) + "/OptimiserProgress.txt";
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
			L0 =1;
	
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
    if height(fullEpoch) > 0
		if height(miniBatches) > 0
	        ender = max(max(xB),max(xF));
		else
			ender = max(xF);
		end
    else
        ender = max(xB);
    end
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
		cz1(end+1) = abs(shrinkLines.dX(j)/shrinkLines.F(j));
		cz2(end+1) = shrinkLines.GradNorm(j);
	end
	scatter(cx,cy,40,cols(i,:),'Filled','HandleVisibility','Off');
% 	xlim([minLim,ender])
	
	ylabel("Complete Epochs");
	xlabel("Time Elapsed (s)");
	grid on;
	
	%         set(gca,'yscale','log')
	%         set(gca,'xscale','log')
	
	subplot(2,2,2);
	hold on;
	plot(xB,smooth(miniBatches.F/L0,miniSmooth),'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
	
	plot(xF,fullEpoch.F/L0,'LineWidth',1.4,'Color',cols(i,:));
	scatter(cy,cz,40,cols(i,:),'Filled','HandleVisibility','Off');
	set(gca,'yscale','log')
% 	        set(gca,'xscale','log')
	xlabel("Complete Epochs");
	ylabel("$L/L_0$");
	hold off;
	grid on;
	xlim([minLim,ender])
	
	
	subplot(2,2,3);
    hold on;
	plot(xB,smooth(miniBatches.dX,miniSmooth),'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
	
	plot(xF,fullEpoch.dX,'LineWidth',1.4,'Color',cols(i,:));
	scatter(cy,cz1,40,cols(i,:),'Filled','HandleVisibility','Off');
	hold off;
	
	xlabel("Complete Epochs");
	ylabel("$|\Delta X|$");
	set(gca,'yscale','log')
	grid on;

	xlim([minLim,ender])

	
	subplot(2,2,4);
	hold on;
	plot(xB,smooth(miniBatches.GradNorm,miniSmooth),'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
	
	plot(xF,fullEpoch.GradNorm,'LineWidth',1.4,'Color',cols(i,:));
	scatter(cy,cz2,40,cols(i,:),'Filled','HandleVisibility','Off');
	set(gca,'yscale','log')
	%         set(gca,'xscale','log')
	xlabel("Complete Epochs");
	ylabel("$|\nabla L|$");
	hold off;
	grid on;
	xlim([minLim,ender])
	
end

legend(files)
end
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

