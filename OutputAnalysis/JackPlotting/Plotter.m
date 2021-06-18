set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

% files = ["Diagnostic61_mScaling","Diagnostic61_nScaling"];
files = ["hometest_noVariance","hometest_noPrior"];
% getData(60);

N1 =60;
N2 = 188;
gap = 10;
progressPlot(files, 10)
% gifPlot(files,N1,N2,gap,"mixed_evolution.gif",false,0,0,213);
temporalPlot(files,N2,100,0,42);



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
function temporalPlot(folders,number,magOffset,mStart,mEnd)
figure(1);
clf;
ny = 2;
nx = 2;
t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
xmin = t(1);
xmax = t(2);
xmin = 2310;%2230;
xmax = 2415;%2248;
ymin = -18;
ymax = 18;
gaps = readtable("edr3_gaps.csv");
map = colororder;
subplot(ny,nx,3);
hold on;
for j = 1:length(folders)
	plot([-1000,-1000],[0,1],'Color',map(j,:));
end
hold off;
for i = 1:length(folders)
	folder = folders(i);
	properties = readtable("../../../CodeOutput/" + folder + "/Optimiser_Properties.dat");
	
	name = "../../../CodeOutput/" + folder + "/TempPositions/TempPosition";
	
	tnumber = number;
	if name == "../../../CodeOutput/InitTest/TempPositions/TempPosition" && number > 92
		tnumber = 92;
	end
	
	if tnumber > -1
		name = name + num2str(tnumber);
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
    pow = 4;
    pop = length(varianceSegment)/(pow+2);
    
    for k = 1:pop
       ps = [];
        for j = 0:pow
           ps(end+1) = varianceSegment(j*pop+k);  
        end
        frac = varianceSegment((1+pow)*pop + k);
        fprintf("Pop %d has fraction %.3f and variance model ",k,frac)
        for j = 0:pow
            s = "%f n^%d ";
            if j < pow
                s = s + " + ";
            end
           fprintf(s,varianceSegment(j*pop+k),j);  
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
              yl = [-2,8];
             
        end
		plot(p(q/Nl),'Color',map(i,:),"HandleVisibility","Off")
        
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
	xlim([0,Nm])
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
cols2 = min(1,cols*1.5);
clf;
hold on;
miniSmooth = 10;

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
	% 		L0 =1;
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
	plot(miniBatches.Elapsed,smooth(miniBatches.F/L0,miniSmooth),'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
	
	plot(fullEpoch.Elapsed,fullEpoch.F/L0,'LineWidth',1.4,'Color',cols(i,:));
	scatter(cx,cz,40,cols(i,:),'Filled','HandleVisibility','Off');
% 	set(gca,'yscale','log')
	%         set(gca,'xscale','log')
	xlabel("Elapsed Time (s)");
	ylabel("$L/L_0$");
	hold off;
	grid on;
	xlim([minLim,ender])
	
	
	subplot(2,2,3);
	
	colormap(map);
	hold on;
	
	if height(miniBatches) > 1
		miniBatches.dF(1) = 0;
		x1 = miniBatches.Elapsed';
		y1 = smooth(miniBatches.dF./miniBatches.F,miniSmooth)';
		zq = zeros(size(x1));
		s1 = abs(y1(abs(y1) > 0 & miniBatches.Elapsed' > minLim));
		minner = min(s1);
		maxer = max(s1);
		col = (y1>0)*(2*i);
		
		
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
		col = (y2>0)*(2*i-1);
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
	plot(miniBatches.Elapsed,smooth(miniBatches.GradNorm,miniSmooth),'Color',cols2(i,:),'LineWidth',0.5,'HandleVisibility','Off');
	
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

legend(files)
end

function magComparison(files,numbers,mStart,mEnd,gap,fileName)
figure(1);
number = numbers(1);
frameTitle = "Frame " + num2str(number);
if number == -1
	frameTitle = "Converged Position";
end
ny = 2;
nx = ceil(gap/ny);
for m = mStart:gap:mEnd
	clf;
	
	xmin = 2320;
	xmax = 2326;
	ymin = 0;
	ymax = 1;
	
	for i = 1:length(files)
		
		
		folder = files(i);
		% 					subplot(ny,nx,i);
		properties = readtable("../../../CodeOutput/" + folder + "/Optimiser_Properties.dat");
		
		name = "../../../CodeOutput/" + folder + "/TempPositions/TempPosition";
		number = numbers(i);
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
		
		magT = z(Nt+Nm*Nl+1:end);
		for j = 0:gap-1
			
			fT = frameTitle + " (" + string(m+j) + ".csv)";
			if m + j < 213
				
				subplot(ny,nx,j+1);
				
				magnitudeFrame(magT,Nm,m+j,m+j+1,xmin,xmax,fT);
			end
		end
	end
	frame = getframe(gcf);
	im = frame2im(frame);
	[imind,cm] = rgb2ind(im,256);
	
	if m == mStart
		imwrite(imind,cm,fileName,'gif', 'Loopcount',inf);
	else
		imwrite(imind,cm,fileName,'gif','WriteMode','append');
	end
end
end

function magnitudeFrame(magT,Nm,mStart,mEnd,xmin,xmax,frameTitle)
Ntm = length(magT)/Nm;
t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
if Ntm > 0
	Nts = reshape(magT,Nm,Ntm);
	
	Ntms = linspace(t(1),t(2),Ntm);
	hold on;
	
	for mm = mStart+1:mEnd
		if Ntm > 1
			
			plot(Ntms,1./(1+exp(-Nts(mm,:))));
		else
			
			plot([t(1),t(2)],[1,1]*Nts(mm));
		end
	end
	%             leg = string([mStart]+magOffset) + ".csv";
	legend("Small Data" ," Big Data","No Batch");
	
	hold off;
	xlim([xmin,xmax])
	ylim([0,1])
	
	title("Temporal-Magnitude Component " + frameTitle);
	ylabel("Magnitude Parameter, $x_{mt}$");
	xlabel("OBMT (Revolutions)")
	
	grid on;
end

gaps = readtable("edr3_gaps.csv");
for i = 1:height(gaps)
	t1 = (gaps.tbeg(i));
	t2 = (gaps.tend(i));
	hold on
	fill([t1,t1,t2,t2],[0,1,1,0],'b','LineStyle','None','FaceAlpha',0.15,"HandleVisibility","Off");
	hold off;
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
function gifPlot(folder,startN,maxN,gap,fileName,includeFinal,magOffset,mStart,mEnd)
for i = startN:gap:maxN
	temporalPlot(folder,i,magOffset,mStart,mEnd);
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

function magGif(folder,number,offset,start,max,jump,fileName)
figure(1);
ny = 2;
nx = ceil(jump/ny);
for m = start:jump:max
	
	for i = 1:jump
		if start + i < 213
			subplot(ny,nx,i)
			temporalPlot(folder,number,offset,start +i,start+i);
			%         subplot(2,1,1)
			
			frame = getframe(gcf);
			im = frame2im(frame);
			[imind,cm] = rgb2ind(im,256);
			
			if m == start
				imwrite(imind,cm,fileName,'gif', 'Loopcount',inf);
			else
				imwrite(imind,cm,fileName,'gif','WriteMode','append');
			end
		end
	end
end

end
