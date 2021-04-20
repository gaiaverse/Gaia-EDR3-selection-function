rooter = "Test2";
figure(2);
clf;
f = readtable(rooter + "/surfaceMap.dat");
global bounds;
xt = f.Var1;
xs = f.Var2;
L = f.Var3;

gx_true = f.Var4;
gy_true = f.Var5;
gx_manual = f.Var6;
gy_manual = f.Var7;

cutNans = isnan(xt) | isnan(xs);

xt = xt(~cutNans);
xs = xs(~cutNans);
L = L(~cutNans);
gx_true = gx_true(~cutNans);
gx_manual = gx_manual(~cutNans);
gy_true = gy_true(~cutNans);
gy_manual = gy_manual(~cutNans);

[~,I] = max(L);

maxPos = [xt(I),xs(I)]

bounds = [0,2,0,2];

tri = delaunay(xt,xs);
t = tiledlayout(2,2);

nexttile(t,[2,1]);
trisurf(tri,xt,xs,L,'LineStyle','None');
view(2);
title("Log Likelihood");
xlabel("Nt term"); ylabel("Ns term");
axis(bounds);
nexttile(t);

magnify = 5;
thin = 1;
n = length(xt);

j = 1:n;
selectorT = (mod(1:n,thin) == 0);
selectorS = (mod(1:n,thin) == 0);


quiver(xt(selectorT),xs(selectorS),gx_true(selectorT),gy_true(selectorS),magnify);
title("Analytical Gradient");
xlabel("Nt term"); ylabel("Ns term");

axis(bounds);
nexttile(t);

quiver(xt(selectorT),xs(selectorS),gx_manual(selectorT),gy_manual(selectorS),magnify,'r');
title("Analytical Gradient");
xlabel("Nt term"); ylabel("Ns term");
axis(bounds);

global root
root = [1.02758, 0.973634];

figure(1);
clf;
t2 = tiledlayout(2,2);

convPlot(tri,xt,xs,gx_manual,gy_manual,"Manual");

convPlot(tri,xt,xs,gx_true,gy_true,"Analytical");



function convPlot(tri,x,y,z1,z2,type)
	global root
	[x1,y1] = findZeroLine(x,y,z1);
	[x2,y2] = findZeroLine(x,y,z2);
	
	
	
	triPlot(tri,x,y,z1,type + " dL/dNt");
	r = sqrt(2);
	maxer = ones(length(x1),1) * (1 + max(max(z1,z2)));
	hold on;
	plot3(x1,y1,maxer,'r');
	plot3(x2,y2,maxer,'k');
	theta = linspace(0,pi/2,length(x1));
	plot3(r*cos(theta),r*sin(theta),maxer,'y');
	
	scatter3(root(1),root(2),maxer(1),'p','MarkerEdgeColor','m','MarkerFaceColor','m');
	hold off;
	
	triPlot(tri,x,y,z2, type + " dL/dNs");

	hold on;
	plot3(x1,y1,maxer,'r');
	plot3(x2,y2,maxer,'k');
	plot3(r*cos(theta),r*sin(theta),maxer,'y');
	scatter3(root(1),root(2),maxer(1),'p','MarkerEdgeColor','m','MarkerFaceColor','m');
	hold off;
end

function triPlot(tri,x,y,z,tit)
	nexttile;
	global bounds;
	trisurf(tri,x,y,z,'LineStyle','None');
	title(tit);
	xlabel("Nt term"); ylabel("Ns term");
	view(2);
	axis(bounds);

	colorbar;
end

function [xu,yu] = findZeroLine(x,y,grid)
	
	xus = unique(x);
	yus = unique(y);
	nx = length(xus);
	ny = length(yus);
	dy = yus(2) - yus(1);
	xu = [];
	yu = [];
	for i = 1:nx-1
		sIdx = (i-1)*ny+1;
		selector = sIdx:sIdx + ny-1;
		line = grid(selector)';
		if any(sign(line) ~= sign(line(1)))
		
			q = find(sign(line) ~= sign(line(1)));

			Idx = q(1);

			if Idx > 1
				zUp = line(Idx);
				zDown = line(Idx - 1);
				
				Delta = dy/(zUp - zDown) * zDown;
				yv = yus(Idx-1) - Delta;
			else

				yv = line(Idx);
			end
			xu(end+1) = xus(i);
			yu(end+1) = yv;

		end
	end


end