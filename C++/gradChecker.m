rooter = "Test";
figure(1);
clf;
f = readtable(rooter + "/surfaceMap.dat");

xt = f.Var1;
xs = f.Var2;
L = f.Var3;

gx_true = f.Var4;
gy_true = f.Var5;
gx_manual = f.Var6;
gy_manual = f.Var7;


[~,I] = max(L);

maxPos = [xt(I),xs(I)]

bounds = [-1,1,-1,1];

tri = delaunay(xt,xs);
t = tiledlayout(2,2);

nexttile(t,[2,1]);
trisurf(tri,xt,xs,L,'LineStyle','None');
view(2);
title("Log Likelihood");
xlabel("Nt term"); ylabel("Ns term");
axis(bounds);
nexttile(t);

magnify = 1;
thin = 10;
n = length(xt);

j = 1:nt;
selectorT = (mod(1:n,thin) == 0);
selectorS = (mod(1:n,thin) == 0);


quiver(xt(selectorT),xs(selectorS),gx_true(selectorT),gy_true(selectorS),magnify);
title("Analytical Gradient");
xlabel("Nt term"); ylabel("Ns term");

axis(bounds);
nexttile(t);

quiver(xt(selectorT),xs(selectorS),gx_manual(selectorT),gy_manual(selectorS),magnify);
title("Analytical Gradient");
xlabel("Nt term"); ylabel("Ns term");
axis(bounds);


figure(2);
clf;
t2 = tiledlayout(2,2);
nexttile(t2);
trisurf(tri,xt,xs,gx_manual,'LineStyle','None');
title("Manual dL/dNt");
xlabel("Nt term"); ylabel("Ns term");
view(2);
colorbar;

nexttile(t2);
trisurf(tri,xt,xs,gx_true,'LineStyle','None');
title("Analytical dL/dNt");
xlabel("Nt term"); ylabel("Ns term");
view(2);
colorbar;


nexttile(t2);
trisurf(tri,xt,xs,gy_manual,'LineStyle','None');
title("Manual dL/dNs");
xlabel("Nt term"); ylabel("Ns term");
view(2);
colorbar;

nexttile(t2);
trisurf(tri,xt,xs,gy_true,'LineStyle','None');
title("Analytical dL/dNs");
xlabel("Nt term"); ylabel("Ns term");
view(2);
colorbar;