f = readtable("gradientTest.txt");

m = (f.x0);

t = tiledlayout('flow');
clf;
nexttile;
scatter(m,f.PriorMu);
%syscale_symlog

[~,i] = max(f.PriorMu);
mMax = m(i)
grid on;


nexttile;
zerOff = 1e-10;
h = plot(m,f.AnalyticalGrad+zerOff);
dtTemplate = get(h, 'DataTipTemplate');
dtTemplate.DataTipRows(1).Format = '%e'; % adjust
dtTemplate.DataTipRows(2).Format = '%e';
hold on;
plot(m,f.NumericalGrad);
hold off;
legend("Analytical", "Numerical");
%`  `yscale_symlog
grid on;

% nexttile;
% diff = (f.AnalyticalGrad - f.NumericalGrad)./(f.AnalyticalGrad + zerOff) + zerOff;
% plot(f.log_m, diff);
% yscale_symlog


function v = log_correct(x)
    %v = x;
    v = -log(abs(x));
end