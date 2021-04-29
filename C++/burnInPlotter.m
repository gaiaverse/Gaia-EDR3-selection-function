colormap(winter);
system('scp "tplxdt071.nat.physics.ox.ac.uk:Documents/Work/GaiaSelectionFunction/Code/C++/burnTest*" ../C++/');
files = "burnTesting_3.txt";

dotSize = 20;

global x y z t success nRuns cx cy cz burn tx cutoff;
x = [];
y = [];
z = [];
t = [];
tx = [];
cx = [];
cy =[];
burn = 0;
cz = [];
success = false(0);
nRuns = 0;
cutoff = 1000;

clf;
for file = files
    file
    fid = fopen(file);
    tline = fgetl(fid);
    while ischar(tline)


        r = string(tline);
        checkPos(r);
        checkDuration(r);
        checkConverge(r);
        tline = fgetl(fid);
    end
end

if length(success) < length(x)
	success(end+1) = 0;
end
T = tiledlayout(2,2);
% title(T,"Base Model","FontSize",25);

% nexttile(T);
% 
% colors = [1,0,0;0,0,0];
% cs = colors(success+1,:);
% scatter3(x,y,z,size,cs,'filled');
% title("Starting positions");
% xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");



nexttile(T);
edges = [0:10:120];
t = t';

r = [];
b = [];
m = [];
for i = 1:width(t)
    p = t(:,i);
    b = [b;ones(height(t),1)*i];
    r =[r;p];
    m(i) = mean(p(~isnan(p) & p > 0));
 end
histogram2(r,b,40)
 view(-12,45);


title("Calculation duration (mean = " + mean(m) + ")");
xlabel("Seconds");
ylabel("Burn in");
zlabel("Frequency");

nexttile(T);

plot([1:length(m)],m);
hold on;
plot([1:length(old_m)],old_m);
hold off;
xlabel("Burn in width");
ylabel("Mean duration");
legend(["Iterative-n Burn ins", "Block-n Burn ins"]);

nexttile(T);
surf(t);
xlabel("burn in value");
ylabel("random seed - 2770");
zlabel("Duration");
% 
% nexttile(T);
% 
% 
% scatter3(x(success),y(success),z(success),size,t(success));
% xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");
% title("Duration by starting position");
% colorbar;



%subplot(1,2,2)
nexttile(T);
s = (~isnan(t(:)'));
s = s(1:length(cx));
q = t(success(:))';
scatter3(cx(s),cy(s),cz(s),dotSize,tx(s));
xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");
title("Duration By Final Position");
view(45,45);
colorbar;


function checkPos(r)
global x y z success nRuns t burn;
	if contains(r, "BEGINNING")
        burn = burn + 1;
        nRuns = 0;
    end
    if contains(r,"Position Vector initialised");
		nRuns = nRuns + 1;
		[j,k] = regexp(r,"-?\d*\.{0,1}\d+");
		
		pos = [0,0,0];
		for i = 1:length(j)
			w = extractBetween(r,j(i),k(i));
			pos(i) = str2double(w);
		end
		x(burn,nRuns) = pos(1);y(burn,nRuns)=pos(2);z(burn,nRuns) = pos(3);
		
%         x(burn,nRuns) = 0; y(burn,nRuns) = 0; z(burn,nRuns)=0;

		success(burn,nRuns) = false;
		t(burn,nRuns) = 0;
        
	end
end

function checkDuration(r)
	global success t nRuns burn tx cutoff;

	if contains(r,"Duration")
		[j,k] = regexp(r,"\d+");

        mins = 0;
        secs= 0;
        selector = 1;
		if contains(lower(r),"minute")
			mins = str2double(extractBetween(r,j(selector),k(selector)));
            selector = 2;
        end
		if contains(lower(r),"second")
            secs =  str2double(extractBetween(r,j(selector),k(selector)));
        end
        t(burn, nRuns) = mins*60 + secs;
        
       
		if t(burn,nRuns) > cutoff
			t(burn,nRuns) = NaN;
        end
        tx(end+1) = t(burn,nRuns);
		success(burn,nRuns) = true;
	end	
end

function checkConverge(r)
	global cx cy cz;
	if contains(r, "Converged roots")
	
		[j,k] = regexp(r,"-?\d*\.{0,1}\d+");
		
		pos = [0,0,0];
		for i = 1:3
			w = extractBetween(r,j(i),k(i));
			pos(i) = str2double(w);
		end
		cx(end+1) = pos(1);cy(end+1)=pos(2);cz(end+1) = pos(3);


	end

end