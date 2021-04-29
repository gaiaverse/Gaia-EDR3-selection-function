% system('scp "tplxdt071.nat.physics.ox.ac.uk:Documents/Work/GaiaSelectionFunction/Code/C++/burnSearcher*" ../C++/');

% files = "burnSearcher_long" + ["1","2","3","4"]+ ".txt";
files = "smallInit2.txt";
dotSize = 20;

global x y z t success nRuns cx cy cz;
x = [];
y = [];
z = [];
t = [];
cx = [];
cy =[];
cz = [];
success = false(0);
nRuns = 0;


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
T = tiledlayout(1,2);
title(T,"Base Model","FontSize",25);

% nexttile(T);
% 
% colors = [1,0,0;0,0,0];
% cs = colors(success+1,:);
% scatter3(x,y,z,size,cs,'filled');
% title("Starting positions");
% xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");




nexttile(T);
edges = [0:10:120];
histogram(t,20);
title("Calculation duration (mean = " + num2str(mean(t(~isnan(t)))) + ")");
xlabel("Seconds");
ylabel("Counts");

% 
% nexttile(T);
% 
% 
% scatter3(x(success),y(success),z(success),size,t(success));
% xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");
% title("Duration by starting position");
% colorbar;




nexttile(T);
s = (cz >0.1);
q = t(success);
scatter3(cx(s),cy(s),cz(s),dotSize,q(s));
xlabel("Nt_1"); ylabel("Nt_2"); zlabel("Ns");
title("Duration By Final Position");
view(45,45);
colorbar;


function checkPos(r)
global x y z success nRuns t;
	if contains(r, "BEGINNING")
		nRuns = nRuns + 1;
% 		[j,k] = regexp(r,"-?\d*\.{0,1}\d+");
% 		
% 		pos = [0,0,0];
% 		for i = 1:length(j)
% 			w = extractBetween(r,j(i),k(i));
% 			pos(i) = str2double(w);
% 		end
% 		x(end+1) = pos(1);y(end+1)=pos(2);z(end+1) = pos(3);
		
        x(nRuns) = 0; y(nRuns) = 0; z(nRuns)=0;

		success(nRuns) = false;
		t(nRuns) = NaN;
	end
end

function checkDuration(r)
	global success t nRuns;

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
        t(nRuns) = mins*60 + secs;
% 		if t(nRuns) > 60
% 			t(nRuns) = NaN;
% 		end
		success(nRuns) = true;
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