function [outputArg1,outputArg2] = catalogue(fileLocation,fileName)

	properties = readtable(fileLocation + "/OptimiserProperties.dat","ReadRowNames",true,"Delimiter","=");
	pData = table2array(properties)';
    vnames = properties.Properties.RowNames;
    properties = array2table(pData,"VariableNames",vnames);

	
	Nt = properties.Nt(1);
	Nl = properties.Nl(1);
	Nm = properties.Nm(1);
	
	t = 1717.6256+(linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4;
	tmin = t(1);
	tmax = t(2);
	
	z= readmatrix(fileLocation + fileName);
	
	xt = z(1:Nt);
	t = linspace(tmin,tmax,Nt);
	p = 1./(1 + exp(-xt));
	

	begins = [];
	ends = [];
	durations = [];
	quality = [];
	enclosed = [];
	gapThreshold = 0.9;
	
	nonSpuriousThreshold = 0.05;
	spuriousCost = 0.25;
	spuriousAlpha = log(spuriousCost)/(nonSpuriousThreshold - gapThreshold);
	spuriousB = exp(-spuriousAlpha * gapThreshold);
	
	inAGap = false;
	gapID = 0;
	
	cutNt = Nt;
	
	for i = 1:cutNt
		
		if p(i) < gapThreshold
			if ~inAGap
				inAGap = true;
				gapID = gapID + 1;
				begins(gapID) = round(t(i),3);
				quality(gapID) = 0;
				enclosed(gapID) = i;
			end
			
			if p(i) > nonSpuriousThreshold
				quality(gapID) = quality(gapID) + spuriousB * exp(spuriousAlpha * p(i));
			end
			
		else
			if inAGap
				inAGap = false;
				ends(gapID) = round(t(i),3);
				durations(gapID) = round(ends(gapID) - begins(gapID),3);
				enclosed(gapID) = i - enclosed(gapID);
			end
		end
			
	end
	
	if inAGap
		ends(gapID) = tmax;
		durations(gapID) = ends(gapID) - begins(gapID);
		enclosed(gapID) = cutNt - enclosed(gapID);
	end

	quality = round(100*(1-quality ./ enclosed),1);
	
	M = [begins', ends', durations', quality'];
	
	thresh = 80;
	goodGaps = array2table(M(quality > thresh,:));
	spuriousGaps = array2table(M(quality <= thresh,:));
	
	vNames = ["Begin","End","Duration","Quality"];
	goodGaps.Properties.VariableNames = vNames;
	spuriousGaps.Properties.VariableNames = vNames;
	
	writetable(goodGaps,fileLocation + "/confidentGapCatalogue.csv");
	writetable(spuriousGaps,fileLocation + "/suspectGapCatalogue.csv");

	
	cutter = t < t(cutNt);
	plot(t(cutter),p(cutter));
end

