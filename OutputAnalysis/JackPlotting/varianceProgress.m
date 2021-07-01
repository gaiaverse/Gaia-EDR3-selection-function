function varianceProgress(folders,startN,stopN, gap, includeFinal)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex');

    textures = ["-","--",":","-."];
    figure(3);
        clf;
    for z = 1:length(folders)
        folder = folders(z);
        propertyFile = "../../../CodeOutput/" + folder + "/OptimiserProperties.dat";
        p = readtable(propertyFile,"ReadRowNames",true,"Delimiter","=");
        pData = table2array(p)';
        vnames = p.Properties.RowNames;
        properties = array2table(pData,"VariableNames",vnames);

        hyperOrder = properties.hyperOrder(1);
        nVariancePopulations= properties.NVariancePopulations(1);

        totalN = properties.totalTransformedParams(1);
        start = totalN - (hyperOrder+2)*nVariancePopulations+1;


        v = [startN:gap:stopN];
        files = "../../../CodeOutput/" + folder + "/TempPositions/TempPosition" + string(v) +"_TransformedParameters.dat";

        if includeFinal == true
            files = [files, "../../../CodeOutput/" + folder + "/FinalPosition_TransformedParameters.dat"];
            v = [v, stopN+gap];
        end

        terms= [];
        fracs = [];
        varRange = 1:nVariancePopulations;
        hypRange = 1:round(hyperOrder) + 1;
        for i = 1:length(files)
            files(i)
            f = readmatrix(files(i));
            variance = f(start:end);
             for j= varRange
                for k = hypRange  
                    pFac = k-1;
                    if k == 1
                        pFac = 1;
                    end
                   terms(i,k,j) = variance((k-1)*nVariancePopulations + j)^pFac; 
                end
                fracs(i,j) = variance((hyperOrder+1)*nVariancePopulations + j);
            end

        end

        
        nx = floor(sqrt(hyperOrder +4));
        ny = ceil((hyperOrder + 3)/nx);

        nF = mod(nx*ny,hyperOrder+2);
        s = [1:nF] + hyperOrder + 2;
        nBar = 100;
        m = colororder;
        for j = 1:nVariancePopulations

            sum = zeros(size(v));
            for k = 1:hyperOrder + 1
               subplot(nx,ny,k)
               hold on;
               plot(v,terms(:,k,j),textures(z),'Color',m(j,:));

               hold off;
               sum = sum + terms(:,k,j)*nBar^(k-1);
            end

            subplot(nx,ny,hyperOrder+2);
            hold on;
            plot(v,sum,textures(z),'Color',m(j,:));
            hold off;


            subplot(nx,ny,s);
            hold on;
            plot(v,fracs(:,j),textures(z),'Color',m(j,:));
            hold off;
        end

        for k = 1:hyperOrder+3
            sp= k;
            if k == hyperOrder + 3
              sp = s;  
            end
            subplot(nx,ny,sp)
            set(gca,'yscale','log')
%             set(gca,'xscale','log')
%             xlim([max(gap,startN)-0.5,stopN])
            if k < hyperOrder+2
               title("$n^"+string(k-1)+"$ term");
            else
                if k < hyperOrder + 3
                   title("Induced Variance for star with n = " + num2str(nBar)); 
                else
                    title("Population Fraction");
%                     legend(folders)
                end
            end
           xlabel("Complete Epochs");
           ylabel("Optimiser Value");
        end
    end
end

