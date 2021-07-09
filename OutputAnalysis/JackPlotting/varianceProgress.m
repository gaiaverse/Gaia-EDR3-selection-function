function varianceProgress(folders,minLim)
    if nargin < 2
        minLim = 0;
    end
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


        v = [0:1:nVariancePopulations-1];
        files = "../../../CodeOutput/" + folder + "/TempPositions/Population" + string(v) +"HyperParams.dat";     

        terms= [];
        fracs = [];
        varRange = 1:nVariancePopulations;
        hypRange = 1:round(hyperOrder) + 1;
        for i = 1:length(files)
            files(i)
            f = readmatrix(files(i));
            
            fracs(:,i) = f(:,2);
            for j = hypRange
               rawterms(:,j,i) = f(:,2+j); 
            end
        end
        
        L = [1:height(f)];
        maxxer = L(end);
        if minLim > maxxer
            minLim = 0;
        end
        bounds = [minLim, maxxer];
        terms = zeros(size(rawterms));
        for i = varRange
            
            terms(:,1,i) = rawterms(:,1,i).^2;

            for j = 1:hyperOrder+1
                v = zeros(size(L))';
                q = j - 1;
                for k = max(1,ceil(q/2)):hyperOrder/2
                   v = v + nchoosek(2*k,2*k-q) * rawterms(:,2*k,i).^(2*k-q) .* rawterms(:,2*k+1,i).^q;
                end

                terms(:,j,i) = terms(:,j,i) + v; 
 
            end
        end
        
%         dim = 2;
%         [terms(:,dim,1),terms(:,dim,2),terms(:,dim,3)]

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
               plot(L,terms(:,k,j),textures(z),'Color',m(j,:));

               hold off;
               sum = sum + terms(:,k,j)*nBar^(k-1);
            end

            subplot(nx,ny,hyperOrder+2);
            hold on;
            plot(L,sum,textures(z),'Color',m(j,:));
            hold off;


            subplot(nx,ny,s);
            hold on;
            plot(L,fracs(:,j),textures(z),'Color',m(j,:));
            hold off;
        end

        for k = 1:hyperOrder+3
            sp= k;
            if k == hyperOrder + 3
              sp = s;  
            end
            subplot(nx,ny,sp)
%             set(gca,'yscale','log')
%             set(gca,'xscale','log')
%             xlim([max(gap,startN)-0.5,stopN])
            if k < hyperOrder+2
               title("$n^"+string(k-1)+"$ term");
            else
                set(gca,'yscale','log')
                if k < hyperOrder + 3
                   title("Induced Variance for star with n = " + num2str(nBar)); 
                   
                else
                    title("Population Fraction");
%                     legend(folders)
                end
            end
            xlim(bounds);
           xlabel("Minibatch Evaluations");
           ylabel("Optimiser Value");
        end
    end
end

