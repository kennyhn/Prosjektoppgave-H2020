function [u, pixels] = readContourPlot(contourfile,xrange,yrange,scalerange,showFigure)
    
    Ncolors = 64; %Choose number of colors in colormap
    Ngray = 3;
        
    RGB = imread(contourfile); %
    ipixels = size(RGB,1);
    jpixels = size(RGB,2);
    pixels = [ipixels jpixels];
    cmap = [jet(Ncolors); gray(Ngray)];
    
    X = zeros(ipixels,jpixels);
   
    for i =1:ipixels
        for j = 1:jpixels
            RGBcell(1:3) = double(RGB(i,j,:))./255;
            error = sqrt(sum(((cmap-RGBcell)./(RGBcell)).^2,2));
            minError = min(error);
            indx = find(error==minError);
            if length(indx) > 1
                indx = indx(1);
            end
            X(i,j) = indx;
        end
    end

    uscaleIncr = (scalerange(2)-scalerange(1))/(Ncolors-1);
    uscale = [scalerange(1):uscaleIncr:scalerange(2) nan(1,Ngray)];
    u = uscale(X);

    a = jpixels/(max([jpixels ipixels]));
    b = ipixels/(max([jpixels ipixels]));
    if showFigure == true
        figure, contourf(flip(u),256,'edgecolor','none'), colorbar, colormap(cmap)
        pbaspect([a b 1])

        xtIncr = floor(jpixels/10);
        ytIncr = floor(ipixels/10);
        xt = 1:xtIncr:jpixels;
        yt = 1:ytIncr:ipixels;
        xticks(xt);
        yticks(yt);

        xincr = (xrange(2)-xrange(1))/(length(xt)-1);
        xvals = xrange(1):xincr:xrange(2);
        yincr = (yrange(2)-yrange(1))/(length(yt)-1);
        yvals = yrange(1):yincr:yrange(2);
        xlabels = {};
        ylabels = {};
        for i=1:length(xt)
           xlabels{1,i} = num2str(xvals(i));
        end
        for i=1:length(yt)
           ylabels{1,i} = num2str(yvals(i));
        end
        xticklabels(xlabels)
        yticklabels(ylabels)
    end
end
