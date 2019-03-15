clear all
close all
addpath('./Manual/');
addpath('../config/');

% load config file. Config contains parameters for some methods
load config;
%define config as global variable
global config;

if config.calc_size == false
    width = config.cropRect(3);
    actualPixelSize = ones(1,width);

    save([config.frontdir,'/H.mat'],'actualPixelSize');
else
   
    for i=1:length(config.satdatefolder)
        frontfolder_i = [config.satdir,config.satdatefolder{i}];
        frontmat_i = [frontfolder_i,'/front.mat'];
        if exist(frontmat_i)
            load(frontmat_i);
            continue
        end

        Date_it = datetime(config.satdatefolder{i},'InputFormat','yyyy-MM-dd');
        if Date_it < config.imageFileDate{1}
            startImagefiledir = config.imagefiledir{1};
        else
            test = find(cellfun(@(d) datenum(d),config.imageFileDate)-datenum(Date_it)<1);
            startImagefiledir = config.imagefiledir{test(end)};
        end
        currentImage0 = imread(startImagefiledir);
        currentImage0 = rgb2gray(currentImage0);
        currentImage = imcrop(currentImage0,config.cropRect);
        [height,width] = size(currentImage);

        % Calculate distance of front on image
        f = figure
        imshow(currentImage)
        hold on
        if i>1
            [xim_1,yim_1] = interp_pts(ximp,yimp,15);
            linetoedit = line(xim_1,yim_1);
        else
            linetoedit = line([10,width],[height,10]);
        end
        [xim,yim] = editpolyline(linetoedit);
    %    d = dist(cat(1,xim,yim),'Euclidean');
        d = distance(xim,yim);
        size_tl = 0;
        for j=1:length(d)-1
            size_tl = size_tl + d(j,j+1);
        end

        [ximp,yimp] = interp_pts(xim,yim,floor(size_tl));

        satFile = dir(fullfile([config.satdir,config.satdatefolder{i}],strcat('*B8.','TIF')));
        filename = satFile.name;
        geoinfo = geotiffinfo(fullfile(config.satdir,config.satdatefolder{i},filename));
        if geoinfo.Zone ~= 33
            satFile = dir(fullfile([config.satdir,config.satdatefolder{i}],strcat('*B8_33.','TIF')));
            filename = satFile.name;
            geoinfo = geotiffinfo(fullfile(config.satdir,config.satdatefolder{i},filename));
        end
        refbox = geoinfo.BoundingBox;
        ximref1 = refbox(1,1); % Lower left corner
        yimref1 = refbox(1,2); % Lower left corner
        ximref2 = refbox(2,1); % Upper right corner
        yimref2 = refbox(2,2); % Upper right corner

        px_x = geoinfo.PixelScale(1);
        px_y = geoinfo.PixelScale(2);
        [sat_im,map] = imread(fullfile(config.satdir,config.satdatefolder{i},filename));

        % Zoom
        xz1 = (config.ximz1-ximref1)/px_x+1;
        xz2 = (config.ximz2-ximref1)/px_x+1;
        yz1 = -(config.yimz1-yimref2)/px_y+1;
        yz2 = -(config.yimz2-yimref2)/px_y+1;

        %%
        % Reference points on each side of the front for start and end of line
        if i == 1
            x1 = (config.xl1-ximref1)/px_x+1;
            x2 = (config.xl2-ximref1)/px_x+1;
            y1 = -(config.yl1-yimref2)/px_y+1;
            y2 = -(config.yl2-yimref2)/px_y+1;
        else
            x1 = (xiref-ximref1)/px_x+1;
            y1 = -(yiref-yimref2)/px_y+1;
            x1 = x1(1:200:end);
            y1 = y1(1:200:end);
        end

        % Plot the point
        imsize = size(sat_im);

        % Camera specs
        focal = config.focalLength*1e3;
        angl = 2*atan(config.xsens/(2*focal));
        yangl = 2*atan(config.ysens/(2*focal));

        % Limits of camera view
        xlim1 = (config.xlimref1-ximref1)/px_x+1;
        ylim1 = -(config.ylimref1-yimref2)/px_y+1;

        xcamim = (config.xcam-ximref1)/px_x+1;
        ycamim = -(config.ycam-yimref2)/px_y+1;

        fun= @(xplot) ((ycamim-ylim1)/(xcamim-xlim1))*(xplot-xcamim)+ycamim;
        xplot = xcamim:imsize(2);
        yplot = fun(xplot);
        rot_mat = @(angl) [cos(angl) -sin(angl); sin(angl) cos(angl)];
        mat = rot_mat(-angl);
        mat2 = rot_mat(-angl/2);
        f_rot = mat*[xplot-xcamim;((ycamim-ylim1)/(xcamim-xlim1))*(xplot-xcamim)+ycamim-ycamim];
        f_rot2 = mat2*[xplot-xcamim;((ycamim-ylim1)/(xcamim-xlim1))*(xplot-xcamim)+ycamim-ycamim];
        f = figure;
        imshow(sat_im,map);hold on
        axis([xz1,xz2,yz2,yz1]);
        plot(xplot,yplot,f_rot(1,:)+xcamim+1,f_rot(2,:)+ycamim+1,'color','r','linewidth',2,'linestyle','--');axis equal
        axis([xz1,xz2,yz2,yz1]);

        % Draw the line
        if i==1
            linetoedit = line([x1 x2],[y1,y2]);
        else
            linetoedit = line(x1,y1);
        end
        [x,y] = editpolyline(linetoedit);

        % Rotate and translate axes
        anglcam = -atan((yplot(2)-yplot(1))/(xplot(2)-xplot(1)))+angl/2;
        rot_matcam = @(angl) [cos(angl) sin(angl); -sin(angl) cos(angl)];
        matcam = rot_matcam(anglcam+pi/2);
        f_rotcam = matcam*[x;y];
        xr = f_rotcam(1,:);
        yr = f_rotcam(2,:);
        [xr,ir] = sort(xr);
        yr = yr(ir);
        pts = floor(size_tl); % Number of points to interpolate
        xir = linspace(min(xr),max(xr),pts);
        yir = interp1(xr,yr,xir);
        matcami = rot_matcam(-pi/2-anglcam);
        fi_rotcam = matcami*[xir;yir];
        xi = fi_rotcam(1,:);
        yi = fi_rotcam(2,:);


        plot(x,y,'marker','x'); axis equal
        hold on
    %%
        line(xi,yi,'marker','o','color','g');
        axis([xz1,xz2,yz2,yz1]);

        % Real coordinates
        xiref = ximref1+px_x*(xi-1);
        yiref = (yimref2-px_y*(yi-1));
        A = cat(1,xiref,yiref)';
        distref = sqrt( sum( abs( diff( A ) ).^2, 2 ) );
        xdistref = movsum(ximp,2)/2;
        xpxsize = interp1(xdistref(2:end-2),distref,1:width);
        % Calculate distance between camera and front for each point
        Hp = sqrt((sqrt((xiref-config.xcam).^2+(yiref-config.ycam).^2)).^2+config.zcam^2);

        % Interpolate for width
        H = interp1(ximp(1:end-2),Hp,1:width);
        ypxsize = (config.imagePixelSize / config.focalLength) * H;

        actualPixelSize = xpxsize.*ypxsize;

        % Take only the values in the cropped area
        %H = H(cropRect(1):cropRect(1)+cropRect(3));

        frontdir_i = [config.frontdir,config.satdatefolder{i}];
        if ~exist(frontdir_i)
            mkdir(frontdir_i);
        end 

        save([frontdir_i,'/H.mat'],'H','xpxsize','ypxsize','actualPixelSize');
        save([frontdir_i,'/front.mat'],'ximp','yimp','xiref','yiref','size_tl');
        dlmwrite([frontdir_i,'/front.txt'],cat(1,xiref,yiref)','precision','%.2f')
    end
end
