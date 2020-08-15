%% important info
% plot only hoop at left and right frames
% m06: left 415 +- 15; right 815 +- 15
% m01: left 375 +- 15; right 855 +- 10
% if errors happen with prev_centers, change the threshold in 'track the
% body' part to 0.5 / 0.6
x_left = 415;
l_oscillate = 15;
x_right = 815;
r_oscillate = 15;

%% variables
% user inputs
% filedir = 'C:\Users\Ristroph\Dropbox\hula hoop\motion tracking\imgs\cone_down\sequence_img\Photos\m06img';
filedir = '/Users/arielzhu/Desktop/hula_hoop_project/hula_hoop_motion_tracking/ariel/cone_m06_normal_test_more/m06img';
% filedir = 'images/hyper_sequence_10/m04img';
begini = 31; 
endi = 200;
% store useful info
all_centers_top = [];
all_majors = [];
all_minors = [];
all_angles = [];
all_centers_body = [];
all_radius_body = [];
all_heights = [];
all_x = [];
all_majoraxis_top = [];
all_minoraxis_top = [];
all_orientation_top = [];
all_centroids = []; % side info

%% arrow 
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );

%% track the hoop
% loop over all images
for i = begini:endi
    % build the file name string
    if i < 10
        file = strcat(filedir,'000',int2str(i),'.bmp');
    elseif i < 100
        file = strcat(filedir,'00',int2str(i),'.bmp');
    elseif i < 1000
        file = strcat(filedir,'0',int2str(i),'.bmp');
    else file = strcat(filedir,int2str(i),'.bmp');
    end
    
    % read in the file, threshold and cut into top and side views
    img = imread(file);
    threshold = 0.3; % originally 0.3
    img = im2bw(img,threshold); % change a image to bianry image to black and white
    imgtop = img(1:1100,:);
    imgside =img(1101:1800,:);
    
    % get the user to build masks of the body in the first image
    if i == begini
        maskside = roipoly(imgside);
        bodyside = immultiply(maskside,imgside);
        bodysidex = sum(bodyside,1);
        bodysidey = sum(bodyside,2);
        
        masktop = roipoly(imgtop);
        bodytop = immultiply(masktop,imgtop);
        bodytopx = sum(bodytop,1);
        bodytopy = sum(bodytop,2);
    end
    
    % calculate body location in current image, both views
    % uses finddelay to get shifts in x and y
    imgsidex = sum(imgside,1);
    imgsidey = sum(imgside,2);
    maxshift = 200;
    xside = finddelay(bodysidex,imgsidex,maxshift);
    yside = finddelay(bodysidey,imgsidey,maxshift);
    
    imgtopx = sum(imgtop,1);
    imgtopy = sum(imgtop,2);
    xtop = finddelay(bodytopx,imgtopx,maxshift);
    ytop = finddelay(bodytopy,imgtopy,maxshift);
    
    % generate new shifted body masks in both views
    newmaskside = imtranslate(maskside,[xside,yside]);
    newmasktop = imtranslate(masktop,[xtop,ytop]);
    
    % subtract out new body mask to get hoop
    hoopside = im2bw(imsubtract(imgside,newmaskside)); 
    hooptop = im2bw(imsubtract(imgtop,newmasktop));

    % filter out blobs based on area and proximity to middle
    minarea = 200; % blobs must be larger than this
%     mindist = 50; % blobs must be further than this from any edge
    mindist = 15;
   
    cc = bwconncomp(hoopside); 
    stats = regionprops(cc,'Area','Centroid'); 
    areas = [stats.Area];
    centroids = cat(1,stats.Centroid);
    idx = find( areas(:) > minarea & min(centroids(:,1),size(hoopside,2)-centroids(:,1)) > mindist & min(centroids(:,2),size(hoopside,1)-centroids(:,2)) > mindist ); 
    hoopside = ismember(labelmatrix(cc),idx); 
    
    cc = bwconncomp(hooptop); 
    stats = regionprops(cc,'Area','Centroid'); 
    areas = [stats.Area];
    centroids = cat(1,stats.Centroid);
    idx = find( areas(:) > minarea & min(centroids(:,1),size(hooptop,2)-centroids(:,1)) > mindist & min(centroids(:,2),size(hooptop,1)-centroids(:,2)) > mindist ); 
    hooptop = ismember(labelmatrix(cc),idx); 
    
    % build convex hull: combine several parts of the hoop together as a
    % complete hoop, and fill up the inside of the hoop
    hoophullside = bwconvhull(hoopside);
    hoophulltop = bwconvhull(hooptop);
    
    % get centroid and eccentricity/angle info
    statside = regionprops(hoophullside,{'Centroid','MajorAxisLength','MinorAxisLength','Orientation'}); 
    stattop = regionprops(hoophulltop,{'Centroid','MajorAxisLength','MinorAxisLength','Orientation'}); 
    % store info
    all_centers_top(i-begini+1,:) = stattop.Centroid;
    all_majors(i-begini+1) = statside.MajorAxisLength;
    all_minors(i-begini+1) = statside.MinorAxisLength;
    all_angles(i-begini+1) = statside.Orientation;
    all_heights(i-begini+1) = statside.Centroid(2);
    all_x(i-begini+1,:) = statside.Centroid(1);
    all_centroids(i-begini+1,:) = statside.Centroid;
    all_majoraxis_top(i-begini+1) = stattop.MajorAxisLength;
    all_minoraxis_top(i-begini+1) = stattop.MinorAxisLength;
    all_orientation_top(i-begini+1) = stattop.Orientation;
end

%% track the body
k = 1;
prev_centers = [];
for i = begini:endi
    % build the file name string
    if i < 10
        file = strcat(filedir,'000',int2str(i),'.bmp');
    elseif i < 100
        file = strcat(filedir,'00',int2str(i),'.bmp');
    elseif i < 1000
        file = strcat(filedir,'0',int2str(i),'.bmp');
    else file = strcat(filedir,int2str(i),'.bmp');
    end
    img = imread(file);
    threshold = 0.6; % if start from 31, change it back to 0.6
    img = im2bw(img,threshold); % change a image to bianry image to black and white
    img = imcomplement(img); % reverse black and white
    imgtop = img(1:1100,:);

% attempt: directly find the circle with suitable size
%     imshow(imgtop);
    [centers, radius] = imfindcircles(imgtop, [10 30]);
    % check the size of centers, i.e., check how many circles are found
    size_centers = size(centers);
    if size_centers(:,1) == 1
        all_centers_body(k,:) = centers;
        prev_centers = centers;
        all_radius_body(k,:) = radius;
    else
        % deal with the case when several circles are detected
        dist = 10000000;
        best_centers = [];
        for c = 1:size_centers(:,1)
            curr_centers = centers(c,:);
            curr_radius = radius(c,:);
            % SPECIAL THINGS NEEDED HERE IF FIRST DETECTION HAS MULTIPLE
            % CIRCLES
            if(norm(curr_centers-prev_centers) < dist)
                all_centers_body(k,:) = curr_centers;
                all_radius_body(k,:) = curr_radius;
                dist = norm(curr_centers-prev_centers);
                best_centers = curr_centers;
            end
        end
        prev_centers = best_centers;
    end
%     viscircles(centers, radius, 'EdgeColor', 'b');
%     pause;
    k = k + 1;
end

%% shift the body
all_body_centers_shift = [];
horizontal_shift = -5; 
vertical_shift = 22; 
for k = 1:endi-begini+1
    all_body_centers_shift(k,:) = [all_centers_body(k,1)+horizontal_shift, all_centers_body(k,2)+vertical_shift];
end

%% processing (plot)
%% Normal Mode
figure('Name', 'Normal Mode Procedure');
% visualize the unshifted gyration
% figure('Name', 'Gyration');
subplot(2,2,1), hold on
plot(all_centers_body(:,1), all_centers_body(:,2), '.'), hold on
title('Gyration & Fitted Circle');
xlabel('x-position');
ylabel('y-position');
axis equal

% fit a circle
[xg, yg, Rg, a] = circfit(all_centers_body(:,1), all_centers_body(:,2));
viscircles([xg, yg], Rg, 'Edgecolor', 'r'), hold on
plot(xg, yg, '+r'), hold off;

% fit the shifted circle
[xg_s, yg_s, Rg_s, a_s] = circfit(all_body_centers_shift(:,1), all_body_centers_shift(:,2));

% % get xb, yb against xg, yg
% body_centers_against_gyration = [];
% tester_angle = [];
% all_theta = [];
% all_star = [];
% for k = 1:(endi - begini + 1)
%     xb_g = all_centers_body(k,1) - xg;
%     yb_g = yg - all_centers_body(k,2);
%     % get theta
%     theta = atan2(yb_g, xb_g);%%%%%CHECK%%%%%     atan2(y,x);
%     all_theta(k) = theta;
%     body_centers_gyration(k,:) = [xb_g yb_g];
%     % get the angle between the "hoop vector" and the "body vector"
%     star = atan2(all_centers_body(k,2)-all_centers_top(k,2), all_centers_top(k,1)-all_centers_body(k,1));
%     all_star(k) = star;
%     % check for some mis-case of atan2
%     if abs(star - theta) > 2
%         tester_angle(k) = star + theta;
%     else
%         tester_angle(k) = star - theta;
%     end
% end

% shifted version
% get xb, yb against xg, yg
body_centers_against_gyration_s = [];
tester_angle_s = [];
all_theta_s = [];
all_star_s = [];
for k = 1:(endi - begini + 1)
    xb_g_s = all_body_centers_shift(k,1) - xg_s;
    yb_g_s = yg_s - all_body_centers_shift(k,2);
    % get theta
    theta_s = atan2(yb_g_s, xb_g_s);%%%%%CHECK%%%%%     atan2(y,x);
    all_theta_s(k) = theta_s;
    body_centers_against_gyration_s(k,:) = [xb_g_s yb_g_s];
    % get the angle between the "hoop vector" and the "body vector"
    star_s = atan2(all_body_centers_shift(k,2)-all_centers_top(k,2), all_centers_top(k,1)-all_body_centers_shift(k,1));
    all_star_s(k) = star_s;
    % check for some mis-case of atan2
    if abs(star_s - theta_s) > 2
        tester_angle_s(k) = theta_s - star_s;
    else
        tester_angle_s(k) = star_s - theta_s;
    end
end

% plot theta
subplot(2,2,2), hold on
% figure('Name', 'Theta');
plot([1:endi-begini+1], all_theta_s, '.');
title('Theta');
xlabel('Frame');
ylabel('Theta');

% plot star
subplot(2,2,3), hold on
% figure('Name', 'Star');
plot([1:endi-begini+1], all_star_s, '.');
title('Star');
xlabel('Frame');
ylabel('Star');

% final plot
subplot(2,2,4), hold on
% figure('Name', 'Normal Mode');
plot([1:endi-begini+1], tester_angle_s, '.');
title('Normal Mode Tester Angle: star - theta');
xlabel('Frame');
ylabel('tester angle in radians');
ylim([-2 2]);

%% z-plot
figure('Name', 'z-plot Procedure');
subplot(2,2,1), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate || abs(all_x(i-begini+1)-x_right) <= r_oscillate
        % see if to be ignored
        if i > begini
            if abs(all_heights(i-begini+1)-all_heights(i-begini)) <= 30
                plot(i, all_heights(i-begini+1), 'b.'), hold on
            end
        end
    end
end
title('z-plot at left and right most frames (equal axis)');
xlabel('frame');
ylabel('z (in Matlab frame)');
axis equal

subplot(2,2,2), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, all_heights(i-begini+1), 'b.', 'MarkerSize', 10), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        plot(i, all_heights(i-begini+1), 'r.', 'MarkerSize', 10), hold on
    end
end
title('z-plot at left and right most frames');
xlabel('frame');
ylabel('z (in Matlab frame)');

% for all
subplot(2,2,3), hold on
plot([begini:endi], all_heights, '.'), hold on
title('z-plot for all frames');
xlabel('frame');
ylabel('z (in Matlab frame)');

% x
% for all, and label the right and left frames
subplot(2,2,4), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate || abs(all_x(i-begini+1)-x_right) <= r_oscillate
        plot(i, all_x(i-begini+1), 'r.'), hold on
    else
        plot(i, all_x(i-begini+1), 'b.'), hold on
    end
end
title('x-plot ');
xlabel('frame');
ylabel('x (in Matlab frame)');

%% sagging angle
figure('Name', 'Sagging Angle Procedure');
subplot(1,2,1), hold on
% only plot those on the left and right frames
% some data comes from observation in 'plot_z.m'
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, -all_angles(i-begini+1), 'b.'), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        % for right frames, plot the positive value of the angle
        plot(i, all_angles(i-begini+1), 'r.'), hold on
    end
end
title('Sagging Angles at Right & Left with Default Axis');
xlabel('Frame');
ylabel('Sagging angle (degrees) (blue left, red right)');

subplot(1,2,2), hold on
plot([begini:endi], all_angles, 'b.');
title('All Sagging Angles');
xlabel('Frame');
ylabel('Sagging angle (degrees)');

%% final plots
figure('Name', 'Motion Tracking Final Plots');
title('Motion Tracking Final Plots');
subplot(2,2,1), hold on
plot([1:endi-begini+1], tester_angle_s, '.');
title('Normal Mode Tester Angle: T = \gamma - \theta');
xlabel('Frame');
ylabel('tester angle T (radians)');
ylim([-2 2]);
subplot(2,2,2), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, all_heights(i-begini+1), 'b.'), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        plot(i, all_heights(i-begini+1), 'r.'), hold on
    end
end
title('z-plot at left and right most frames');
xlabel('frame');
ylabel('z (in Matlab frame)');
subplot(2,2,3), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, -all_angles(i-begini+1), 'b.'), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        % for right frames, plot the positive value of the angle
        plot(i, all_angles(i-begini+1), 'r.'), hold on
    end
end
title('Sagging Angles at Right & Left with Default Axis');
xlabel('Frame');
ylabel('Sagging angle (degrees) (blue left, red right)');

%% general plot
figure('Name', 'General Plot');
title('Motion Tracking General Plot');
subplot(3,1,1), hold on
tester_angle_degree = [];
for i = 1:endi-begini+1
    tester_angle_degree(i) = tester_angle_s(i) * 57.2958;
end
plot([1:endi-begini+1], tester_angle_degree, '.');
% title('Normal Mode Tester Angle: \delta = \gamma - \theta');
xlabel('Frame');
ylabel({'tester angle';'\delta'});
set(gca, 'ylim', [-80 80]);
set(gca,'position',[0.15 0.65 0.75 0.27]);
set(gca,'box','on');
subplot(3,1,2), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, all_heights(i-begini+1), 'b.'), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        plot(i, all_heights(i-begini+1), 'r.'), hold on
    end
end
% title('z-plot at left and right most frames');
xlabel('frame');
ylabel({'z';'(Matlab Frame'});
set(gca,'position',[0.15 0.38 0.75 0.27]);
set(gca,'box','on');
subplot(3,1,3), hold on
for i = begini:endi
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        plot(i, -all_angles(i-begini+1), 'b.'), hold on
    elseif abs(all_x(i-begini+1)-x_right) <= r_oscillate
        % for right frames, plot the positive value of the angle
        plot(i, all_angles(i-begini+1), 'r.'), hold on
    end
end
% title('Sagging Angles at Right & Left with Default Axis');
xlabel('Frame');
ylabel('Sagging angle');
set(gca,'position',[0.15 0.11 0.75 0.27]);
set(gca,'box','on');

%% visualization
% visiualize the tracking of hoop and body
for i = begini:endi
    if i < 10
        file = strcat(filedir,'000',int2str(i),'.bmp');
    elseif i < 100
        file = strcat(filedir,'00',int2str(i),'.bmp');
    elseif i < 1000
        file = strcat(filedir,'0',int2str(i),'.bmp');
    else file = strcat(filedir,int2str(i),'.bmp');
    end
    img = imread(file);
    threshold = 0.3; % originally 0.3
    img = im2bw(img,threshold); % change a image to bianry image to black and white
    imgtop = img(1:1100,:);
    imgside =img(1101:1800,:);
    figure(1), hold on
    %% top
    subplot(2,1,2), imshow(imgtop), hold on
    % visualize the gyration i.e., the fixed circle
    viscircles([xg_s, yg_s], Rg_s, 'Edgecolor', 'c'), hold on
    t = linspace(0,2*pi,50);
    a = all_majoraxis_top(i-begini+1) / 2;
    b = all_minoraxis_top(i-begini+1) / 2;
    Xc = all_centers_top(i-begini+1,1);
    Yc = all_centers_top(i-begini+1,2);
    phi = deg2rad(-all_orientation_top(i-begini+1));
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    
    % visualize the outer hoop
    plot(x,y,'r','Linewidth',1), hold on
    
    % visualize the inner hoop
    inner = 17
    ao = all_majoraxis_top(i-begini+1) / 2 - inner;
    bo = all_minoraxis_top(i-begini+1) / 2 - inner;
    xo = Xc + ao*cos(t)*cos(phi) - bo*sin(t)*sin(phi);
    yo = Yc + ao*cos(t)*sin(phi) + bo*sin(t)*cos(phi);
    plot(xo,yo,'-.m','Linewidth',1.5), hold on
    
    body_centers = all_centers_body(i-begini+1, :);
    body_centers_shift = all_body_centers_shift(i-begini+1, :);
    body_radius = all_radius_body(i-begini+1);
    viscircles(body_centers, body_radius, 'EdgeColor', 'b'), hold on
    
    % visualize the vector connecting body center and hoop center
    plot(Xc, Yc, '+m'), hold on
    drawArrow([body_centers_shift(1),Xc], [body_centers_shift(2), Yc], 'color', 'm'), hold on
    
    % visualize the gyration center
%     plot(xg, yg, '+r'), hold on
   
    % vector connecting the gyration center and body center
    drawArrow([xg_s,body_centers_shift(1)], [yg_s, body_centers_shift(2)], 'color', 'r'), hold on
    
    % visualize the shifted body center and gyration center
    plot(xg_s, yg_s, '+m'), hold on
    plot(all_body_centers_shift(i-begini+1,1), all_body_centers_shift(i-begini+1,2), 'om'), hold on
    
    
    %% side
    subplot(2,1,1), imshow(imgside), hold on
    
    % plots the extracted ellipse
    major = all_majors(i-begini+1);
    minor = all_minors(i-begini+1);
    orientation = all_angles(i-begini+1);
    centroid = all_centroids(i-begini+1,:);
    t = linspace(0,2*pi,50);
    a = major/2;
    b = minor/2;
    Xc = centroid(1);
    Yc = centroid(2);
    phi = deg2rad(-orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'r','Linewidth',1)    
    % visualize the center
    plot(centroid(1), centroid(2), 'ro'), hold on
    %% sagging angle
    % for the selected plot frame, visualize the angle
    if abs(all_x(i-begini+1)-x_right) <= r_oscillate
%         edge_len = Xc * tand(orientation);
        line([Xc+a*cosd(orientation), Xc-a*cosd(orientation)], [Yc-a*sind(orientation), Yc+a*sind(orientation)], 'Color', 'm'), hold on
    end
    if abs(all_x(i-begini+1)-x_left) <= l_oscillate
        edge_len = Xc * tand(orientation);
        line([Xc-a*cosd(orientation), Xc+a*cosd(orientation)], [Yc+a*sind(orientation), Yc-a*sind(orientation)], 'Color', 'm'), hold on
    end
%     if fix(i-begini+1) == 1
%         % visualize the fixed centroid as well
%         fixed_centroid = all_new_centroids(i-begini+1,:);
%         plot(fixed_centroid(1), fixed_centroid(2), 'mo'), hold on
%     end
%     
%     if i >= 620
%         pause
%     end
    pause
end


