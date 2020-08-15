%% define
phi0 = pi;
g=980;
Rg=1;
Rh=7.412625;
M=15.64;
pi=3.1415926;

Rw=0.5;
a=Rw;
b=Rw;
alpha1=5;
alpha2=10;
alpha3=15;
alpha4=20;
alpha5=25;

% equilibrium1=[];
% stability1=[];
% equilibrium2=[];
% stability2=[];
% equilibrium3=[];
% stability3=[];
% equilibrium4=[];
% stability4=[];
% equilibrium5=[];
% stability5=[];
% color1=zeros(5,3);
% color2=zeros(10,3);
% color3=zeros(11,3);
% color4=zeros(11,3);
% color5=zeros(12,3);
sta_1 = [];
equ_1 = [];
sta_2 = [];
equ_2 = [];
sta_3 = [];
equ_3 = [];
sta_4 = [];
equ_4 = [];
sta_5 = [];
equ_5 = [];
serr_1 = [];
eerr_1 = [];
serr_2 = [];
eerr_2 = [];
serr_3 = [];
eerr_3 = [];
serr_4 = [];
eerr_4 = [];
serr_5 = [];
eerr_5 = [];

%% alpha=5
c=Rw/tand(alpha1);
f_test=[8.96,9.02, 8.47,8.53, 7.987,8, 7.47,7.53, 6.99,7.03];
z_test=[-3.1,-3.3, -4.45,-4.55, -5.1,-5.25, -7.4,-7.55, -8.7,-8.6];

for time=1:1:5
    f_mid=(f_test(time*2-1)+f_test(time*2))*0.5;
    z_mid=(z_test(time*2-1)+z_test(time*2))*0.5;
    deltaF = f_test(time*2) - f_test(time*2-1);
    deltaZ = z_test(time*2) - z_test(time*2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equilibrium error bar: y -> f; x -> z
    EquiError5 = vpa(deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid));
    eerr_1(end+1) = 0.5 * EquiError5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stability error bar:
    tangent5 = tand(alpha1);
    StaError5 = vpa(deltaYAxis(alpha1,g,Rh,Rg,Rw,c,tangent5,deltaF,deltaZ,z_mid,f_mid));  
    serr_1(end+1) = 0.5 * StaError5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base
    f=f_test(time*2-1);
    w=2*pi*f;
    z=z_test(time*2-1);
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
    Reff=Rh+Rg-Rb;
    % find local angle
    Fx=(2*Rb)/(Rw^2);
    Fy=0;
    Fz=((-2)*z)/(c^2);
    z_intersect=(Fx*Rb)/Fz+z;
    % curvature
    angle=atand(Rb/(z_intersect+abs(z)));
    k=((tand(alpha1))^2*(-(Rw^2+(tand(alpha1))^2*z^2)^(-3/2)*(tand(alpha1))^2*z^2+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)))/(1+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)*(tand(alpha1))^2*z)^(3/2);
    % equilibrium condition
    baseEqu = (w^2*Reff*sind(angle))/g;
    baseSta = (k*Reff)/(sind(angle)*tand(angle));
%     equilibrium1(end+1) = baseEqu;
%     equilibrium1(end+1) = baseEqu+EquiError5;
%     stability1(end+1) = baseSta;
%     stability1(end+1) = baseSta+StaError5;
    equ_1(end+1) = baseEqu + 0.5 * EquiError5;
    sta_1(end+1) = baseSta + 0.5 * StaError5;
end

%% alpha=10
c=Rw/tand(alpha2);
f_test=[9.52,9.60, 8.96,9.02, 8.47,8.53, 7.947,8.0, 7.49,7.53, 6.99,7.03, 6.48,6.54, 5.93,6.03, 5.46,5.52, 5.02,5.07];
z_test=[-0.55,-0.7, -0.75,-0.55, -0.65,-0.85, -0.65,-0.8, -0.95,-1.05, -1.25,-1.05, -1.55,-1.65, -1.75,-1.85, -2.45,-2.65, -5.45,-5.65];

for time=1:1:10
    f_mid=(f_test(time*2-1)+f_test(time*2))*0.5;
    z_mid=(z_test(time*2-1)+z_test(time*2))*0.5;
    deltaF = f_test(time*2) - f_test(time*2-1);
    deltaZ = z_test(time*2) - z_test(time*2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error bar: y -> f; x -> z
    EquiError10 = vpa(deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid));
    eerr_2(end+1) = 0.5 * EquiError10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stability error bar:
    tangent10 = tand(alpha2);
    StaError10 = vpa(deltaYAxis(alpha2,g,Rh,Rg,Rw,c,tangent10,deltaF,deltaZ,z_mid,f_mid));  
    serr_2(end+1) = 0.5 * StaError10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base
    f=f_test(time*2-1);
    w=2*pi*f;
    z=z_test(time*2-1);
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
    Reff=Rh+Rg-Rb;
    % find local angle
    Fx=(2*Rb)/(Rw^2);
    Fy=0;
    Fz=((-2)*z)/(c^2);
    z_intersect=(Fx*Rb)/Fz+z;
    % curvature
    angle=atand(Rb/(z_intersect+abs(z)));
    k=((tand(alpha1))^2*(-(Rw^2+(tand(alpha1))^2*z^2)^(-3/2)*(tand(alpha1))^2*z^2+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)))/(1+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)*(tand(alpha1))^2*z)^(3/2);
    % equilibrium condition
    baseEqu = (w^2*Reff*sind(angle))/g;
    baseSta = (k*Reff)/(sind(angle)*tand(angle));
%     equilibrium2(end+1) = baseEqu;
%     equilibrium2(end+1) = baseEqu+EquiError10;
%     stability2(end+1) = baseSta;
%     stability2(end+1) = baseSta+StaError10;
    equ_2(end+1) = baseEqu + 0.5 * EquiError10;
    sta_2(end+1) = baseSta + 0.5 * StaError10;
end



%% alpha=15
c=Rw/tand(alpha3);
f_test=[9.41,9.52, 8.98,9.09, 8.49,8.57, 8,8.001, 7.521,7.531, 6.95,7.01, 6.52,6.58, 6.01,6.07, 5.561,5.58, 5.01,5.07, 4.57,4.59];
z_test=[-0.2,-0.35, -0.3,-0.4, -0.3,-0.45, -0.35,-0.45, -0.4,-0.5, -0.5,-0.7, -0.5,-0.7, -0.6,-0.8, -0.9,-1, -1.1,-1.2, -2,-2.2];

for time=1:1:11
    f_mid=(f_test(time*2-1)+f_test(time*2))*0.5;
    z_mid=(z_test(time*2-1)+z_test(time*2))*0.5;
    deltaF = f_test(time*2) - f_test(time*2-1);
    deltaZ = z_test(time*2) - z_test(time*2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error bar: y -> f; x -> z
    EquiError15 = vpa(deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid));
    eerr_3(end+1) = 0.5 * EquiError15;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stability error bar:
    tangent15 = tand(alpha3);
    StaError15 = vpa(deltaYAxis(alpha3,g,Rh,Rg,Rw,c,tangent15,deltaF,deltaZ,z_mid,f_mid));
    serr_3(end+1) = 0.5 * StaError15;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base
    f=f_test(time*2-1);
    w=2*pi*f;
    z=z_test(time*2-1);
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
    Reff=Rh+Rg-Rb;
    % find local angle
    Fx=(2*Rb)/(Rw^2);
    Fy=0;
    Fz=((-2)*z)/(c^2);
    z_intersect=(Fx*Rb)/Fz+z;
    % curvature
    angle=atand(Rb/(z_intersect+abs(z)));
    k=((tand(alpha1))^2*(-(Rw^2+(tand(alpha1))^2*z^2)^(-3/2)*(tand(alpha1))^2*z^2+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)))/(1+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)*(tand(alpha1))^2*z)^(3/2);
    % equilibrium condition
    baseEqu = (w^2*Reff*sind(angle))/g;
    baseSta = (k*Reff)/(sind(angle)*tand(angle));
%     equilibrium3(end+1) = baseEqu;
%     equilibrium3(end+1) = baseEqu+EquiError15;
%     stability3(end+1) = baseSta;
%     stability3(end+1) = baseSta+StaError15;
    equ_3(end+1) = baseEqu + 0.5 * EquiError15;
    sta_3(end+1) = baseSta + 0.5 * StaError15;
end

%% alpha=20
c=Rw/tand(alpha4);
f_test=[8.94,9.02, 8.49,8.55, 8.02,8.10, 7.49,7.55, 7.01,7.08, 6.46,6.52, 6.03,6.07, 5.48,5.52, 4.98,5.05, 4.49,4.53, 4.02,4.06];
z_test=[-0.1,-0.2, -0.1,-0.3, -0.1,-0.3, -0.1,-0.3, -0.1,-0.3, -0.15,-0.35, -0.3,-0.6, -0.4,-0.65, -0.5,-0.8, -0.8,-1.0, -1.2,-1.4];

for time=1:1:11
    f_mid=(f_test(time*2-1)+f_test(time*2))*0.5;
    z_mid=(z_test(time*2-1)+z_test(time*2))*0.5;
    deltaF = f_test(time*2) - f_test(time*2-1);
    deltaZ = z_test(time*2) - z_test(time*2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error bar: y -> f; x -> z
    EquiError20 = vpa(deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid));
    eerr_4(end+1) = 0.5 * EquiError20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stability error bar:
    tangent20 = tand(alpha4);
    StaError20 = vpa(deltaYAxis(alpha4,g,Rh,Rg,Rw,c,tangent20,deltaF,deltaZ,z_mid,f_mid));
    serr_4(end+1) = 0.5 * StaError20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base
    f=f_test(time*2-1);
    w=2*pi*f;
    z=z_test(time*2-1);
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
    Reff=Rh+Rg-Rb;
    % find local angle
    Fx=(2*Rb)/(Rw^2);
    Fy=0;
    Fz=((-2)*z)/(c^2);
    z_intersect=(Fx*Rb)/Fz+z;
    % curvature
    angle=atand(Rb/(z_intersect+abs(z)));
    k=((tand(alpha1))^2*(-(Rw^2+(tand(alpha1))^2*z^2)^(-3/2)*(tand(alpha1))^2*z^2+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)))/(1+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)*(tand(alpha1))^2*z)^(3/2);
    % equilibrium condition
    baseEqu = (w^2*Reff*sind(angle))/g;
    baseSta = (k*Reff)/(sind(angle)*tand(angle));
%     equilibrium4(end+1) = baseEqu;
%     equilibrium4(end+1) = baseEqu+EquiError20;
%     stability4(end+1) = baseSta;
%     stability4(end+1) = baseSta+StaError20;
    equ_4(end+1) = baseEqu + 0.5 * EquiError20;
    sta_4(end+1) = baseSta + 0.5 * StaError20;
end

%% alpha=25
c=Rw/tand(alpha5);
f_test=[8.92,9.04, 8.49,8.52, 8.01,8.08, 7.49,7.57, 6.97,7.03, 6.52,6.58, 5.98,6.04, 5.48,5.52, 4.98,5.07, 4.49,4.53, 4.02,4.08, 3.523,3.539];
z_test=[-0.1,-0.3, 0,-0.2, -0.1,-0.3, -0.1,-0.2, -0.1,-0.25, -0.1,-0.2, -0.2,-0.4, -0.25,-0.4, -0.3,-0.4, -0.45,-0.6, -0.45,-0.8, -1.1,-1.4];

for time=1:1:12
    f_mid=(f_test(time*2-1)+f_test(time*2))*0.5;
    z_mid=(z_test(time*2-1)+z_test(time*2))*0.5;
    deltaF = f_test(time*2) - f_test(time*2-1);
    deltaZ = z_test(time*2) - z_test(time*2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error bar: y -> f; x -> z
    EquiError25 = vpa(deltaXAxis(g,Rh,Rg,Rw,c,deltaF,deltaZ,z_mid,f_mid));
    eerr_5(end+1) = 0.5 * EquiError25;
    % stability error bar:
    tangent25 = tand(alpha5);
    StaError25 = vpa(deltaYAxis(alpha5,g,Rh,Rg,Rw,c,tangent25,deltaF,deltaZ,z_mid,f_mid)); 
    serr_5(end+1) = 0.5 * StaError25;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base
    f=f_mid;
    w=2*pi*f;
    z=z_mid;
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
    Reff=Rh+Rg-Rb;
    % find local angle
    Fx=(2*Rb)/(Rw^2);
    Fy=0;
    Fz=((-2)*z)/(c^2);
    z_intersect=(Fx*Rb)/Fz+z;
    % curvature
    angle=atand(Rb/(z_intersect+abs(z)));
    k=((tand(alpha1))^2*(-(Rw^2+(tand(alpha1))^2*z^2)^(-3/2)*(tand(alpha1))^2*z^2+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)))/(1+(Rw^2+(tand(alpha1))^2*z^2)^(-1/2)*(tand(alpha1))^2*z)^(3/2);
    % equilibrium condition
    baseEqu = (w^2*Reff*sind(angle))/g;
    baseSta = (k*Reff)/(sind(angle)*tand(angle));
%     equilibrium5(end+1) = baseEqu-0.5*EquiError25;
%     equilibrium5(end+1) = baseEqu+0.5*EquiError25;
%     stability5(end+1) = baseSta-0.5*StaError25;
%     stability5(end+1) = baseSta+0.5*StaError25;
    equ_5(end+1) = baseEqu + 0.5 * EquiError25;
    sta_5(end+1) = baseSta + 0.5 * StaError25;
end

% %% scatter alpha=5
% for time=1:1:5
%     rectangle('Position',[equilibrium1(2*time-1) stability1(2*time) abs(equilibrium1(2*time)-equilibrium1(2*time-1)) abs(stability1(2*time-1)-stability1(2*time))],'Edgecolor','b');
% end
% hold on
% 
% %% scatter alpha=10
% for time=1:1:10
%     rectangle('Position',[equilibrium2(2*time-1) stability2(2*time) abs(equilibrium2(2*time)-equilibrium2(2*time-1)) abs(stability2(2*time-1)-stability2(2*time))],'Edgecolor','m');
% end
% hold on
% 
% %% scatter alpha=15
% for time=1:1:11
%     rectangle('Position',[equilibrium3(2*time-1) stability3(2*time) abs(equilibrium3(2*time)-equilibrium3(2*time-1)) abs(stability3(2*time-1)-stability3(2*time))],'Edgecolor','g');
% end
% hold on
% 
% %% scatter alpha=20
% for time=1:1:11
%     rectangle('Position',[equilibrium4(2*time-1) stability4(2*time) abs(equilibrium4(2*time)-equilibrium4(2*time-1)) abs(stability4(2*time-1)-stability4(2*time))],'Edgecolor','c');
% end
% hold on
% 
% %% scatter alpha=25
% for time=1:1:12
%     rectangle('Position',[equilibrium5(2*time-1) stability5(2*time) abs(equilibrium5(2*time)-equilibrium5(2*time-1)) abs(stability5(2*time-1)-stability5(2*time))],'Edgecolor','r');
% end
% hold on
% % equilibrium5
% % stability5

%% scatter
errorbar(equ_1, sta_1, eerr_1, 'horizontal', 'bo', 'Linewidth', 1.3), hold on
errorbar(equ_1, sta_1, serr_1, 'vertical', 'bo', 'Linewidth', 1.3), hold on
errorbar(equ_2, sta_2, eerr_2, 'horizontal', 'mo', 'Linewidth', 1.3), hold on
errorbar(equ_2, sta_2, serr_2, 'vertical', 'mo', 'Linewidth', 1.3), hold on
errorbar(equ_3, sta_3, eerr_3, 'horizontal', 'go', 'Linewidth', 1.3), hold on
errorbar(equ_3, sta_3, serr_3, 'vertical', 'go', 'Linewidth', 1.3), hold on
errorbar(equ_4, sta_4, eerr_4, 'horizontal', 'co', 'Linewidth', 1.3), hold on
errorbar(equ_4, sta_4, serr_4, 'vertical', 'co', 'Linewidth', 1.3), hold on
errorbar(equ_5, sta_5, eerr_5, 'horizontal', 'ro', 'Linewidth', 1.3), hold on
errorbar(equ_5, sta_5, serr_5, 'vertical', 'ro', 'Linewidth', 1.3), hold on


xlabel("Equilibrium Tester: w^2 * R_{eff} * sin(\alpha) / g");
ylabel("Stability Tester: k * R_{eff} / (sin(\a)tan(\a))");
title('General Condition Testing Graph for 5 Hyperboloids');
% legend("alpha=5","alpha=10","alpha=15","alpha=20","alpha=25");