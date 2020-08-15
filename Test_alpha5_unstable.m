%% simulation part

% data
alpha=5;
phi0 = pi;
yes_animate = 0;
yesplot = 1;
MaxStep = 0.1;
maxTime =400;
g=980;
Rg = 1;
Rh = 7.412625;
M = 15.64;
Rw=0.5;
a=Rw;
b=Rw;

% graph
f_list=[];
z_list=[];
f_defined=[];
f_udefined=[];
%normal_force_list=[];

% calculate
c=Rw/tand(alpha);
%for loop for f
for f=1:0.001:10
    w=2*pi*f;
    phidot0= w;
    %for loop for z
    for z=-10:0.01:0
        Rb=sqrt(Rw^2+(Rw^2*z^2)/c^2);
%         N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
%         normal_force=N(:,1);
        normal_force=abs(-M*w^2*(Rh+Rg-Rb));
%         disp([num2str(normal_force)])
%         normal_force(end+1)=normal_force;
         
        vertical_force=normal_force*sind(angle);
        %find next vertical force
        z_next=z+0.1;
        Rb_next=sqrt(Rw^2+(Rw^2*(z_next)^2)/c^2);
%         N_next=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb_next,M,maxTime);
%         normal_force_next=N_next(:,1);
        normal_force_next=abs(-M*w^2*(Rh+Rg-Rb_next));
        Fx_next=(2*Rb_next)/(Rw^2);
        Fy_next=0;
        Fz_next=((-2)*z_next)/(c^2);
        z_intersect_next=(Fx_next*Rb_next)/Fz_next+z_next;
        if z_intersect_next>=0
            angle_next=atand(Rb_next/(z_intersect_next+abs(z_next)));
        end
        if z_intersect_next<0
            angle_next=atand(Rb_next/(abs(z_next)-abs(z_intersect_next)));
        end
        vertical_force_next=normal_force_next*sind(angle_next);
        %find previous vertical force
        z_previous=z-0.1;
        Rb_previous=sqrt(Rw^2+(Rw^2*(z_previous)^2)/c^2);
%         N_previous=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb_previous,M,maxTime);
%         normal_force_previous=N_previous(:,1);
        normal_force_previous=abs(-M*w^2*(Rh+Rg-Rb_previous));
        Fx_previous=(2*Rb_previous)/(Rw^2);
        Fy_previous=0;
        Fz_previous=((-2)*z_previous)/(c^2);
        z_intersect_previous=(Fx_previous*Rb_previous)/Fz_previous+z_previous;
        if z_intersect_previous>=0
            angle_previous=atand(Rb_previous/(z_intersect_previous+abs(z_previous)));
        end
        if z_intersect_previous<0
            angle_previous=atand(Rb_previous/(abs(z_previous)-abs(z_intersect_previous)));
        end
        vertical_force_previous=normal_force_previous*sind(angle_previous);
        %check if this is an equilibrium pt
        if abs(vertical_force-M*g)<=3000
            %check stable
            if vertical_force_next<M*g && vertical_force_previous>M*g
                if length(f_defined)==0
                    f_defined(end+1)=f;
                    z_list(end+1)=z;
                end
                if f_defined(end)==f
                    z_compare=z_list(end);
                    Rb_compare=sqrt(Rw^2+(Rw^2*(z_compare)^2)/c^2);
    %                 N_compare=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb_compare,M,maxTime);
    %                 normal_force_compare=N_compare(:,1);
                    normal_force_compare=abs(-M*w^2*(Rh+Rg-Rb_compare));
                    Fx_compare=(2*Rb_compare)/(Rw^2);
                    Fy_compare=0;
                    Fz_compare=((-2)*z_compare)/(c^2);
                    z_intersect_compare=(Fx_compare*Rb_compare)/Fz_compare+z_compare;
                    if z_intersect_compare>=0
                        angle_compare=atand(Rb_compare/(z_intersect_compare+abs(z_compare)));
                    end
                    if z_intersect_compare<0
                        angle_compare=atand(Rb_compare/(abs(z_compare)-abs(z_intersect_compare)));
                    end
                    vertical_force_compare=normal_force_compare*sind(angle_compare);
                    if abs(vertical_force-M*g)<abs(vertical_force_compare-M*g)
                        z_list(end)=z;
                    end
                end
                if f_defined(end)~=f
                    f_defined(end+1)=f;
                    z_list(end+1)=z;
                end
            end
            %check unstable
            if vertical_force_next>M*g && vertical_force_previous<M*g
                if length(f_udefined)==0
                    f_udefined(end+1)=f;
                    z_ulist(end+1)=z;
                end
                if f_udefined(end)==f
                    z_compare=z_ulist(end);
                    Rb_compare=sqrt(Rw^2+(Rw^2*(z_compare)^2)/c^2);
    %                 N_compare=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb_compare,M,maxTime);
    %                 normal_force_compare=N_compare(:,1);
                    normal_force_compare=abs(-M*w^2*(Rh+Rg-Rb_compare));
                    Fx_compare=(2*Rb_compare)/(Rw^2);
                    Fy_compare=0;
                    Fz_compare=((-2)*z_compare)/(c^2);
                    z_intersect_compare=(Fx_compare*Rb_compare)/Fz_compare+z_compare;
                    if z_intersect_compare>=0
                        angle_compare=atand(Rb_compare/(z_intersect_compare+abs(z_compare)));
                    end
                    if z_intersect_compare<0
                        angle_compare=atand(Rb_compare/(abs(z_compare)-abs(z_intersect_compare)));
                    end
                    vertical_force_compare=normal_force_compare*sind(angle_compare);
                    if abs(vertical_force-M*g)<abs(vertical_force_compare-M*g)
                        z_ulist(end)=z;
                    end
                end
                if f_udefined(end)~=f
                    f_udefined(end+1)=f;
                    z_ulist(end+1)=z;
                end
            end
        end
    end
end


%% plot part

%normal_force_list
% simulation
plot_list=[];
for f=1:0.1:10
    f_list(end+1)=f;
    plot_list(end+1)=0;
end
% z_list
plot(f_defined,z_list,'r', 'Linewidth', 1.5);
hold on
plot(f_udefined,z_ulist,'--', 'Linewidth', 1.5);
hold on
plot(f_list,plot_list,'k');
hold on

%% experiment part

f_test=[8.96,9.02, 8.47,8.53, 7.987,8, 7.47,7.53, 6.99,7.03];
z_test=[-3.1,-3.3, -4.45,-4.55, -5.1,-5.25, -7.4,-7.55, -8.7,-8.6];
for time=1:1:5
%     disp(num2str(f_test(2*time)));
%     x=[f_test(2*time-1),f_test(2*time)]
%     y=[z_test(2*time-1),z_test(2*time)]
%     disp(num2str(f_test(2*time-1)));
    rectangle('Position',[f_test(2*time-1) z_test(2*time) abs(f_test(2*time)-f_test(2*time-1)) abs(z_test(2*time-1)-z_test(2*time))], 'Linewidth', 1);
    axis([0 10 -10 0]);
end
xlabel("f (\alpha = 5)");
ylabel("stable/unstable point");
title('Equilibrium position & Frequence for \alpha = 5 hyperboloid');