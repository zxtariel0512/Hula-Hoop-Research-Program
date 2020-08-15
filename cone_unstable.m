%% define
%this value can change
f=4.5;
w=2*pi*f;
phi0 = pi;
phidot0= w;
yes_animate = 0;
yesplot = 1;
MaxStep = 0.1;
maxTime =1000;
Rg = 1;
Rh = 7.412625;
M = 15.64;
%set 5 cones
alpha1=4;
alpha2=10;
alpha3=15;
alpha4=20;
alpha5=25;

%% plot lists
z_list=[];
Fv1_list=[];
Fv2_list=[];
Fv3_list=[];
Fv4_list=[];
Fv5_list=[];
%set z_list first
for z=-14:0.1:-4
    z_list(end+1)=z;
end

%% alpha1=5
for z=-14:0.1:-4
    Rb=abs(z)*tand(alpha1);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    normal_force=N(:,1);
    Fv1=normal_force*sind(alpha1);
    Fv1_list(end+1)=Fv1;
end

%% alpha2=10
for z=-14:0.1:-4
    Rb=abs(z)*tand(alpha2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    normal_force=N(:,1);
    Fv2=normal_force*sind(alpha2);
    Fv2_list(end+1)=Fv2;
end

%% alpha3=15
for z=-14:0.1:-4
    Rb=abs(z)*tand(alpha3);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    normal_force=N(:,1);
    Fv3=normal_force*sind(alpha3);
    Fv3_list(end+1)=Fv3;
end

%% alpha4=20
for z=-14:0.1:-4
    Rb=abs(z)*tand(alpha4);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    normal_force=N(:,1);
    Fv4=normal_force*sind(alpha4);
    Fv4_list(end+1)=Fv4;
end

%% alpha5=25
for z=-14:0.1:-4
    Rb=abs(z)*tand(alpha5);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    normal_force=N(:,1);
    Fv5=normal_force*sind(alpha5);
    Fv5_list(end+1)=Fv5;
end

%% plot
plot(z_list,Fv1_list,'r', 'Linewidth', 1.5);
hold on
plot(z_list,Fv2_list,'g', 'Linewidth', 1.5);
hold on
plot(z_list,Fv3_list,'b', 'Linewidth', 1.5);
hold on
plot(z_list,Fv4_list,'c', 'Linewidth', 1.5);
hold on
plot(z_list,Fv5_list,'m', 'Linewidth', 1.5);
hold on
plot([-14,-4],[M*980,M*980],'--', 'Linewidth', 3);
xlabel('z (f=4.5)');
ylabel('Fv (Vertical Force)');
title('Vertical Force F_v & Hoop Position (height) at a given frequence');
legend('\alpha=5','\alpha=10','\alpha=15','\alpha=20','\alpha=25', 'Location','northwest');