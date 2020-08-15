%% define
%this f can change
f=4;
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

%special for hyperboloid
Rw=0.5;
a=Rw;
b=Rw;
alpha1=5;
alpha2=10;
alpha3=15;
alpha4=20;
alpha5=25;

%% tangent plane
%F(x,y,z)=x^2/Rw^2+y^2/Rw^2-z^2/c^2-1=0
%Fx=2*x/Rw^2
%Fy=2*y/Rw^2
%Fz=-2*z/c^2
%point:(Rb,0,-z) down, this can change, not a big deal

%% alpha = 5
c1=Rw/tand(alpha1);
vertical_force_list1=[];
z_list_negative=[];
%from -10 to 10, with z=0 on the wrist
for z=-10:0.1:0
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c1^2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    current_normal_force=N(:,1);
    if z<=0
        %calculate the tangent plane
        Fx=(2*Rb)/(Rw^2);
        Fy=0;
        Fz=((-2)*z)/(c1^2);
        %so the equation for tangent plane is:Fx*(x-Rb)+Fy*y+Fz*(z-z')=0
        %now let x=0,y=0, get the intersection point
        z_intersect=(Fx*Rb)/Fz+z;
        if z_intersect>=0
            current_angle=atand(Rb/(z_intersect+abs(z)));
        end
        if z_intersect<0
            current_angle=atand(Rb/(abs(z)-abs(z_intersect)));
        end
        current_vertical_force=current_normal_force*sind(current_angle);
        vertical_force_list1(end+1)=current_vertical_force;
        z_list_negative(end+1)=z;
    end
    
%     current_angle=atand(abs(z)/Rb);
%     angle_list(end+1)=current_angle;
%     current_vertical_force=current_normal_force*sind(current_angle);
%     vertical_force_list(end+1)=current_vertical_force;
    
end

%% alpha = 10
c2=Rw/tand(alpha2);
vertical_force_list2=[];
%from -10 to 10, with z=0 on the wrist
for z=-10:0.1:0
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c2^2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    current_normal_force=N(:,1);
    if z<=0
        %calculate the tangent plane
        Fx=(2*Rb)/(Rw^2);
        Fy=0;
        Fz=((-2)*z)/(c2^2);
        %so the equation for tangent plane is:Fx*(x-Rb)+Fy*y+Fz*(z-z')=0
        %now let x=0,y=0, get the intersection point
        z_intersect=(Fx*Rb)/Fz+z;
        if z_intersect>=0
            current_angle=atand(Rb/(z_intersect+abs(z)));
        end
        if z_intersect<0
            current_angle=atand(Rb/(abs(z)-abs(z_intersect)));
        end
        current_vertical_force=current_normal_force*sind(current_angle);
        vertical_force_list2(end+1)=current_vertical_force;
    end
    
%     current_angle=atand(abs(z)/Rb);
%     angle_list(end+1)=current_angle;
%     current_vertical_force=current_normal_force*sind(current_angle);
%     vertical_force_list(end+1)=current_vertical_force;
    
end

%% alpha = 15
c3=Rw/tand(alpha3);
vertical_force_list3=[];
%from -10 to 10, with z=0 on the wrist
for z=-10:0.1:0
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c3^2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    current_normal_force=N(:,1);
    if z<=0
        %calculate the tangent plane
        Fx=(2*Rb)/(Rw^2);
        Fy=0;
        Fz=((-2)*z)/(c3^2);
        %so the equation for tangent plane is:Fx*(x-Rb)+Fy*y+Fz*(z-z')=0
        %now let x=0,y=0, get the intersection point
        z_intersect=(Fx*Rb)/Fz+z;
        if z_intersect>=0
            current_angle=atand(Rb/(z_intersect+abs(z)))
        end
        if z_intersect<0
            current_angle=atand(Rb/(abs(z)-abs(z_intersect)))
        end
        current_vertical_force=current_normal_force*sind(current_angle);
        vertical_force_list3(end+1)=current_vertical_force;
    end
    
%     current_angle=atand(abs(z)/Rb);
%     angle_list(end+1)=current_angle;
%     current_vertical_force=current_normal_force*sind(current_angle);
%     vertical_force_list(end+1)=current_vertical_force;
    
end

%% alpha = 20
c4=Rw/tand(alpha4);
vertical_force_list4=[];
%from -10 to 10, with z=0 on the wrist
for z=-10:0.1:0
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c4^2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    current_normal_force=N(:,1);
    if z<=0
        %calculate the tangent plane
        Fx=(2*Rb)/(Rw^2);
        Fy=0;
        Fz=((-2)*z)/(c4^2);
        %so the equation for tangent plane is:Fx*(x-Rb)+Fy*y+Fz*(z-z')=0
        %now let x=0,y=0, get the intersection point
        z_intersect=(Fx*Rb)/Fz+z;
        if z_intersect>=0
            current_angle=atand(Rb/(z_intersect+abs(z)));
        end
        if z_intersect<0
            current_angle=atand(Rb/(abs(z)-abs(z_intersect)));
        end
        current_vertical_force=current_normal_force*sind(current_angle);
        vertical_force_list4(end+1)=current_vertical_force;
    end
    
%     current_angle=atand(abs(z)/Rb);
%     angle_list(end+1)=current_angle;
%     current_vertical_force=current_normal_force*sind(current_angle);
%     vertical_force_list(end+1)=current_vertical_force;
    
end

%% alpha = 25
c5=Rw/tand(alpha5);
vertical_force_list5=[];
%from -10 to 10, with z=0 on the wrist
for z=-10:0.1:0
    Rb=sqrt(Rw^2+(Rw^2*z^2)/c5^2);
    N=solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime);
    current_normal_force=N(:,1);
    if z<=0
        %calculate the tangent plane
        Fx=(2*Rb)/(Rw^2);
        Fy=0;
        Fz=((-2)*z)/(c5^2);
        %so the equation for tangent plane is:Fx*(x-Rb)+Fy*y+Fz*(z-z')=0
        %now let x=0,y=0, get the intersection point
        z_intersect=(Fx*Rb)/Fz+z;
        if z_intersect>=0
            current_angle=atand(Rb/(z_intersect+abs(z)));
        end
        if z_intersect<0
            current_angle=atand(Rb/(abs(z)-abs(z_intersect)));
        end
        current_vertical_force=current_normal_force*sind(current_angle);
        vertical_force_list5(end+1)=current_vertical_force;
    end
    
%     current_angle=atand(abs(z)/Rb);
%     angle_list(end+1)=current_angle;
%     current_vertical_force=current_normal_force*sind(current_angle);
%     vertical_force_list(end+1)=current_vertical_force;
    
end
%% plot
plot(z_list_negative,vertical_force_list1,'r', 'Linewidth', 1.5);
hold on
plot(z_list_negative,vertical_force_list2,'g', 'Linewidth', 1.5);
hold on
plot(z_list_negative,vertical_force_list3,'b', 'Linewidth', 1.5);
hold on
plot(z_list_negative,vertical_force_list4,'c', 'Linewidth', 1.5);
hold on
plot(z_list_negative,vertical_force_list5,'k', 'Linewidth', 1.5);
hold on
plot([-10,0],[M*980,M*980],'--', 'Linewidth', 3);
xlabel('z (origin at the waist) (f=4hz R_w=0.5)','FontWeight','bold','FontSize',12');
ylabel('Vertical Force F_v','FontSize',12);
title('Hyperboloid: Vertical Force F_v & Hoop Position (height) at a given frequence');
legend('alpha=5','alpha=10','alpha=15','alpha=20','alpha=25', 'Location', 'northwest');


