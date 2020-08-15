
function [normal_force] = solve_hulahoop_normalForce(phi0,phidot0,w,Rg,Rh,Rb,M,maxTime,MaxStep)
if nargin< 8
    maxTime = 10;
end
% if nargin < 9
%     yesplot = 0;
% end
% if nargin < 10
%     yes_animate = 0;
% end
if nargin< 11 || MaxStep==0
    MaxStep = 0.05;
end

% if yes_animate>1
%     body_frame = 1;
% else
%     body_frame = 0;
% end

% colors = {[0    0.4470    0.7410],[0.4660    0.6740    0.1880],[0.6350    0.0780    0.1840],...
%     [0.8500    0.3250    0.0980],[0.9290    0.6940    0.1250],[0.4940    0.1840    0.5560],[0,1,1],[1,0,1],[1,1,1]};
% 
% cmap = summer(20);
% green_colors = {cmap(1,:),cmap(4,:),cmap(7,:)};

% video_output = ['~/Desktop/hulahoop bodyframe_ ' num2str(body_frame) ' w_' num2str(w) ' Rg_' num2str(Rg)  ' Rh_' num2str(Rh)  ' Rb_' num2str(Rb) ' M_' num2str(M)  ' phi_' num2str(phi0,3) ' phidot_' num2str(phidot0,3) ' maxTime_' num2str(maxTime,3)];
% filename = video_output;
% filename2 = [video_output ' plot phi vid'];
% filename3 = [video_output ' plot norm vid'];


% disp(['w = ' num2str(w)])
% disp(['Rg = ' num2str(Rg)])
% disp(['Rh = ' num2str(Rh)])
% disp(['Rb = ' num2str(Rb)])
% disp(['M = ' num2str(M)])


options = odeset('RelTol',1e-7,'MaxStep',MaxStep);
tspan = [0:maxTime/w];
 if length(tspan) > 4

XYlim = round(Rh + Rg + Rb +.5);
num_points = 4*Rh;
pointSz = 10*XYlim;


init_phi = [phi0; phidot0];

%init_psi = (Rh - Rb)*init_phi/Rh;

%% solve ODE for Phi
sol_phi = ode45(@(t,y) phiODE_v2(t,y,w,Rg,Rb,Rh),tspan,init_phi,options);

tt = sol_phi.x;
phi = sol_phi.y(1,:);
phidot = sol_phi.y(2,:);
phidotdot = w^2*sin(w*tt-phi)*Rg/(2*(Rb-Rh));

%% solve ODE for Psi
%psi = psi_eq(tt,tt,phi,vernum);
%psi = psi_eq_v2(tt,tt,phi,Rg,Rb,Rh)
%psidot = psi_eq(tt,tt,phidot,vernum);
%psidot = psi_eq_v2(tt,tt,phidot,Rg,Rb,Rh);

%sol_psi = ode45(@(t,y) psiODE(t,y,tt1,phi_t,vernum),tspan,psi_0,options);
%old_tt2 = sol_psi.x;
%old_psi_t = sol_psi.y(1,:);
%old_psi_dt = sol_psi.y(2,:);

%% Solve for Lambda

%lmbda = lambda_eq(tt,tt,phi,vernum);
%lmbda = lambda_eq_v2(tt,tt,phi,Rg,Rb,Rh,M,w)

%% Normal force equal to hoop_mass times normal component of hoop acceleration
normal_force = M*(Rh-Rb)*phidot.^2 - M*(w^2)*Rg*(cos(w*tt).*cos(phi)+sin(w*tt).*sin(phi));
%normal_force = M*(Rh-Rb)*phidot.^2 + M*(w^2)*Rg*(cos(w*tt).*cos(phi)-sin(w*tt).*sin(phi)); %from Belyokov paper

% fric_force = -M*(Rh-Rb)*phidotdot - M*(w^2)*Rg*(cos(w*tt).*sin(phi)-sin(w*tt).*cos(phi));
%% Convert to Cartesian and animate, there is certainly a better way to animate
% if yes_animate==1
%     %close all
%     
%     max_col = length(colors);
%     fignum = round(1+1000*init_phi(1));
%     figure(fignum)
%     set(gcf,'Position',[300 30 880 680])
%     xlim([-10 10])
%     ylim([-10 10])
%     axis equal
%     init_points = [0:2*pi/num_points:2*pi];
%     
%     for i=1:length(tt)
%         cla
%         hold off
%         for j=i%max(i-(max_col-1),1):i
%             t = tt(j);
%             phi_step = phi(j);
%             psi_step = psi(j);
%             dot_psi_step = psidot(j);
%             
%             if body_frame
%                 X_body_cm = 0;
%                 Y_body_cm = 0;
%             else
%                 X_body_cm = Rg*cos(w*t);
%                 Y_body_cm = Rg*sin(w*t);
%             end
%             
%             
%             X_contact = X_body_cm + Rb*cos(phi_step);
%             Y_contact = Y_body_cm + Rb*sin(phi_step);
%             
%             X_hoop_cm = X_body_cm - (Rh-Rb)*cos(phi_step);
%             Y_hoop_cm = Y_body_cm - (Rh-Rb)*sin(phi_step);
%             
%             X_point = X_hoop_cm + Rh*cos(psi_step+init_points);
%             Y_point = Y_hoop_cm + Rh*sin(psi_step+init_points);
%             
%             
%             centers = [X_body_cm,Y_body_cm; X_hoop_cm, Y_hoop_cm];
%             radii = [Rb,Rh];
%             viscircles([0,0],Rg,'Color','k','LineStyle',':');
%             hold on
%             viscircles(centers(1,:),radii(1),'Color',colors{1},'LineStyle','-');
%             %plot([X_contact,X_hoop_cm],[Y_contact,Y_hoop_cm],'k-','linewidth',3)
%             viscircles(centers(2,:),radii(2),'Color',colors{2},'LineStyle','-');
%             %scatter(X_point(1:3:end),Y_point(1:2:end),pointSz/2,colors{3},'filled')
%             scatter(X_point(1:2:end),Y_point(1:2:end),pointSz/2,green_colors{1},'filled')
%             scatter(X_point(2:2:end),Y_point(2:2:end),pointSz/2,green_colors{3},'filled')
%             scatter(X_body_cm,Y_body_cm,80,colors{1},'filled')
%             scatter(X_hoop_cm,Y_hoop_cm,80,colors{2},'filled')
%             %plot([0,3*X_body_cm],[0,3*Y_body_cm],'Color',colors{4},'linewidth',1.5)
%             %plot(0,0,'xk','MarkerSize',10,'LineWidth',2)
%             scatter(0,0,80,'k','filled')
            
%             %% visualize phi
%             %plot([X_body_cm,X_body_cm+Rb],[Y_body_cm,Y_body_cm],'-','Color',colors{3},'linewidth',1.5)
%             %plot([X_contact,X_body_cm],[Y_contact,Y_body_cm],'-','Color',colors{3},'linewidth',1.5)
%             plot([X_body_cm,X_body_cm+Rb],[Y_body_cm,Y_body_cm],'-','Color',colors{4},'linewidth',1.5)
%             plot([X_contact,X_body_cm],[Y_contact,Y_body_cm],'-','Color',colors{4},'linewidth',1.5)
            
%             %% visualize psi
%             %plot([X_hoop_cm,X_hoop_cm+Rh],[Y_hoop_cm,Y_hoop_cm],'-','Color',colors{4},'linewidth',1.5)
%             %plot([X_hoop_cm,X_point(1)],[Y_hoop_cm,Y_point(1)],'-','Color',colors{4},'linewidth',1.5)
%             plot([X_hoop_cm,X_hoop_cm+Rh],[Y_hoop_cm,Y_hoop_cm],'-','Color',green_colors{2},'linewidth',1.5)
%             plot([X_hoop_cm,X_point(1)],[Y_hoop_cm,Y_point(1)],'-','Color',green_colors{2},'linewidth',1.5)
            
%             %% visualize theta
%             %plot([0,X_body_cm],[0,Y_body_cm],'-','Color','k','linewidth',1.5)
%             %plot([0,Rg],[0,0],'-','Color','k','linewidth',1.5)
%             plot([0,X_body_cm],[0,Y_body_cm],'-','Color',colors{6},'linewidth',1.5)
%             plot([0,Rg],[0,0],'-','Color',colors{6},'linewidth',1.5)
%         end
        
%         hold off
%         xlim([-XYlim XYlim])
%         ylim([-XYlim XYlim])
%         %box on
%         %set(gca,'YTickLabel',[]);
%         %set(gca,'XTickLabel',[]);
%         axis off
%         shg;
%         %a2 = getframe;
%         %F = im2frame(gcf);
%         Mov(i) = getframe(gca);
%         %pause(.0005)
%         
%     end
    
%     figure(fignum)
%     if body_frame
%         title('body frame')
%     else
%         title('lab frame')
%     end
% end



%if yes_animate == 2
    %for i=1:length(tt)
        %t = tt(i);
        %phi_step = phi(i);
        %figure(10)
        %hold off
        %plot(tt,phi-w*tt,'-','LineWidth',2.5,'Color',colors{1})
        %hold on
        %plot(t,phi_step-w*t,'.','MarkerSize',20,'Color',colors{3})
        %hold off
        %ylabel('\phi - \theta')
        %xlabel('t')
        %set(gca,'FontSize',16)
        %set(gca,'LineWidth',2)
        %set(gcf,'Position',[20 20 600 200])
        %Mov2(i) = getframe(gcf);
    %end
    
%end

%if yes_animate == 3
    %for i=1:length(tt)
        %t = tt(i);
        %norm_step = normal_force(i);
        %figure(10)
        %hold off
        %plot(tt,normal_force,'-','LineWidth',2.5,'Color',colors{1})
        %hold on
        %plot(t,norm_step,'.','MarkerSize',20,'Color',colors{3})
        %hold off
        %ylabel('Normal Force')
        %xlabel('t')
        %set(gca,'FontSize',16)
        %set(gca,'LineWidth',2)
        %set(gcf,'Position',[20 20 600 200])
        %Mov3(i) = getframe(gcf);
    %end
    
%end


% figure
% if body_frame
%     title('body frame')
% else
%     title('lab frame')
% end
% axes('Position',[0 0 1 1])
% movie(Mov,1)

%if yes_animate== 1
    %v = VideoWriter(filename);
    %open(v)
    %writeVideo(v,Mov)
    %close(v)
%end


%if yes_animate == 2
    %v = VideoWriter(filename2);
    %open(v)
    %writeVideo(v,Mov2)
    %close(v)
%end

%if yes_animate == 3
    %v = VideoWriter(filename3);
    %open(v)
    %writeVideo(v,Mov3)
    %close(v)
%end

%if yesplot
    %figure(101)
    %plot(tt,phi,'-b')
    %ylabel('phi')
    %xlabel('t')
    
    %figure(102)
    %plot(tt,phidot,'-b')
    %ylabel('phi dot')
    %xlabel('t')
    
    %figure(102)
    %plot(tt,lmbda)
    %ylabel('lambda aka torque')
    %xlabel('t')
    
    %figure(103)
    %plot(tt,normal_force)
    %ylabel('normal force')
    %xlabel('t')
    
    %figure(104)
    %plot(tt,fric_force)
    %ylabel('friction force')
    %xlabel('t')
    
    %figure(105)
    %plot(tt,phi-w*tt)
    %ylabel('phi - theta')
    %xlabel('t')
    
    
%     figure(10)
%     plot(tt,phi-w*tt,'LineWidth',2.5)
%     ylabel('\phi - \theta')
%     xlabel('t')
%     set(gca,'FontSize',16)
%     set(gca,'LineWidth',2)
%     set(gcf,'Position',[20 20 600 200])
    %print(gcf,[filename ' plot phi.png'],'-dpng')
    
%     figure(11)
%     plot(tt,normal_force,'LineWidth',2.5)
%     ylabel('Normal Force')
%     xlabel('t')
%     set(gca,'FontSize',16)
%     set(gca,'LineWidth',2)
%     set(gcf,'Position',[20 20 600 200])
    %print(gcf,[filename ' plot normal.png'],'-dpng')
    
end

% else
%     warning('maxTime must be at least 5 times larger than w')
% 
% end
%phi