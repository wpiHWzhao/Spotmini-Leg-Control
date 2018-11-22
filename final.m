% The contorl of the spotmini leg by Haowei Zhao using Computed Torque
% control law.

clc
clear;
close all;
fprintf('I am going to use Computed Torque Controller. It can control the motion of the leg. Meanwhile, by mannually cancelling the M C G matrix, the simulation would be less computationally expensive and more accurate.');
global l1 l2 lc1 lc2 m1 m2  I1 I2 g x0 p1 p2

% The parameter of Pelican Robot
l1=0.26;
l2=0.26;
lc1=0.0983;
lc2=0.0229;
m1=6.5225;
m2=2.0458;
I1=0.1213;
I2=0.0116;
g=9.81;

% Some referance point that would help to track the trajectory
p1=[0.1,-0.3];
p2=[-0.1,-0.3];
v1=[0,0];
v2=[0,0];

%Implement Inverse Kinematic to calculate initial position
Q1=InverseKinematics(p1);
Q2=InverseKinematics(p2);

x0=[Q1(1,:),Q1(2,:),v1];

% The main function 
Pelicanleg

function []=Pelicanleg()
global l1 l2 lc1 lc2 m1 m2 I1 I2 g x0 p1 p2

tf=15;

global torque
torque=[];

% Solve the problem with ode45
[T,X]=ode45(@(t,x)LegODE(t,x),[0 tf],x0);
v=size(X);
v=v(1,1);

% Plot the data
figure('Name','The theta1 of the final');
plot(T,X(:,1),'-r');
title('The theta1 of the final');
print('-dpng','The theta1 of the final.png');
hold on

figure('Name','The theta2 of the final');
plot(T,X(:,2),'--r');
title('The theta2 of the final');
print('-dpng','The theta2 of the final.png');
hold on

figure('Name','The dtheta1 of the final');
plot(T,X(:,3),'-r');
title('The dtheta1 of the final');
set(gca,'XTick',0:1:15);
print('-dpng','The dtheta1 of the final.png');
hold on

figure('Name','The dtheta2 of the final');
plot(T,X(:,4),'--r');
title('The dtheta2 of the final');
set(gca,'XTick',0:1:15);
print('-dpng','The dtheta2 of the final.png');
hold on

figure('Name','Input of controller of final');
plot(T, torque(1,1:size(T,1)),'-' );
hold on
plot(T, torque(2,1:size(T,1)),'r--');
title('Input of controller of final');
print('-dpng','Input of controller of final.png');

z=1;
Xpath=[];
Ypath=[];
while z<v+1
Xpath=[Xpath,l1*cos(X(z,1))+l2*cos(X(z,1)-X(z,2))];
Ypath=[Ypath,l1*sin(-X(z,1))-l2*sin(X(z,1)-X(z,2))];
z=z+1;
end

figure('Name','Path of end-effector');
plot(Xpath,Ypath,'-r');
title('Path of end-effector');
axis([-0.15 0.15 -0.35 -0.15]);
print('-dpng','Path of end-effector,png');
hold on

figure('Name','The x position of the end-effector');
plot(T,l1*cos(X(:,1))+l2*cos(X(:,1)-X(:,2)),'-r');
title('The x position of the end-effector');
print('-dpng','The x position of the end-effector.png');
hold on

figure('Name','The y position of the end-effector');
plot(T,l1*sin(-X(:,1))-l2*sin(X(:,1)-X(:,2)),'--r');
title('The y position of the end-effector');
print('-dpng','The y position of the end-effector.png');
hold on 

z=1;
Pdot=[];
while z<v+1
Jacop=[-l2*sin(-X(z,1)+X(z,2))-l1*sin(-X(z,1)),-l2*sin(-X(z,1)+X(z,2));
       l2*cos(-X(z,1)+X(z,2))+l1*cos(-X(z,1)),l2*cos(-X(z,1)+X(z,2))];
Pdot=[Pdot,Jacop*[-X(z,3);X(z,4)]];
z=z+1;
end

figure('Name','The horizontal velocity of end-effector');
plot(T,Pdot(1,:),'-r');
title('The horizontal velocity of end-effector')
set(gca,'XTick',0:1:15);
print('-dpng','The horizontal velocity of end-effector.png');
hold on

figure('Name','The vertical velocity of end-effector');
plot(T,Pdot(2,:),'--r');
title('The vertical velocity of end-effector')
set(gca,'XTick',0:1:15);
print('-dpng','The vertical velocity of end-effector.png');
hold on

% Generate the animation
filename = 'final.gif';
h = figure;
j=1;
XP=[];
YP=[];
while  j<v+1
    XP=[XP,l1*cos(X(j,1))+l2*cos(X(j,1)-X(j,2))];
    YP=[YP,l1*sin(-X(j,1))-l2*sin(X(j,1)-X(j,2))];
    plot([0,l1*cos(X(j,1)),l1*cos(X(j,1))+l2*cos(X(j,1)-X(j,2))],[0,l1*sin(-X(j,1)),l1*sin(-X(j,1))-l2*sin(X(j,1)-X(j,2))],'-ob',0.1,-0.3,'-or',XP,YP,'-r');
    xlabel('x position')
    ylabel('y position')
    axis([-0.5 0.5 -0.6 0.1])
    drawnow
    j=j+1;
    ax = gca;
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    marg = 30;
    rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    F = getframe(gca,rect);
    ax.Units = 'normalized';
    frame=getframe;
    im=frame2im(frame);
    [I,map]=rgb2ind(im,256);
if j==2
     imwrite(I,map,'Final.gif','gif','Loopcount',1,'DelayTime',0.05);
else
     imwrite(I,map,'Final.gif','gif','WriteMode','append','DelayTime',0.05);  
end
end

    function dx= LegODE(t,x)
       
        global pv tau ddtheta_d M C G
        
        % First part of the trajectory - horizontal line
        if t<=10
        pv=[p1(1,1)-0.02*t,p1(1,2)];
        if pv(1,1)>-0.1
          Qphase1=InverseKinematics(pv);
        end
        if pv(1,1)<=-0.1
          Qphase1=InverseKinematics(p2); 
        end
        theta_d=Qphase1;
        Jaco=[-l2*sin(-Qphase1(1,1)+Qphase1(2,1))-l1*sin(-Qphase1(1,1)),-l2*sin(-Qphase1(1,1)+Qphase1(2,1));
              l2*cos(-Qphase1(1,1)+Qphase1(2,1))+l1*cos(-Qphase1(1,1)),l2*cos(-Qphase1(1,1)+Qphase1(2,1))];
        dtheta_d=Jaco\[-0.02;0];% Don't use inv(Jaco), which would increase error. 
        ddtheta_d=[0;0];
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control1(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
        end
        
       % Second part of trajectory - vertical line
        if t>10 && t<=11
          a=0.1;
          pv=[p2(1,1),p2(1,2)+(1/2)*a*(t-10)^2];
          Qphase2=InverseKinematics(pv);
          theta_d=Qphase2;
          Jaco=[-l2*sin(-Qphase2(1,1)+Qphase2(2,1))-l1*sin(-Qphase2(1,1)),-l2*sin(-Qphase2(1,1)+Qphase2(2,1));
                 l2*cos(-Qphase2(1,1)+Qphase2(2,1))+l1*cos(-Qphase2(1,1)),l2*cos(-Qphase2(1,1)+Qphase2(2,1))];
          dtheta_d=Jaco\[0;a*(t-10)];
          JacoDot=[-l2*cos(-Qphase2(1,1)+Qphase2(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*cos(-Qphase2(1,1))*(dtheta_d(1,1)),-l2*cos(-Qphase2(1,1)+Qphase2(2,1))*(dtheta_d(1,1)+dtheta_d(2,1));
                   -l2*sin(-Qphase2(1,1)+Qphase2(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*sin(-Qphase2(1,1))*(dtheta_d(1,1)),-l2*sin(-Qphase2(1,1)+Qphase2(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))];
          ddtheta_d=Jaco\([0;a]-JacoDot*dtheta_d);% Don't use inv(Jaco), which would increase error. 
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control2(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
        end
        
        % Third part of trajectory - curve 
        if t>11 && t<=12
         a=0.1;
          pv=[-0.1+(1/2)*a*(t-11)^2,-0.25+0.1*(t-11)-(1/2)*a*(t-11)^2];
          Qphase3=InverseKinematics(pv);
          theta_d=Qphase3;
          Jaco=[-l2*sin(-Qphase3(1,1)+Qphase3(2,1))-l1*sin(-Qphase3(1,1)),-l2*sin(-Qphase3(1,1)+Qphase3(2,1));
                 l2*cos(-Qphase3(1,1)+Qphase3(2,1))+l1*cos(-Qphase3(1,1)),l2*cos(-Qphase3(1,1)+Qphase3(2,1))];
          dtheta_d=Jaco\[a*(t-11);0.1-a*(t-11)];
          JacoDot=[-l2*cos(-Qphase3(1,1)+Qphase3(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*cos(-Qphase3(1,1))*(dtheta_d(1,1)),-l2*cos(-Qphase3(1,1)+Qphase3(2,1))*(dtheta_d(1,1)+dtheta_d(2,1));
                   -l2*sin(-Qphase3(1,1)+Qphase3(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*sin(-Qphase3(1,1))*(dtheta_d(1,1)),-l2*sin(-Qphase3(1,1)+Qphase3(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))];
          ddtheta_d=Jaco\([a;-a]-JacoDot*dtheta_d);% Don't use inv(Jaco), which would increase error. 
         theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control3(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
        end
        
        % Fourth part of trajectory - horizontal line
         if t>12 && t<=13
         a=0.1;
          pv=[-0.05+0.1*(t-12),-0.2];
          Qphase4=InverseKinematics(pv);
          theta_d=Qphase4;
          Jaco=[-l2*sin(-Qphase4(1,1)+Qphase4(2,1))-l1*sin(-Qphase4(1,1)),-l2*sin(-Qphase4(1,1)+Qphase4(2,1));
                 l2*cos(-Qphase4(1,1)+Qphase4(2,1))+l1*cos(-Qphase4(1,1)),l2*cos(-Qphase4(1,1)+Qphase4(2,1))];
          dtheta_d=Jaco\[0.1;0];
          JacoDot=[-l2*cos(-Qphase4(1,1)+Qphase4(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*cos(-Qphase4(1,1))*(dtheta_d(1,1)),-l2*cos(-Qphase4(1,1)+Qphase4(2,1))*(dtheta_d(1,1)+dtheta_d(2,1));
                   -l2*sin(-Qphase4(1,1)+Qphase4(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*sin(-Qphase4(1,1))*(dtheta_d(1,1)),-l2*sin(-Qphase4(1,1)+Qphase4(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))];
          ddtheta_d=Jaco\([0;0]-JacoDot*dtheta_d);% Don't use inv(Jaco), which would increase error. 
         theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control4(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
         end
        
         % Fifth part of trajectory - curve
         if t>13 && t<=14
         a=0.1;
          pv=[0.05+0.1*(t-13)-(1/2)*a*(t-13)^2,-0.2-(1/2)*a*(t-13)^2];
          Qphase5=InverseKinematics(pv);
          theta_d=Qphase5;
          Jaco=[-l2*sin(-Qphase5(1,1)+Qphase5(2,1))-l1*sin(-Qphase5(1,1)),-l2*sin(-Qphase5(1,1)+Qphase5(2,1));
                 l2*cos(-Qphase5(1,1)+Qphase5(2,1))+l1*cos(-Qphase5(1,1)),l2*cos(-Qphase5(1,1)+Qphase5(2,1))];
          dtheta_d=Jaco\[0.1-a*(t-13);-a*(t-13)];
          JacoDot=[-l2*cos(-Qphase5(1,1)+Qphase5(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*cos(-Qphase5(1,1))*(dtheta_d(1,1)),-l2*cos(-Qphase5(1,1)+Qphase5(2,1))*(dtheta_d(1,1)+dtheta_d(2,1));
                   -l2*sin(-Qphase5(1,1)+Qphase5(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*sin(-Qphase5(1,1))*(dtheta_d(1,1)),-l2*sin(-Qphase5(1,1)+Qphase5(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))];
          ddtheta_d=Jaco\([-a;-a]-JacoDot*dtheta_d);% Don't use inv(Jaco), which would increase error. 
         theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control5(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
         end
        
        % Sixth part of trajectory - vertical line
         if t>14 && t<=15
         a=0.1;
          pv=[0.1,-0.25-(0.1*(t-14)-(1/2)*a*(t-14)^2)];
          Qphase6=InverseKinematics(pv);
          theta_d=Qphase6;
          Jaco=[-l2*sin(-Qphase6(1,1)+Qphase6(2,1))-l1*sin(-Qphase6(1,1)),-l2*sin(-Qphase6(1,1)+Qphase6(2,1));
                 l2*cos(-Qphase6(1,1)+Qphase6(2,1))+l1*cos(-Qphase6(1,1)),l2*cos(-Qphase6(1,1)+Qphase6(2,1))];
          dtheta_d=Jaco\[0;-0.1+a*(t-14)];
          JacoDot=[-l2*cos(-Qphase6(1,1)+Qphase6(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*cos(-Qphase6(1,1))*(dtheta_d(1,1)),-l2*cos(-Qphase6(1,1)+Qphase6(2,1))*(dtheta_d(1,1)+dtheta_d(2,1));
                   -l2*sin(-Qphase6(1,1)+Qphase6(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))-l1*sin(-Qphase6(1,1))*(dtheta_d(1,1)),-l2*sin(-Qphase6(1,1)+Qphase6(2,1))*(dtheta_d(1,1)+dtheta_d(2,1))];
          ddtheta_d=Jaco\([0;a]-JacoDot*dtheta_d);% Don't use inv(Jaco), which would increase error. 
         theta= x(1:2,1);
        dtheta= x(3:4,1);
        tau=Control6(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
         end
        
        M11=m1*lc1^2+m2*(l1^2+lc2^2+2*l1*lc2*cos(x(2)))+I1+I2;
        M12=m2*(lc2^2+l1*lc2*cos(x(2)))+I2;
        M21=m2*(lc2^2+l1*lc2*cos(x(2)))+I2;
        M22=m2*lc2^2+I2;
        
        C11=-m2*l1*lc2*sin(x(2))*x(4);
        C12=-m2*l1*lc2*sin(x(2))*(x(3)+x(4));
        C21=m2*l1*lc2*sin(x(2))*x(3);
        C22=0;
        
        g1=(m1*lc1+m2*l1)*g*sin(x(1))+m2*lc2*g*sin(x(1)+x(2));
        g2=m2*lc2*g*sin(x(1)+x(2));
        
        M=[M11,M12;
           M21,M22];
        C=[C11,C12;
           C21,C22];
        G=[g1;g2];
        
        invM=inv(M);
        invMC=invM*C; 
        
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = ddtheta_d+tau;
        % ATTENTION - IT IS NOT A PID CONTROLLER. 
        %Based on Computed Torque Contoller, I implemented the control law 
        %and then cancelled the M C and G. This will make the code less 
        %computationally expensive. 
        %dx(3:4) = (-M\C)*x(3:4)+M\tau-M\G is the orignal equation 
        
        
    end

 function tau = Control1(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=100*eye(2);
        Kv=10*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
 end
function tau = Control2(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=350*eye(2);
        Kv=5*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
end
function tau = Control3(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=300*eye(2);
        Kv=5*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
end
function tau = Control4(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=300*eye(2);
        Kv=5*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
end
function tau = Control5(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=300*eye(2);
        Kv=5*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
end
function tau = Control6(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        
        Kp=200*eye(2);
        Kv=5*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e+Kv*de;
        % tau = M*(Kp*e+Kv*de)+C*dtheta+G+M*ddtheta_d is the original
        % control law
 end

end

% Function of inverse kinematics
function [q]=InverseKinematics(p)
global l1 l2
th2=acos((p(:,1)^2+p(:,2)^2-l1^2-l2^2)/(2*l1*l2));
th=acos((l1^2+p(:,1)^2+p(:,2)^2-l2^2)/(2*l1*(p(:,1)^2+p(:,2)^2)^(1/2)));
th1=abs(atan2(p(:,2),p(:,1)))+th;
q=[th1;th2];
end
