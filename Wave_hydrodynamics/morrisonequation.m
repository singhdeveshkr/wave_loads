clear all
clc %% Next line has inputs for the problem
D = input ('Enter the diameter of the pile (m)...');
d = input ('Enter the water depth (m)...');
rho = 1025;
g = 9.81;
H = input ('Enter the wave height (m)...');
T = input ('Enter the wave period (sec)...');
Vcs = input ('Enter the surface current velocity (m/s)...');
Vcb = input ('Enter the seabed current velocity (m/s)...');
Vc = 0.5 * (Vcs + Vcb);
fprintf('Depth Averaged current is %4.2f m/s\n ',Vc);
depth_ratio = d/(g*T^2); %% to calculate Tapp (We can put Multiplication factor as 1 if Tapp is not required)
x = Vc/(g*T);
fprintf('d/gT^2 = %5.3f \n',depth_ratio);
fprintf('Vc/gTw = %5.3f \n',x);
r = input ('Enter the multiplication factor from fig 5.2 for Vc/gTw...');
Tapp = r * T;
fprintf('Apparent Wave Period is %5.2f m/s \n',Tapp);
L0 = g/(2*pi) * Tapp^2; %% for calculating the actual wave length (Deep or Intermediate)
Cd = input ('Enter the Drag Coefficient...');
Cm = input ('Enter the Inertia Coefficient...');
q = d/L0;
fprintf('d/L0 = %5.3f \n',q);
r = input('Enter d/L from Wave Table...'); %% If deep water condition is to be used then d/L can be same as d/Lo.
L = d/r;
k = 2*pi/L;
fprintf('Wavelength for the given water depth is %5.3f m \n',L);
fprintf('Wave Number is %5.3f cycles/m \n',k);
sec = input('Enter the number of sections per pile...'); %% Code starts
omega = 2*pi/Tapp;
y = 0:-d/sec:-d;
theta = 0:5:360;
a = length(theta);
b = length(y);
for i = 1:a
   for j = 1:b
        u(i,j) = (((0.5 * H * omega * cosh(k*(d+y(j)))) * (cosd(theta(i))))/(cosh(k*d)))+((Vcs * (d+y(j)))/d);
        acc(i,j) = - (0.5 * omega^2 * H * cosh(k*(d+y(j))) * sind(theta(i)))/(cosh(k*d));
        Fd(i,j) = (0.5 * Cd * rho * D * u(i,j) * abs(u(i,j)))/1000;
        Fi(i,j) = (pi * D^2 * Cm * rho * acc(i,j) * 0.25)/1000;
        F(i,j) = Fd(i,j) + Fi(i,j);
   end
end
for i = 1:a
    Drag(i) = 0.5*((Fd(i,1)+Fd(i,end))+(2*(sum(Fd(i,:))-Fd(i,1)-Fd(i,end))))*(d/sec);
    Inertia(i) = 0.5*((Fi(i,1)+Fi(i,end))+(2*(sum(Fi(i,:))-Fi(i,1)-Fi(i,end))))*(d/sec);
    Ft(i) = 0.5*((F(i,1)+F(i,end))+(2*(sum(F(i,:))-F(i,1)-F(i,end))))*(d/sec);
end
figure
hold on
grid on
grid minor
plot (theta,Drag,"LineWidth",3)
plot (theta,Inertia,"LineWidth",3)
plot (theta,Ft,"LineWidth",3)
xlabel('Phase Angle')
ylabel('Force (Drag,Inertia,Total)')
title ('Force on a pile')
legend('Drag Force','Inertia Force','Total Force')
F_Max = max(Ft);
fprintf ('Maximum force on the pile is %5.2f kN \n',F_Max);