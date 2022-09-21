%% Lecture 11, Question 1
clear; clc;

R_E= 6378; %radius of Earth, in km 
mu= 3.986*10^5; %Gravitational parameter of Earth, in km^3/s^2
r_alt= 720; %Altitude of the S/C, in km
r= R_E+r_alt; %S/C distance from the center of the Earth, in km
n= sqrt(mu/r^3); %Mean motion of first orbit, in rad/s

y_obs= [14.17; 14.54; 16.96; 21.09; 27.80; 36.56];

t0= 0;
t1= 4*60;
t2= 8*60;
t3= 12*60;
t4= 16*60;
t5= 20*60;

Phi_t0= [4-3*cos(n*t0) 0 sin(n*t0)/n 2*(1-cos(n*t0))/n; 6*(sin(n*t0)-n*t0) 1 -2*(1-cos(n*t0))/n ((4*sin(n*t0)/n) - 3*t0); 3*n*sin(n*t0) 0 cos(n*t0) 2*sin(n*t0); -6*n*(1-cos(n*t0)) 0 -2*sin(n*t0) (4*cos(n*t0) - 3)];%State transformation matrix which multipies initial conditions [x; y; z; x_dot; y_dot; z_dot] to give the final positions

Phi_t1= [4-3*cos(n*t1) 0 sin(n*t1)/n 2*(1-cos(n*t1))/n; 6*(sin(n*t1)-n*t1) 1 -2*(1-cos(n*t1))/n ((4*sin(n*t1)/n) - 3*t1); 3*n*sin(n*t1) 0 cos(n*t1) 2*sin(n*t1); -6*n*(1-cos(n*t1)) 0 -2*sin(n*t1) (4*cos(n*t1) - 3)];

Phi_t2= [4-3*cos(n*t2) 0 sin(n*t2)/n 2*(1-cos(n*t2))/n; 6*(sin(n*t2)-n*t2) 1 -2*(1-cos(n*t2))/n ((4*sin(n*t2)/n) - 3*t2); 3*n*sin(n*t2) 0 cos(n*t2) 2*sin(n*t2); -6*n*(1-cos(n*t2)) 0 -2*sin(n*t2) (4*cos(n*t2) - 3)];

Phi_t3= [4-3*cos(n*t3) 0 sin(n*t3)/n 2*(1-cos(n*t3))/n; 6*(sin(n*t3)-n*t3) 1 -2*(1-cos(n*t3))/n ((4*sin(n*t3)/n) - 3*t3); 3*n*sin(n*t3) 0 cos(n*t3) 2*sin(n*t3); -6*n*(1-cos(n*t3)) 0 -2*sin(n*t3) (4*cos(n*t3) - 3)];

Phi_t4= [4-3*cos(n*t4) 0 sin(n*t4)/n 2*(1-cos(n*t4))/n; 6*(sin(n*t4)-n*t4) 1 -2*(1-cos(n*t4))/n ((4*sin(n*t4)/n) - 3*t4); 3*n*sin(n*t4) 0 cos(n*t4) 2*sin(n*t4); -6*n*(1-cos(n*t4)) 0 -2*sin(n*t4) (4*cos(n*t4) - 3)];

Phi_t5= [4-3*cos(n*t5) 0 sin(n*t5)/n 2*(1-cos(n*t5))/n; 6*(sin(n*t5)-n*t5) 1 -2*(1-cos(n*t5))/n ((4*sin(n*t5)/n) - 3*t5); 3*n*sin(n*t5) 0 cos(n*t5) 2*sin(n*t5); -6*n*(1-cos(n*t5)) 0 -2*sin(n*t5) (4*cos(n*t5) - 3)];


%x= [-8.022; 11.68; -0.0026125; 1.67e-4];
x= [-8.022; 11.68; 0.0024125; -1.67e-4];

z=0;
RMS=100;
delta_RMS=10;
%%{
H0= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
H1= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
H2= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
H3= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
H4= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
H5= [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
%}

%{
H0= [x(1)/y_obs(1) x(2)/y_obs(1) 0 0];
H1= [x(1)/y_obs(2) x(2)/y_obs(2) 0 0];
H2= [x(1)/y_obs(3) x(2)/y_obs(3) 0 0];
H3= [x(1)/y_obs(4) x(2)/y_obs(4) 0 0];
H4= [x(1)/y_obs(5) x(2)/y_obs(5) 0 0];
H5= [x(1)/y_obs(6) x(2)/y_obs(6) 0 0];
%}
H=[H0; H1; H2; H3; H4; H5];

A0= H0*Phi_t0;
A1= H1*Phi_t1;
A2= H2*Phi_t2;
A3= H3*Phi_t3;
A4= H4*Phi_t4;
A5= H5*Phi_t5;
A= [A0; A1; A2; A3; A4; A5];

%%{
while (RMS >= 1)
y_pred= A*x;

e= y_obs - y_pred;

P= A'*A; 
delta_x= P\A'*e;

x= x + delta_x; 
%delta_RMS= abs(RMS - sqrt((e' * eye(6) * e) / 6));
RMS= sqrt((e' * eye(6) * e) / 6);

z=z+1;
end

fprintf('Difference in observed range and predicted range within a micron: \n')
RMS

fprintf('Number of times it ran: \n')
z


fprintf('\n \nEstimated x (m): \n')
x
%}