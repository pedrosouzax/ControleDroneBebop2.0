function [d2q,dq,q,d,B] = droneDynamic(k,t,dt,q,dq,v,tempo_voo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = 0.5;
a= 1;
b = tempo_voo;

dx = 0.05*exp(-(t-3)^4)*sin(1.3*pi*t);
dy = 0.1*exp(-(t-6)^4)*sin(1.3*pi*t);
dz = 0.15*exp(-(t-10)^4)*sin(1.3*pi*t);
%dpsi = 0.20*exp(-(t-(a + (b-a).*rand()))^4)*sin(1.3*pi*t);
%dpsi = 0.15*exp(-(t-15)^4)*sin(1.3*pi*t);
dpsi = 0;

d = [dx;dy;dz;dpsi]; %Vetor disturbio
B = 1/m*eye(4); %Matriz de disturbio

%B = 0;

q(4) = normalize_angle_f(q(4), -pi);
%q(4) = wrapTo2Pi(q(4));

fi =  q(4);%normalize_angle_f(q(4), -pi);

M = [k(1)*cos(fi), -k(3)*sin(fi), 0, 0
     k(1)*sin(fi), k(3)*cos(fi), 0, 0
     0, 0, k(5), 0
     0, 0 , 0, k(7)];

N = [k(2)*cos(fi), -k(4)*sin(fi), 0, 0
      k(2)*sin(fi), k(4)*cos(fi), 0, 0
      0, 0, k(6), 0
      0, 0, 0, k(8)];

Rt = [cos(fi), -sin(fi), 0, 0
      sin(fi), cos(fi), 0, 0
      0, 0, 1, 0
      0, 0, 0, 1];

d2q = M*v - N*(Rt')*dq + (B*d);
dq = d2q*dt + dq;    
q = dq*dt + q;

end

