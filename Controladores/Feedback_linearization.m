function [ v_fb ] = Feedback_linearization(k,q,d2q_d,dq_d,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%u(1:4)= tanh(u(1:4));

fi = normalize_angle_f(q(4),-pi);
%fi = wrapTo2Pi(q(4));

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

M_inv = inv(M);

v_fb = M_inv*(u + d2q_d - N*Rt'*dq_d);

end

