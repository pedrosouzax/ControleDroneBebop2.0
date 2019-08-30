function [ F0 ] = droneNominalModel(k, d2q_d, dq_til, dq_d, q_til, q, ru, t_11, t_12)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% ru = 0.7;
% gama = 0.77;
% t_11 = -0.35;
% t_12 = 1.3;


T11 = t_11 * eye(4);
T12 = t_12 * eye(4);

q(4) = normalize_angle_f(q(4), -pi);

%q(4) = wrapTo2Pi(q(4));

fi = q(4); %normalize_angle_f(q(4), -pi);

%fi = 0.5;

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

T11_inv = inv(T11);
M_inv = inv(M);

% F0 =  M_inv * (d2q_d - (T11_inv * T12 * dq_til)) + (M_inv * N*Rt') * (dq_d - (T11_inv * T12 * q_til'));

F0 =  M_inv * (d2q_d - T11_inv * T12 * dq_til) + (M_inv*N*Rt') * (dq_d - T11_inv * T12 * q_til);


%F0 =  M_inv * (d2q_d - (T11_inv * T12 * dq_til)) + (N*Rt.'/M_inv) * (dq_d - (T11_inv * T12 * q_til'));
%F0 =  M_inv * (d2q_d - (T11_inv * T12 * dq_til)) + (N*Rt'*M) * (dq_d - (T11_inv * T12 * q_til'));

%
end

