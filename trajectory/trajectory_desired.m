%% Planejamento de trajetória
%Função para planejar a trajetória para o drone percorrer
%A trajetória 1 é a elipse
%A trajetória 2 é a trajetória de um ponto inial a um ponto final

function [qd,dqd,d2qd,d3qd,to,tf,ax,ay,az,afi] = trajectory_desired(traj_des,q,dq,d2q, qf, dqf, d2qf, t, dt, traj_opt, vel_d, A, omg, to,tf,ax,ay,az,afi)

if(traj_des == 1) % Trajetória de elipse
        xd = A*cos(omg*t);
        yd= A*sin(omg*t);
        zd=(A/2)*sin(omg*t);
        fid = 0;
        dxd= -omg*A*sin(omg*t);
        dyd= omg*A*cos(omg*t);
        dzd= (omg*A/2)*cos(omg*t);
        dfid = 0;
        d2xd= -(omg^2)*A*cos(omg*t);
        d2yd= -(omg^2)*A*sin(omg*t);
        d2zd= -(omg^2)*(A/2)*sin(omg*t);
        d2fid = 0;
        d3xd = 0;
        d3yd = 0;
        d3zd = 0;
        d3fid = 0;

elseif(traj_des == 2)% Trajetória ponto final

             if ((t/dt)==1)
                to = t;
                td = norm(qf-q)/vel_d;
                tf = to+td;
                [ax,ay,az,afi] = trajectory_par_pol_f(dt,to,tf,q,dq,d2q,qf,dqf,d2qf);
            end
    
            if t<=tf
                [xd,dxd,d2xd,d3xd] = traj_pol_f(ax,t,to);
                [yd,dyd,d2yd,d3yd] = traj_pol_f(ay,t,to);
                [zd,dzd,d2zd,d3zd] = traj_pol_f(az,t,to);
                [fid,dfid,d2fid,d3fid] = traj_pol_f(afi,t,to);
            end

else
disp('trajetória inexistente');    
end

qd = [xd; yd; zd; fid];
dqd = [dxd; dyd; dzd; dfid];
d2qd = [d2xd; d2yd; d2zd; d2fid];
d3qd = [d3xd; d3yd; d3zd; d3fid];