%%%%%%%%%%%%%%%%%% Federal University of Sao Carlos - UFSCar %%%%%%%%%%%%%%
%%%%%% Author: Isabella Cristina Souza Faria.                %%%%%%%%%%%%%%
%%%%%% e-mail: isamoreno2009@gmail.com                       %%%%%%%%%%%%%%
%%%%%% Professor Advisor: Roberto Santos Inoue.              %%%%%%%%%%%%%%
%%%%%% E-mail: rsinoue@ufscar.br                             %%%%%%%%%%%%%%
%%%%%% Date: January 20,2015                                 %%%%%%%%%%%%%%
%%%%%% Reference: Luukonen, Teppo. Modelling and control of quadcopter.
%%%%%% Aalto Universisty School of Science
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Revision 1: Roberto Santos Inoue - January 29, 2015   %%%%%%%%%%%%%%


function [ax,ay,az,afi] = trajectory_par_pol_f(dt,to,tf,so,dso,d2so,sf,dsf,d2sf);


xo= so(1);
yo= so(2);
zo= so(3);
fio= so(4);
xf= sf(1);
yf= sf(2);
zf= sf(3);
fif=sf(4);

dxo= dso(1);
dyo= dso(2);
dzo= dso(3);
dfio= dso(4);
dxf= dsf(1);
dyf= dsf(2);
dzf= dsf(3);
dfif= dsf(4);

d2xo= d2so(1);
d2yo= d2so(2);
d2zo= d2so(3);
d2fio= d2so(4);
d2xf= d2sf(1);
d2yf= d2sf(2);
d2zf= d2sf(3);
d2fif= d2sf(4);


[ax] = par_pol_f(to,tf,xo,dxo,d2xo,xf,dxf,d2xf);
[ay] = par_pol_f(to,tf,yo,dyo,d2yo,yf,dyf,d2yf);
[az] = par_pol_f(to,tf,zo,dzo,d2zo,zf,dzf,d2zf);
[afi] = par_pol_f(to,tf,fio,dfio,d2fio,fif,dfif,d2fif);