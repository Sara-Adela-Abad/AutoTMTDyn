function rj = rj(q1,q2,q3,q4)
%RJ
%    RJ = RJ(Q1,Q2,Q3,Q4)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    02-Jul-2015 14:28:00

t2 = sin(q2);
t3 = cos(q2);
t4 = cos(q4);
t5 = sin(q4);
t6 = sin(q3);
rj = reshape([0.0,-1.0,-t2-1.0,-t2-t2.*t5.*(1.0./2.0)+t3.*t4.*t6.*(1.0./2.0)-1.0,0.0,0.0,0.0,t4.*cos(q3).*(-1.0./2.0),q1,q1,q1-t3,q1-t3-t3.*t5.*(1.0./2.0)-t2.*t4.*t6.*(1.0./2.0),0.0,0.0,0.0,3.0],[4, 4]);
