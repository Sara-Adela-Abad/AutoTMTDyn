function rj = rj(q1,q2,q3)
%RJ
%    RJ = RJ(Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    24-Jul-2014 13:59:43

t2 = cos(q1);
t3 = t2.*2.0;
t4 = q1+q2;
t5 = sin(q1);
t6 = cos(t4);
t7 = t6.*2.0;
t8 = q1+q2+q3;
t9 = sin(t4);
rj = reshape([0.0,t3,t3+t7,t3+t7+cos(t8),0.0,0.0,0.0,0.0,0.0,t5.*-2.0,t5.*-2.0-t9.*2.0,t5.*-2.0-t9.*2.0-sin(t8)],[4, 3]);
