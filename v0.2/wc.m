function wc = wc(q1,q2,q3,q4,u1,u2,u3,u4)
%WC
%    WC = WC(Q1,Q2,Q3,Q4,U1,U2,U3,U4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    03-Oct-2014 05:35:22

t2 = cos(q1);
t3 = sin(q1);
t4 = sin(q2);
t5 = t2.*t4;
t6 = cos(q2);
t7 = t3.*t6;
t8 = t5+t7;
t9 = t3.*t4;
t11 = t2.*t6;
t10 = t9-t11;
t12 = cos(q3);
t13 = t8.*t12;
t14 = sin(q3);
t16 = t10.*t14;
t15 = t13-t16;
t18 = t10.*t12;
t19 = t8.*t14;
t17 = -t18-t19;
t20 = sin(q4);
t21 = t2.*t20;
t22 = cos(q4);
t23 = t3.*t22;
t24 = t21+t23;
t25 = t3.*t20;
t27 = t2.*t22;
t26 = t25-t27;
wc = reshape([0.0,0.0,0.0,0.0,t2.^2.*u1+t3.^2.*u1,t8.*(t8.*u1+t8.*u2)+(t9-t11).*(t10.*u1+t10.*u2),t15.*(t15.*u1+t15.*u2+t15.*u3)+t17.*(t17.*u1+t17.*u2+t17.*u3),t24.*(t24.*u1+t24.*u4)+(t25-t27).*(t26.*u1+t26.*u4),0.0,0.0,0.0,0.0],[4, 3]);