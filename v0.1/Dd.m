function Dd = Dd(q1,q2,q3,u1,u2,u3)
%DD
%    DD = DD(Q1,Q2,Q3,U1,U2,U3)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    24-Jul-2014 13:59:43

t2 = cos(q1);
t3 = q1+q2;
t4 = cos(t3);
t5 = sin(q1);
t6 = sin(t3);
t7 = u1+u2;
t8 = q1+q2+q3;
t9 = cos(t8);
t10 = t4.*2.0;
t11 = t2.*2.0;
t12 = t9+t10;
t13 = sin(t8);
t14 = t6.*2.0;
t15 = t5.*2.0;
t16 = t13+t14;
t17 = t16.*u2;
t18 = t13.*u3;
t19 = u1+u2+u3;
Dd = reshape([-t2.*u1,0.0,t5.*u1,-u1.*(t4+t11)-t4.*u2,0.0,u1.*(t6+t15)+t6.*u2,-t9.*u3-t12.*u2-u1.*(t9+t10+t11),0.0,t17+t18+u1.*(t13+t14+t15),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t4.*t7,0.0,t6.*t7,-t9.*u3-t12.*u1-t12.*u2,0.0,t17+t18+t16.*u1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t9.*t19,0.0,t13.*t19,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[18, 3]);