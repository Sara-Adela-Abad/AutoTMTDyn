function ref = ref(q1,q2,q3)
%REF
%    REF = REF(Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    03-Oct-2014 05:35:22

t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q1);
t5 = sin(q2);
t6 = cos(q3);
t7 = t2.*t5;
t8 = t3.*t4;
t9 = t7+t8;
t10 = sin(q3);
t11 = t2.*t3;
t12 = t11-t4.*t5;
ref = [t2.*2.0+t2.*t3.*2.0-t4.*t5.*2.0+t6.*t12-t9.*t10,0.0,t4.*-2.0-t2.*t5.*2.0-t3.*t4.*2.0-t6.*t9-t10.*t12];
