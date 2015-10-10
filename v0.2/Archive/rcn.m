function rcn = rcn(q1,q2,q3,q4,q5,q6,q7,q8,q10,q11,q12)
%RCN
%    RCN = RCN(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q10,Q11,Q12)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    01-Oct-2014 21:43:12

t2 = cos(q1);
t3 = cos(q4);
t4 = cos(q3);
t5 = cos(q2);
t6 = sin(q1);
t7 = sin(q3);
t8 = cos(q5);
t9 = sin(q2);
t10 = sin(q4);
t11 = t9.*t10;
t19 = t3.*t4.*t5;
t12 = t11-t19;
t13 = t6.*t12;
t20 = t2.*t3.*t7;
t14 = t13-t20;
t15 = sin(q5);
t16 = t2.*t4;
t22 = t5.*t6.*t7;
t17 = t16-t22;
t18 = cos(q10);
t21 = t14.*t15;
t23 = t8.*t17;
t24 = t21+t23;
t25 = t18.*t24.*(1.0./2.0);
t26 = sin(q10);
t27 = t8.*t14;
t38 = t15.*t17;
t28 = t27-t38;
t29 = sqrt(3.0);
t30 = t3.*t9;
t31 = t4.*t5.*t10;
t32 = t30+t31;
t33 = t6.*t32;
t34 = t2.*t7.*t10;
t35 = t33+t34;
t36 = t18.*t29.*t35.*(1.0./2.0);
t37 = sin(q12);
t39 = cos(q12);
t40 = t5.*t10;
t41 = t3.*t4.*t9;
t42 = t40+t41;
t43 = t3.*t5;
t52 = t4.*t9.*t10;
t44 = t43-t52;
t45 = t15.*t42;
t46 = t7.*t8.*t9;
t47 = t45+t46;
t48 = t18.*t47.*(1.0./2.0);
t49 = t8.*t42;
t55 = t7.*t9.*t15;
t50 = t49-t55;
t51 = sin(q11);
t53 = t18.*t29.*t44.*(1.0./2.0);
t54 = cos(q11);
t56 = t2.*t12;
t57 = t3.*t6.*t7;
t58 = t56+t57;
t59 = t4.*t6;
t60 = t2.*t5.*t7;
t61 = t59+t60;
t62 = t8.*t61;
t70 = t15.*t58;
t71 = t62-t70;
t63 = t18.*t71.*(1.0./2.0);
t64 = t8.*t58;
t65 = t15.*t61;
t66 = t64+t65;
t67 = t26.*t66;
t68 = t2.*t32;
t72 = t6.*t7.*t10;
t69 = t68-t72;
t73 = cos(q6);
t74 = sin(q6);
t75 = t28.*t74;
t76 = t6.*t32.*(1.0./2.0);
t77 = t2.*t7.*t10.*(1.0./2.0);
t78 = t29.*t35.*t73.*(1.0./2.0);
t79 = sin(q8);
t80 = cos(q8);
t81 = t50.*t74;
t82 = sin(q7);
t83 = t29.*t47.*(1.0./2.0);
t84 = t4.*t9.*t10.*(1.0./2.0);
t85 = t29.*t44.*t73.*(1.0./2.0);
t86 = cos(q7);
t87 = t66.*t74;
t88 = t8.*t61.*(3.0./2.0);
t89 = t2.*t32.*(1.0./2.0);
t90 = t29.*t71.*(1.0./2.0);
t91 = t29.*t69.*t73.*(1.0./2.0);
rcn = reshape([t2-t16+t22-t25-t36-t8.*t17.*(3.0./2.0)-t14.*t15.*(3.0./2.0)+t26.*t28-t29.*t35.*(1.0./2.0)+t39.*(t25+t36-t26.*t28)+t37.*t51.*(t76+t77-t24.*t29.*(1.0./2.0))-t37.*t54.*(t18.*t28+t24.*t26.*(1.0./2.0)+t26.*t29.*t35.*(1.0./2.0))+1.0,t2-t16+t22+t75+t78-t8.*t17.*(3.0./2.0)-t14.*t15.*(3.0./2.0)+t29.*t35.*(1.0./2.0)-t24.*t73.*(1.0./2.0)-t80.*(t75+t78-t24.*t73.*(1.0./2.0))+t79.*t82.*(t76+t77+t24.*t29.*(1.0./2.0))-t79.*t86.*(t24.*t74.*(1.0./2.0)+t28.*t73-t29.*t35.*t74.*(1.0./2.0))+1.0,-t48-t53-t7.*t9-t15.*t42.*(3.0./2.0)-t29.*t44.*(1.0./2.0)+t26.*t50+t39.*(t48+t53-t26.*t50)-t7.*t8.*t9.*(3.0./2.0)-t37.*t51.*(t83+t84-t3.*t5.*(1.0./2.0))-t37.*t54.*(t18.*t50+t26.*t47.*(1.0./2.0)+t26.*t29.*t44.*(1.0./2.0)),t81+t85-t7.*t9-t15.*t42.*(3.0./2.0)+t29.*t44.*(1.0./2.0)-t47.*t73.*(1.0./2.0)-t80.*(t81+t85-t47.*t73.*(1.0./2.0))+t79.*t82.*(t83-t84+t3.*t5.*(1.0./2.0))-t7.*t8.*t9.*(3.0./2.0)-t79.*t86.*(t47.*t74.*(1.0./2.0)+t50.*t73-t29.*t44.*t74.*(1.0./2.0)),-t6+t59+t60+t63+t67+t88-t15.*t58.*(3.0./2.0)-t29.*t69.*(1.0./2.0)-t39.*(t63+t67-t18.*t29.*t69.*(1.0./2.0))+t37.*t51.*(t89+t90-t6.*t7.*t10.*(1.0./2.0))-t18.*t29.*t69.*(1.0./2.0)-t37.*t54.*(t18.*t66-t26.*t71.*(1.0./2.0)+t26.*t29.*t69.*(1.0./2.0)),-t6+t59+t60+t87+t88+t91-t15.*t58.*(3.0./2.0)+t29.*t69.*(1.0./2.0)-t80.*(t87+t91+t71.*t73.*(1.0./2.0))+t73.*(t62-t70).*(1.0./2.0)-t79.*t82.*(-t89+t90+t6.*t7.*t10.*(1.0./2.0))+t79.*t86.*(-t66.*t73+t71.*t74.*(1.0./2.0)+t29.*t69.*t74.*(1.0./2.0))],[2, 3]);
