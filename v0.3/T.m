function T = T(q1,q2,q3,q4,q5)
%T
%    T = T(Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    02-Oct-2014 17:35:53

t2 = sin(q1);
t3 = cos(q3);
t4 = cos(q1);
t5 = t2.*t3;
t6 = cos(q2);
t7 = sin(q3);
t8 = t4.*t6.*t7;
t9 = sin(q5);
t10 = sin(q4);
t11 = cos(q4);
t12 = sin(q2);
t13 = cos(q5);
t14 = t4.*t7;
t15 = t2.*t3.*t6;
t16 = t14+t15;
t17 = t3.*t4;
t21 = t2.*t6.*t7;
t18 = t17-t21;
t19 = t11.*t16;
t29 = t2.*t10.*t12;
t20 = t19-t29;
t22 = t2.*t7;
t25 = t3.*t4.*t6;
t23 = t22-t25;
t24 = t5+t8;
t26 = t11.*t23;
t27 = t4.*t10.*t12;
t28 = t26+t27;
t30 = t13.*t18;
t31 = t6.*t10;
t32 = t3.*t11.*t12;
t33 = t31+t32;
t34 = t9.*t20;
t35 = -t30+t34;
t36 = t7.*t9.*t12;
t37 = t13.*t20;
t38 = t9.*t18;
t39 = t37+t38;
t58 = t13.*t33;
t40 = t36-t58;
t41 = t10.*t16;
t42 = t2.*t11.*t12;
t43 = t41+t42;
t44 = t6.*t11;
t54 = t3.*t10.*t12;
t45 = t44-t54;
t46 = t9.*t33;
t47 = t7.*t12.*t13;
t48 = t46+t47;
t49 = t10.*t23;
t65 = t4.*t11.*t12;
t50 = t49-t65;
t59 = t13.*t28;
t60 = t9.*t24;
t51 = t59+t60;
t52 = t13.*t24;
t64 = t9.*t28;
t53 = t52-t64;
t55 = t4.*t6.*t10;
t56 = t3.*t4.*t11.*t12;
t57 = t55+t56;
t61 = t2.*t6.*t10;
t62 = t2.*t3.*t11.*t12;
t63 = t61+t62;
t66 = t10.*t12;
t68 = t3.*t6.*t11;
t67 = t66-t68;
t69 = t9.*t67;
t70 = t13.*t16;
t71 = t9.*t11.*t18;
t72 = t7.*t9.*t11.*t12;
t73 = t50.^2;
t74 = t51.^2;
t75 = t52-t64;
t76 = t39.*t40;
t77 = t43.*t45;
t146 = t35.*t48;
t78 = t76+t77-t146;
t79 = t13.*t63;
t151 = t2.*t7.*t9.*t12;
t80 = t79-t151;
t81 = t9.*t63;
t82 = t2.*t7.*t12.*t13;
t83 = t81+t82;
t84 = t53.*t83;
t85 = t2.*t6.*t11;
t153 = t2.*t3.*t10.*t12;
t86 = t85-t153;
t87 = t50.*t86;
t152 = t51.*t80;
t88 = t84+t87-t152;
t89 = t4.*t6.*t11;
t148 = t3.*t4.*t10.*t12;
t90 = t89-t148;
t91 = t45.*t90;
t92 = t9.*t57;
t93 = t4.*t7.*t12.*t13;
t94 = t92+t93;
t95 = t48.*t94;
t96 = t13.*t57;
t149 = t4.*t7.*t9.*t12;
t97 = t96-t149;
t150 = t40.*t97;
t98 = t91+t95-t150;
t99 = t13.*t67;
t100 = t6.*t7.*t9;
t101 = t99+t100;
t102 = t39.*t101;
t154 = t6.*t7.*t13;
t103 = t69-t154;
t104 = t35.*t103;
t105 = t11.*t12;
t106 = t3.*t6.*t10;
t107 = t105+t106;
t155 = t43.*t107;
t108 = t102+t104-t155;
t109 = t3.*t9.*t12;
t110 = t7.*t11.*t12.*t13;
t111 = t109+t110;
t112 = t39.*t111;
t157 = t3.*t12.*t13;
t113 = t35.*(t72-t157);
t114 = t7.*t10.*t12.*t43;
t115 = t112+t113+t114;
t116 = t13.*t23;
t117 = t9.*t11.*t24;
t118 = t116+t117;
t119 = t48.*t118;
t120 = t9.*t23;
t158 = t11.*t13.*t24;
t121 = t120-t158;
t122 = t40.*t121;
t123 = t9.*t16;
t156 = t11.*t13.*t18;
t124 = t123-t156;
t125 = t51.*t124;
t126 = t70+t71;
t127 = t53.*t126;
t128 = t33.*t43;
t129 = t13.*t39.*t45;
t130 = t9.*t35.*t45;
t131 = t128+t129+t130;
t132 = t20.*t50;
t133 = t9.*t43.*t53;
t160 = t13.*t43.*t51;
t134 = t132+t133-t160;
t135 = t28.*t45;
t136 = t9.*t48.*t50;
t159 = t13.*t40.*t50;
t137 = t135+t136-t159;
t138 = t35.*t51;
t139 = t39.*(t52-t64);
t140 = t138+t139;
t141 = t40.*t53;
t161 = t48.*t51;
t142 = t141-t161;
t143 = t39.*t48;
t144 = t35.*t40;
t145 = t143+t144;
t147 = t52-t64;
T = reshape([-t2+t5+t8+t52-t9.*t28,0.0,-t4+t17-t21+t30-t9.*t20,t35.*t78+t48.*(t73+t74+t53.^2),-t43.*t78+t45.*(t73+t74+t75.^2),(t36-t58).*(t73+t74+t147.^2)-t39.*t78,-t9.*t63-t2.*t7.*t12-t2.*t7.*t12.*t13,t69-t6.*t7-t6.*t7.*t13,-t9.*t57-t4.*t7.*t12-t4.*t7.*t12.*t13,-t35.*t98-t48.*t88-t53.*t108,-t45.*t88+t43.*t98-t50.*t108,-t40.*t88+t39.*t98-t51.*t108,t14+t15+t70+t71,t72-t3.*t12-t3.*t12.*t13,-t22+t25-t13.*t23-t9.*t11.*t24,-t53.*t115-t35.*(t119+t122-t10.*t24.*t45)+t48.*(t125+t127-t10.*t18.*t50),-t50.*t115+t43.*(t119+t122-t10.*t24.*t45)+t45.*(t125+t127-t10.*t18.*t50),(t36-t58).*(t125+t127-t10.*t18.*t50)-t51.*t115+t39.*(t119+t122-t10.*t24.*t45),-t9.*t43,-t9.*t45,t9.*t50,t35.*t137-t48.*t134+t131.*(t52-t64),-t45.*t134-t43.*t137+t50.*t131,-t40.*t134-t39.*t137+t51.*t131,t39,t40,-t9.*t24-t13.*t28,t35.*t142+t48.*t140-t53.*t145,-t43.*t142+t45.*t140-t50.*t145,-t39.*t142-t51.*t145+t140.*(t36-t58),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6, 15]);