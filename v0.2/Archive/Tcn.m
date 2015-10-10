function Tcn = Tcn(q1,q2,q3,q4,q5,q6,q7,q8,q10,q11,q12)
%TCN
%    TCN = TCN(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q10,Q11,Q12)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    01-Oct-2014 21:43:38

t2 = cos(q4);
t3 = sin(q1);
t4 = cos(q3);
t5 = cos(q1);
t6 = cos(q2);
t7 = sin(q3);
t8 = cos(q5);
t9 = sin(q2);
t10 = sin(q4);
t11 = t9.*t10;
t21 = t2.*t4.*t6;
t12 = t11-t21;
t13 = t5.*t12;
t14 = t2.*t3.*t7;
t15 = t13+t14;
t16 = sin(q5);
t17 = t3.*t4;
t18 = t5.*t6.*t7;
t19 = t17+t18;
t20 = cos(q10);
t22 = t8.*t19;
t35 = t15.*t16;
t36 = t22-t35;
t23 = t20.*t36.*(1.0./2.0);
t24 = sin(q10);
t25 = t8.*t15;
t26 = t16.*t19;
t27 = t25+t26;
t28 = t24.*t27;
t29 = sqrt(3.0);
t30 = t2.*t9;
t31 = t4.*t6.*t10;
t32 = t30+t31;
t33 = t5.*t32;
t37 = t3.*t7.*t10;
t34 = t33-t37;
t38 = sin(q12);
t39 = cos(q6);
t40 = sin(q6);
t41 = t27.*t40;
t42 = t8.*t19.*(3.0./2.0);
t43 = t5.*t32.*(1.0./2.0);
t44 = t29.*t36.*(1.0./2.0);
t45 = t29.*t34.*t39.*(1.0./2.0);
t46 = sin(q8);
t47 = t6.*t10;
t48 = t2.*t4.*t9;
t49 = t47+t48;
t50 = cos(q12);
t51 = t3.*t16.*t49;
t52 = t3.*t7.*t8.*t9;
t53 = t51+t52;
t54 = t3.*t8.*t49;
t60 = t3.*t7.*t9.*t16;
t55 = t54-t60;
t56 = t24.*t55;
t57 = t2.*t6;
t61 = t4.*t9.*t10;
t58 = t57-t61;
t59 = cos(q11);
t62 = sin(q11);
t63 = t3.*t20.*t29.*t58.*(1.0./2.0);
t64 = cos(q8);
t65 = t40.*t55;
t66 = cos(q7);
t67 = sin(q7);
t68 = t3.*t58.*(1.0./2.0);
t69 = t3.*t29.*t39.*t58.*(1.0./2.0);
t70 = t2.*t4.*t5;
t73 = t2.*t3.*t6.*t7;
t71 = t70-t73;
t72 = t5.*t7;
t74 = t3.*t4.*t6;
t75 = t72+t74;
t76 = t16.*t71;
t77 = t8.*t75;
t78 = t76+t77;
t79 = t8.*t71;
t84 = t16.*t75;
t80 = t79-t84;
t81 = t24.*t80;
t82 = t4.*t5.*t10;
t85 = t3.*t6.*t7.*t10;
t83 = t82-t85;
t86 = t20.*t29.*t83.*(1.0./2.0);
t87 = t16.*t71.*(3.0./2.0);
t88 = t39.*t78.*(1.0./2.0);
t89 = t8.*t75.*(3.0./2.0);
t90 = t29.*t78.*(1.0./2.0);
t91 = t4.*t5.*t10.*(1.0./2.0);
t92 = t29.*t39.*t83.*(1.0./2.0);
t93 = t3.*t32;
t94 = t5.*t7.*t10;
t95 = t93+t94;
t96 = t3.*t12;
t99 = t2.*t5.*t7;
t97 = t96-t99;
t98 = t8.*t24.*t95;
t100 = t20.*t29.*t97.*(1.0./2.0);
t101 = t29.*t97.*(1.0./2.0);
t102 = t16.*t39.*t95.*(1.0./2.0);
t103 = t29.*t39.*t97.*(1.0./2.0);
t104 = t3.*t12.*(1.0./2.0);
t105 = t16.*t29.*t95.*(1.0./2.0);
t106 = t4.*t5;
t109 = t3.*t6.*t7;
t107 = t106-t109;
t114 = t8.*t97;
t115 = t16.*t107;
t108 = t114-t115;
t110 = t16.*t97;
t111 = t8.*t107;
t112 = t110+t111;
t113 = t24.*t112;
t116 = t40.*t112;
t117 = t16.*t107.*(3.0./2.0);
t118 = t39.*t108;
t119 = t40.*t112.*(1.0./2.0);
t128 = t29.*t40.*t95.*(1.0./2.0);
t120 = t118+t119-t128;
t121 = t40.*t108;
t122 = t29.*t39.*t95.*(1.0./2.0);
t199 = t39.*t112.*(1.0./2.0);
t123 = t121+t122-t199;
t124 = t3.*t32.*(1.0./2.0);
t125 = t29.*t112.*(1.0./2.0);
t126 = t5.*t7.*t10.*(1.0./2.0);
t127 = t124+t125+t126;
t129 = t20.*t108;
t130 = t24.*t112.*(1.0./2.0);
t131 = t24.*t29.*t95.*(1.0./2.0);
t132 = t129+t130+t131;
t133 = t20.*t112.*(1.0./2.0);
t134 = t20.*t29.*t95.*(1.0./2.0);
t198 = t24.*t108;
t135 = t133+t134-t198;
t136 = t124-t125+t126;
t137 = t12.*t16;
t143 = t6.*t7.*t8;
t138 = t137-t143;
t139 = t20.*t138.*(1.0./2.0);
t140 = t8.*t12;
t141 = t6.*t7.*t16;
t142 = t140+t141;
t144 = t20.*t29.*t32.*(1.0./2.0);
t145 = t29.*t32.*(1.0./2.0);
t146 = t12.*t16.*(3.0./2.0);
t147 = t39.*t138.*(1.0./2.0);
t148 = t2.*t9.*(1.0./2.0);
t149 = t4.*t6.*t10.*(1.0./2.0);
t150 = t29.*t32.*t39.*(1.0./2.0);
t151 = t4.*t8.*t9;
t156 = t2.*t7.*t9.*t16;
t152 = t151-t156;
t153 = t4.*t9.*t16;
t154 = t2.*t7.*t8.*t9;
t155 = t153+t154;
t157 = t7.*t9.*t10.*t20.*t29.*(1.0./2.0);
t158 = t29.*t152.*(1.0./2.0);
t159 = t2.*t7.*t9.*t16.*(3.0./2.0);
t160 = t8.*t24.*t58;
t161 = t20.*t29.*t49.*(1.0./2.0);
t162 = t29.*t49.*(1.0./2.0);
t163 = t16.*t39.*t58.*(1.0./2.0);
t164 = t6.*t10.*(1.0./2.0);
t165 = t16.*t29.*t58.*(1.0./2.0);
t166 = t2.*t4.*t9.*(1.0./2.0);
t167 = t29.*t39.*t49.*(1.0./2.0);
t168 = t8.*t49;
t175 = t7.*t9.*t16;
t169 = t168-t175;
t170 = t20.*t169.*(1.0./2.0);
t171 = t16.*t49;
t172 = t7.*t8.*t9;
t173 = t171+t172;
t174 = t24.*t173;
t176 = t39.*t169.*(1.0./2.0);
t177 = t40.*t173;
t178 = t7.*t9.*t16.*(3.0./2.0);
t179 = t39.*t169;
t180 = t40.*t173.*(1.0./2.0);
t188 = t29.*t40.*t58.*(1.0./2.0);
t181 = t179+t180-t188;
t182 = t40.*t169;
t183 = t29.*t39.*t58.*(1.0./2.0);
t184 = t182+t183-t39.*t173.*(1.0./2.0);
t185 = t2.*t6.*(1.0./2.0);
t186 = t29.*t173.*(1.0./2.0);
t193 = t4.*t9.*t10.*(1.0./2.0);
t187 = t185+t186-t193;
t189 = t20.*t169;
t190 = t24.*t173.*(1.0./2.0);
t191 = t24.*t29.*t58.*(1.0./2.0);
t192 = t189+t190+t191;
t194 = t20.*t173.*(1.0./2.0);
t195 = t20.*t29.*t58.*(1.0./2.0);
t196 = t194+t195-t24.*t169;
t197 = -t185+t186+t193;
t200 = t16.*t97.*(3.0./2.0);
t201 = t8.*t107.*(3.0./2.0);
t202 = t29.*t95.*(1.0./2.0);
t203 = t5.*t16.*t49;
t204 = t5.*t7.*t8.*t9;
t205 = t203+t204;
t206 = t20.*t205.*(1.0./2.0);
t207 = t5.*t8.*t49;
t209 = t5.*t7.*t9.*t16;
t208 = t207-t209;
t210 = t5.*t20.*t29.*t58.*(1.0./2.0);
t211 = t40.*t208;
t212 = t5.*t58.*(1.0./2.0);
t213 = t5.*t29.*t39.*t58.*(1.0./2.0);
t214 = t2.*t3.*t4;
t215 = t2.*t5.*t6.*t7;
t216 = t214+t215;
t217 = t3.*t7;
t220 = t4.*t5.*t6;
t218 = t217-t220;
t219 = t16.*t216;
t221 = t8.*t218;
t222 = t219+t221;
t223 = t8.*t216;
t229 = t16.*t218;
t224 = t223-t229;
t225 = t24.*t224;
t226 = t3.*t4.*t10;
t227 = t5.*t6.*t7.*t10;
t228 = t226+t227;
t230 = t20.*t29.*t228.*(1.0./2.0);
t231 = t39.*t222.*(1.0./2.0);
t232 = t29.*t228.*(1.0./2.0);
t233 = t29.*t222.*(1.0./2.0);
t234 = t3.*t4.*t10.*(1.0./2.0);
t235 = t5.*t6.*t7.*t10.*(1.0./2.0);
t236 = t29.*t39.*t228.*(1.0./2.0);
t237 = t8.*t24.*t34;
t238 = t15.*t20.*t29.*(1.0./2.0);
t239 = t15.*t29.*(1.0./2.0);
t240 = t16.*t34.*t39.*(1.0./2.0);
t241 = t15.*t29.*t39.*(1.0./2.0);
t242 = t5.*t12.*(1.0./2.0);
t243 = t16.*t29.*t34.*(1.0./2.0);
t244 = t2.*t3.*t7.*(1.0./2.0);
t245 = t36.*t40.*(1.0./2.0);
t246 = t29.*t34.*t40.*(1.0./2.0);
t247 = t39.*(t22-t35).*(1.0./2.0);
t248 = t27.*t39;
t249 = t3.*t7.*t10.*(1.0./2.0);
t250 = -t43+t44+t249;
t251 = t41+t45+t247;
t252 = t20.*t27;
t253 = t24.*t29.*t34.*(1.0./2.0);
t255 = t24.*t36.*(1.0./2.0);
t254 = t252+t253-t255;
Tcn = reshape([-t3+t17+t18+t23+t28+t42-t15.*t16.*(3.0./2.0)-t29.*t34.*(1.0./2.0)-t50.*(t23+t28-t20.*t29.*t34.*(1.0./2.0))+t38.*t62.*(t43+t44-t3.*t7.*t10.*(1.0./2.0))-t20.*t29.*t34.*(1.0./2.0)-t38.*t59.*t254,0.0,-t5+t106-t109+t133+t134-t198+t200+t201+t202-t50.*t135+t38.*t59.*t132-t38.*t62.*t136,t56-t63-t20.*t53.*(1.0./2.0)+t50.*(-t56+t63+t20.*t53.*(1.0./2.0))+t38.*t62.*(t68-t29.*t53.*(1.0./2.0))-t38.*t59.*(t20.*t55+t24.*t53.*(1.0./2.0)+t3.*t24.*t29.*t58.*(1.0./2.0))-t3.*t7.*t9-t3.*t16.*t49.*(3.0./2.0)-t3.*t29.*t58.*(1.0./2.0)-t3.*t7.*t8.*t9.*(3.0./2.0),t139+t144+t145+t146-t6.*t7-t24.*t142-t50.*(t139+t144-t24.*t142)-t6.*t7.*t8.*(3.0./2.0)-t38.*t62.*(t148+t149-t29.*t138.*(1.0./2.0))+t38.*t59.*(t20.*t142+t24.*t138.*(1.0./2.0)+t24.*t29.*t32.*(1.0./2.0)),-t206-t210+t24.*t208+t50.*(t206+t210-t24.*t208)+t38.*t62.*(t212-t29.*t205.*(1.0./2.0))-t38.*t59.*(t20.*t208+t24.*t205.*(1.0./2.0)+t5.*t24.*t29.*t58.*(1.0./2.0))-t5.*t7.*t9-t5.*t16.*t49.*(3.0./2.0)-t5.*t29.*t58.*(1.0./2.0)-t5.*t7.*t8.*t9.*(3.0./2.0),t72+t74-t81-t86+t87+t89+t20.*t78.*(1.0./2.0)-t29.*t83.*(1.0./2.0)+t50.*(t81+t86-t20.*t78.*(1.0./2.0))+t38.*t62.*(t90+t91-t3.*t6.*t7.*t10.*(1.0./2.0))+t38.*t59.*(t20.*t80+t24.*t78.*(1.0./2.0)-t24.*t29.*t83.*(1.0./2.0)),-t157+t159+t50.*(t157+t20.*t152.*(1.0./2.0)+t24.*t155)-t4.*t9-t20.*t152.*(1.0./2.0)-t24.*t155-t4.*t8.*t9.*(3.0./2.0)-t38.*t62.*(t158-t7.*t9.*t10.*(1.0./2.0))-t38.*t59.*(-t20.*t155+t24.*t152.*(1.0./2.0)+t7.*t9.*t10.*t24.*t29.*(1.0./2.0))-t7.*t9.*t10.*t29.*(1.0./2.0),t220+t225+t230+t232-t3.*t7-t8.*t218.*(3.0./2.0)-t16.*t216.*(3.0./2.0)-t20.*t222.*(1.0./2.0)-t50.*(t225+t230-t20.*t222.*(1.0./2.0))-t38.*t62.*(t233+t234+t235)-t38.*t59.*(t20.*t224+t24.*t222.*(1.0./2.0)-t24.*t29.*t228.*(1.0./2.0)),t98+t100+t101-t16.*t95.*(3.0./2.0)-t50.*(t98+t100-t16.*t20.*t95.*(1.0./2.0))-t38.*t62.*(t104+t105-t2.*t5.*t7.*(1.0./2.0))-t16.*t20.*t95.*(1.0./2.0)-t38.*t59.*(t8.*t20.*t95+t16.*t24.*t95.*(1.0./2.0)-t24.*t29.*t97.*(1.0./2.0)),t160+t161+t162-t16.*t58.*(3.0./2.0)-t50.*(t160+t161-t16.*t20.*t58.*(1.0./2.0))-t38.*t62.*(t164+t165+t166)-t16.*t20.*t58.*(1.0./2.0)-t38.*t59.*(t8.*t20.*t58+t16.*t24.*t58.*(1.0./2.0)-t24.*t29.*t49.*(1.0./2.0)),t237+t238+t239-t16.*t34.*(3.0./2.0)-t50.*(t237+t238-t16.*t20.*t34.*(1.0./2.0))-t38.*t62.*(t242+t243+t244)-t16.*t20.*t34.*(1.0./2.0)-t38.*t59.*(t8.*t20.*t34-t15.*t24.*t29.*(1.0./2.0)+t16.*t24.*t34.*(1.0./2.0)),-t113+t117-t8.*t97.*(3.0./2.0)-t20.*t108.*(1.0./2.0)+t50.*(t113+t20.*t108.*(1.0./2.0))+t38.*t59.*(t20.*t112-t24.*t108.*(1.0./2.0))-t29.*t38.*t62.*t108.*(1.0./2.0),-t170-t174+t178+t50.*(t170+t174)-t8.*t49.*(3.0./2.0)+t38.*t59.*(t20.*t173-t24.*t169.*(1.0./2.0))-t29.*t38.*t62.*t169.*(1.0./2.0),t8.*t15.*(-3.0./2.0)-t16.*t19.*(3.0./2.0)-t20.*t27.*(1.0./2.0)+t24.*(t22-t35)+t50.*(t20.*t27.*(1.0./2.0)-t24.*t36)-t38.*t59.*(t24.*t27.*(1.0./2.0)+t20.*t36)-t27.*t29.*t38.*t62.*(1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t129+t130+t131-t50.*t132-t38.*t59.*t135,t189+t190+t191-t50.*t192-t38.*t59.*t196,t252+t253-t255-t50.*t254+t38.*t59.*(t23+t28-t20.*t29.*t34.*(1.0./2.0)),t38.*t62.*t132+t38.*t59.*t136,t38.*t62.*t192-t38.*t59.*t197,t38.*t59.*(t43+t44-t249)+t38.*t62.*t254,-t38.*t135-t50.*t59.*t132+t50.*t62.*t136,-t38.*t196-t50.*t59.*t192-t50.*t62.*t197,t38.*(t23+t28-t20.*t29.*t34.*(1.0./2.0))+t50.*t62.*(t43+t44-t249)-t50.*t59.*t254,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t3+t17+t18+t41+t42+t45+t247-t15.*t16.*(3.0./2.0)+t29.*t34.*(1.0./2.0)-t64.*(t41+t45+t36.*t39.*(1.0./2.0))-t46.*t67.*t250+t46.*t66.*(t245+t246-t27.*t39),0.0,-t5+t106-t109-t121-t122+t199+t200+t201-t202+t64.*t123+t46.*t66.*t120-t46.*t67.*t127,t65+t69-t39.*t53.*(1.0./2.0)-t64.*(t65+t69-t39.*t53.*(1.0./2.0))+t46.*t67.*(t68+t29.*t53.*(1.0./2.0))-t46.*t66.*(t40.*t53.*(1.0./2.0)+t39.*t55-t3.*t29.*t40.*t58.*(1.0./2.0))-t3.*t7.*t9-t3.*t16.*t49.*(3.0./2.0)+t3.*t29.*t58.*(1.0./2.0)-t3.*t7.*t8.*t9.*(3.0./2.0),-t145+t146+t147-t150-t6.*t7-t40.*t142+t64.*(-t147+t150+t40.*t142)-t6.*t7.*t8.*(3.0./2.0)-t46.*t67.*(t148+t149+t29.*t138.*(1.0./2.0))+t46.*t66.*(t40.*t138.*(1.0./2.0)+t39.*t142-t29.*t32.*t40.*(1.0./2.0)),t211+t213-t39.*t205.*(1.0./2.0)-t64.*(t211+t213-t39.*t205.*(1.0./2.0))+t46.*t67.*(t212+t29.*t205.*(1.0./2.0))-t46.*t66.*(t40.*t205.*(1.0./2.0)+t39.*t208-t5.*t29.*t40.*t58.*(1.0./2.0))-t5.*t7.*t9-t5.*t16.*t49.*(3.0./2.0)+t5.*t29.*t58.*(1.0./2.0)-t5.*t7.*t8.*t9.*(3.0./2.0),t72+t74+t87+t88+t89+t92+t29.*t83.*(1.0./2.0)-t40.*t80-t64.*(t88+t92-t40.*t80)-t46.*t67.*(t90-t91+t3.*t6.*t7.*t10.*(1.0./2.0))+t46.*t66.*(t40.*t78.*(1.0./2.0)+t39.*t80+t29.*t40.*t83.*(1.0./2.0)),t159-t4.*t9-t39.*t152.*(1.0./2.0)-t40.*t155+t64.*(t39.*t152.*(1.0./2.0)+t40.*t155-t7.*t9.*t10.*t29.*t39.*(1.0./2.0))-t4.*t8.*t9.*(3.0./2.0)+t46.*t67.*(t158+t7.*t9.*t10.*(1.0./2.0))+t46.*t66.*(t40.*t152.*(-1.0./2.0)+t39.*t155+t7.*t9.*t10.*t29.*t40.*(1.0./2.0))+t7.*t9.*t10.*t29.*(1.0./2.0)+t7.*t9.*t10.*t29.*t39.*(1.0./2.0),-t217+t220-t231-t232-t236-t8.*t218.*(3.0./2.0)-t16.*t216.*(3.0./2.0)+t40.*t224+t64.*(t231+t236-t40.*t224)-t46.*t67.*(-t233+t234+t235)-t46.*t66.*(t40.*t222.*(1.0./2.0)+t39.*t224+t29.*t40.*t228.*(1.0./2.0)),-t101-t102-t103-t16.*t95.*(3.0./2.0)+t64.*(t102+t103-t8.*t40.*t95)+t46.*t67.*(-t104+t105+t2.*t5.*t7.*(1.0./2.0))+t8.*t40.*t95-t46.*t66.*(t8.*t39.*t95+t16.*t40.*t95.*(1.0./2.0)+t29.*t40.*t97.*(1.0./2.0)),-t162-t163-t167-t16.*t58.*(3.0./2.0)+t64.*(t163+t167-t8.*t40.*t58)-t46.*t67.*(t164-t165+t166)+t8.*t40.*t58-t46.*t66.*(t8.*t39.*t58+t16.*t40.*t58.*(1.0./2.0)+t29.*t40.*t49.*(1.0./2.0)),-t239-t240-t241-t16.*t34.*(3.0./2.0)+t64.*(t240+t241-t8.*t34.*t40)-t46.*t67.*(t242-t243+t244)+t8.*t34.*t40-t46.*t66.*(t8.*t34.*t39+t15.*t29.*t40.*(1.0./2.0)+t16.*t34.*t40.*(1.0./2.0)),-t116+t117-t8.*t97.*(3.0./2.0)-t39.*t108.*(1.0./2.0)+t64.*(t116+t39.*t108.*(1.0./2.0))-t46.*t66.*(t40.*t108.*(1.0./2.0)-t39.*t112)+t29.*t46.*t67.*t108.*(1.0./2.0),-t176-t177+t178+t64.*(t176+t177)-t8.*t49.*(3.0./2.0)-t46.*t66.*(t40.*t169.*(1.0./2.0)-t39.*t173)+t29.*t46.*t67.*t169.*(1.0./2.0),t8.*t15.*(-3.0./2.0)-t16.*t19.*(3.0./2.0)-t27.*t39.*(1.0./2.0)+t40.*(t22-t35)+t64.*(t27.*t39.*(1.0./2.0)-t36.*t40)-t46.*t66.*(t27.*t40.*(1.0./2.0)+t36.*t39)+t27.*t29.*t46.*t67.*(1.0./2.0),t118+t119-t64.*t120-t29.*t40.*t95.*(1.0./2.0)+t46.*t66.*t123,t179+t180-t64.*t181-t29.*t40.*t58.*(1.0./2.0)+t46.*t66.*t184,-t245-t246+t248+t64.*(t245+t246-t27.*t39)+t46.*t66.*t251,t46.*t67.*t120+t46.*t66.*t127,t46.*t67.*t181+t46.*t66.*t187,-t46.*t67.*(t245+t246-t248)-t46.*t66.*t250,t46.*t123-t64.*t66.*t120+t64.*t67.*t127,t46.*t184-t64.*t66.*t181+t64.*t67.*t187,t46.*t251+t64.*t66.*(t245+t246-t248)-t64.*t67.*t250,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3, 15, 2]);