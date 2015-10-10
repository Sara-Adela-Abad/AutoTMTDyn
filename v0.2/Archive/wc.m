function wc = wc(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13)
%WC
%    WC = WC(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12,U13)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    01-Oct-2014 22:00:28

t2 = cos(q1);
t3 = sin(q1);
t4 = cos(q2);
t5 = cos(q3);
t6 = sin(q2);
t7 = sin(q4);
t8 = cos(q4);
t9 = sin(q3);
t10 = t4.*t7;
t11 = t5.*t6.*t8;
t12 = t10+t11;
t13 = t6.*t7;
t23 = t4.*t5.*t8;
t14 = t13-t23;
t15 = t6.*t8;
t16 = t4.*t5.*t7;
t17 = t15+t16;
t18 = t4.*t8;
t26 = t5.*t6.*t7;
t19 = t18-t26;
t20 = t2.*t5;
t22 = t3.*t4.*t9;
t21 = t20-t22;
t24 = t3.*t14;
t38 = t2.*t8.*t9;
t25 = t24-t38;
t27 = t3.*t17;
t28 = t2.*t7.*t9;
t29 = t27+t28;
t30 = t3.*t5;
t31 = t2.*t4.*t9;
t32 = t30+t31;
t33 = t2.*t14;
t34 = t3.*t8.*t9;
t35 = t33+t34;
t36 = t2.*t17;
t39 = t3.*t7.*t9;
t37 = t36-t39;
t40 = t32.*u1;
t41 = t2.*t9;
t42 = t3.*t4.*t5;
t43 = t41+t42;
t44 = t43.*u3;
t96 = t3.*t6.*t9.*u2;
t45 = t40+t44-t96;
t46 = t32.*t45;
t47 = t35.*u1;
t48 = t29.*u4;
t49 = t2.*t5.*t8;
t97 = t3.*t4.*t8.*t9;
t50 = t49-t97;
t51 = t3.*t12.*u2;
t98 = t50.*u3;
t52 = t47+t48+t51-t98;
t53 = t35.*t52;
t54 = t2.*t5.*t7;
t99 = t3.*t4.*t7.*t9;
t55 = t54-t99;
t56 = t55.*u3;
t57 = t37.*u1;
t58 = t3.*t19.*u2;
t100 = t25.*u4;
t59 = t56+t57+t58-t100;
t60 = t37.*t59;
t61 = t46+t53+t60;
t62 = t4.*t9.*u2;
t63 = t5.*t6.*u3;
t64 = t62+t63;
t65 = t14.*u2;
t66 = t6.*t8.*t9.*u3;
t102 = t19.*u4;
t67 = t65+t66-t102;
t68 = t25.*t67;
t69 = t17.*u2;
t70 = t12.*u4;
t103 = t6.*t7.*t9.*u3;
t71 = t69+t70-t103;
t72 = t29.*t71;
t101 = t21.*t64;
t73 = t68+t72-t101;
t74 = t37.*u4;
t75 = t3.*t5.*t8;
t76 = t2.*t4.*t8.*t9;
t77 = t75+t76;
t78 = t77.*u3;
t79 = t2.*t12.*u2;
t104 = t25.*u1;
t80 = t74+t78+t79-t104;
t81 = t12.*t80;
t82 = t3.*t5.*t7;
t83 = t2.*t4.*t7.*t9;
t84 = t82+t83;
t85 = t84.*u3;
t86 = t35.*u4;
t87 = t29.*u1;
t105 = t2.*t19.*u2;
t88 = t85+t86+t87-t105;
t89 = t3.*t9;
t108 = t2.*t4.*t5;
t90 = t89-t108;
t91 = t90.*u3;
t92 = t2.*t6.*t9.*u2;
t107 = t21.*u1;
t93 = t91+t92-t107;
t94 = t6.*t9.*t93;
t106 = t19.*t88;
t95 = t81+t94-t106;
t109 = cos(q5);
t110 = sin(q5);
t111 = t25.*t110;
t112 = t21.*t109;
t113 = t111+t112;
t114 = t12.*t110;
t115 = t6.*t9.*t109;
t116 = t114+t115;
t117 = t35.*t109;
t118 = t32.*t110;
t119 = t117+t118;
t120 = t32.*t109;
t121 = t35.*t110;
t122 = -t120+t121;
t123 = t21.*t110;
t125 = t25.*t109;
t124 = t123-t125;
t126 = t12.*t109;
t165 = t6.*t9.*t110;
t127 = t126-t165;
t128 = t119.*u1;
t129 = t43.*t110;
t194 = t50.*t109;
t217 = t129-t194;
t130 = t217.*u3;
t131 = t3.*t12.*t109;
t195 = t3.*t6.*t9.*t110;
t132 = t131-t195;
t133 = t132.*u2;
t134 = t29.*t109.*u4;
t193 = t113.*u5;
t135 = t128+t130+t133+t134-t193;
t136 = t119.*t135;
t137 = t122.*u1;
t138 = t50.*t110;
t139 = t43.*t109;
t140 = t138+t139;
t141 = t3.*t12.*t110;
t142 = t3.*t6.*t9.*t109;
t143 = t141+t142;
t144 = t143.*u2;
t145 = t29.*t110.*u4;
t196 = t140.*u3;
t146 = t137+t144+t145-t196-t124.*u5;
t147 = t122.*t146;
t148 = t60+t136+t147;
t149 = t5.*t6.*t110;
t150 = t6.*t8.*t9.*t109;
t151 = t149+t150;
t152 = t151.*u3;
t153 = t14.*t109;
t154 = t4.*t9.*t110;
t155 = t153+t154;
t156 = t155.*u2;
t157 = t116.*u5;
t197 = t19.*t109.*u4;
t158 = t152+t156+t157-t197;
t159 = t124.*t158;
t160 = t5.*t6.*t109;
t198 = t6.*t8.*t9.*t110;
t161 = t160-t198;
t162 = t161.*u3;
t163 = t14.*t110;
t199 = t4.*t9.*t109;
t164 = t163-t199;
t166 = t127.*u5;
t167 = t19.*t110.*u4;
t200 = t164.*u2;
t168 = t162+t166+t167-t200;
t169 = t113.*t168;
t170 = t77.*t109;
t190 = t90.*t110;
t171 = t170-t190;
t172 = t171.*u3;
t173 = t2.*t12.*t109;
t191 = t2.*t6.*t9.*t110;
t174 = t173-t191;
t175 = t174.*u2;
t176 = t37.*t109.*u4;
t177 = t119.*u5;
t178 = t77.*t110;
t179 = t90.*t109;
t180 = t178+t179;
t181 = t180.*u3;
t182 = t2.*t12.*t110;
t183 = t2.*t6.*t9.*t109;
t184 = t182+t183;
t185 = t184.*u2;
t186 = t37.*t110.*u4;
t192 = t113.*u1;
t187 = t177+t181+t185+t186-t192;
t188 = t116.*t187;
t189 = -t123+t125;
t201 = cos(q6);
t202 = sin(q6);
t203 = sqrt(3.0);
t204 = t189.*t202;
t205 = t29.*t201.*t203.*(1.0./2.0);
t218 = t113.*t201.*(1.0./2.0);
t206 = t204+t205-t218;
t207 = t119.*t202;
t208 = t37.*t201.*t203.*(1.0./2.0);
t210 = t122.*t201.*(1.0./2.0);
t209 = t207+t208-t210;
t211 = t189.*t201;
t212 = t113.*t202.*(1.0./2.0);
t226 = t29.*t202.*t203.*(1.0./2.0);
t213 = t211+t212-t226;
t214 = t119.*t201;
t215 = t122.*t202.*(1.0./2.0);
t219 = t37.*t202.*t203.*(1.0./2.0);
t216 = t214+t215-t219;
t220 = t2.*t17.*(1.0./2.0);
t221 = t122.*t203.*(1.0./2.0);
t286 = t3.*t7.*t9.*(1.0./2.0);
t222 = t220+t221-t286;
t223 = t127.*t202;
t224 = t19.*t201.*t203.*(1.0./2.0);
t227 = t116.*t201.*(1.0./2.0);
t225 = t223+t224-t227;
t228 = t3.*t17.*(1.0./2.0);
t229 = t113.*t203.*(1.0./2.0);
t230 = t2.*t7.*t9.*(1.0./2.0);
t231 = t228+t229+t230;
t232 = t127.*t201;
t233 = t116.*t202.*(1.0./2.0);
t282 = t19.*t202.*t203.*(1.0./2.0);
t234 = t232+t233-t282;
t235 = t12.*t202.*t203.*(1.0./2.0);
t236 = t19.*t109.*t201;
t237 = t19.*t110.*t202.*(1.0./2.0);
t238 = t235+t236+t237;
t239 = t151.*t201;
t240 = t6.*t7.*t9.*t202.*t203.*(1.0./2.0);
t421 = t161.*t202.*(1.0./2.0);
t241 = t239+t240-t421;
t242 = t241.*u3;
t243 = t155.*t201;
t244 = t164.*t202.*(1.0./2.0);
t422 = t17.*t202.*t203.*(1.0./2.0);
t245 = t243+t244-t422;
t246 = t245.*u2;
t247 = t225.*u6;
t248 = t116.*t201;
t423 = t127.*t202.*(1.0./2.0);
t249 = t248-t423;
t250 = t249.*u5;
t420 = t238.*u4;
t251 = t242+t246+t247+t250-t420;
t252 = t213.*t251;
t253 = t161.*t203.*(1.0./2.0);
t254 = t6.*t7.*t9.*(1.0./2.0);
t255 = t253+t254;
t256 = t255.*u3;
t257 = t4.*t7.*(1.0./2.0);
t258 = t5.*t6.*t8.*(1.0./2.0);
t424 = t19.*t110.*t203.*(1.0./2.0);
t259 = t257+t258-t424;
t260 = t6.*t8.*(1.0./2.0);
t261 = t164.*t203.*(1.0./2.0);
t262 = t4.*t5.*t7.*(1.0./2.0);
t263 = t260+t261+t262;
t264 = t127.*t203.*u5.*(1.0./2.0);
t425 = t259.*u4;
t426 = t263.*u2;
t265 = t256+t264-t425-t426;
t266 = t19.*t110.*t201.*(1.0./2.0);
t267 = t12.*t201.*t203.*(1.0./2.0);
t428 = t19.*t109.*t202;
t268 = t266+t267-t428;
t269 = t268.*u4;
t270 = t161.*t201.*(1.0./2.0);
t271 = t151.*t202;
t429 = t6.*t7.*t9.*t201.*t203.*(1.0./2.0);
t272 = t270+t271-t429;
t273 = t272.*u3;
t274 = t155.*t202;
t275 = t17.*t201.*t203.*(1.0./2.0);
t430 = t164.*t201.*(1.0./2.0);
t276 = t274+t275-t430;
t277 = t276.*u2;
t278 = t127.*t201.*(1.0./2.0);
t279 = t116.*t202;
t280 = t278+t279;
t281 = t280.*u5;
t431 = t234.*u6;
t283 = t269+t273+t277+t281-t431;
t284 = t206.*t283;
t427 = t231.*t265;
t285 = t252+t284-t427;
t287 = t119.*t201.*(1.0./2.0);
t288 = t122.*t202;
t289 = t287+t288;
t290 = t289.*u5;
t291 = t206.*u1;
t292 = t174.*t202;
t293 = t2.*t19.*t201.*t203.*(1.0./2.0);
t393 = t184.*t201.*(1.0./2.0);
t294 = t292+t293-t393;
t295 = t37.*t110.*t201.*(1.0./2.0);
t296 = t35.*t201.*t203.*(1.0./2.0);
t395 = t37.*t109.*t202;
t297 = t295+t296-t395;
t298 = t297.*u4;
t299 = t180.*t201.*(1.0./2.0);
t300 = t84.*t201.*t203.*(1.0./2.0);
t396 = t171.*t202;
t301 = t299+t300-t396;
t302 = t301.*u3;
t392 = t216.*u6;
t394 = t294.*u2;
t303 = t290+t291+t298+t302-t392-t394;
t304 = t225.*t303;
t305 = t4.*t8.*(1.0./2.0);
t306 = t116.*t203.*(1.0./2.0);
t340 = t5.*t6.*t7.*(1.0./2.0);
t307 = t305+t306-t340;
t308 = t3.*t5.*t7.*(1.0./2.0);
t309 = t2.*t4.*t7.*t9.*(1.0./2.0);
t397 = t180.*t203.*(1.0./2.0);
t310 = t308+t309-t397;
t311 = t310.*u3;
t312 = t184.*t203.*(1.0./2.0);
t313 = t2.*t19.*(1.0./2.0);
t314 = t312+t313;
t315 = t2.*t14.*(1.0./2.0);
t316 = t3.*t8.*t9.*(1.0./2.0);
t399 = t37.*t110.*t203.*(1.0./2.0);
t317 = t315+t316-t399;
t318 = t317.*u4;
t319 = t231.*u1;
t398 = t314.*u2;
t400 = t119.*t203.*u5.*(1.0./2.0);
t320 = t311+t318+t319-t398-t400;
t321 = t307.*t320;
t322 = t37.*t109.*t201;
t323 = t37.*t110.*t202.*(1.0./2.0);
t324 = t35.*t202.*t203.*(1.0./2.0);
t325 = t322+t323+t324;
t326 = t325.*u4;
t327 = t119.*t202.*(1.0./2.0);
t401 = t122.*t201;
t475 = t327-t401;
t328 = t475.*u5;
t329 = t171.*t201;
t330 = t180.*t202.*(1.0./2.0);
t331 = t84.*t202.*t203.*(1.0./2.0);
t332 = t329+t330+t331;
t333 = t332.*u3;
t334 = t174.*t201;
t335 = t184.*t202.*(1.0./2.0);
t404 = t2.*t19.*t202.*t203.*(1.0./2.0);
t336 = t334+t335-t404;
t337 = t336.*u2;
t402 = t209.*u6;
t403 = t213.*u1;
t338 = t326+t328+t333+t337-t402-t403;
t405 = t234.*t338;
t339 = t304+t321-t405;
t341 = t189.*t201.*(1.0./2.0);
t342 = t113.*t202;
t343 = t341+t342;
t344 = t132.*t202;
t345 = t3.*t19.*t201.*t203.*(1.0./2.0);
t407 = t143.*t201.*(1.0./2.0);
t346 = t344+t345-t407;
t347 = t346.*u2;
t348 = t209.*u1;
t349 = t213.*u6;
t350 = t29.*t110.*t201.*(1.0./2.0);
t351 = t25.*t201.*t203.*(1.0./2.0);
t408 = t29.*t109.*t202;
t352 = t350+t351-t408;
t353 = t140.*t201.*(1.0./2.0);
t354 = t202.*(t129-t194);
t355 = t55.*t201.*t203.*(1.0./2.0);
t356 = t353+t354+t355;
t357 = t356.*u3;
t406 = t343.*u5;
t409 = t352.*u4;
t358 = t347+t348+t349+t357-t406-t409;
t359 = t209.*t358;
t360 = t29.*t109.*t201;
t361 = t29.*t110.*t202.*(1.0./2.0);
t362 = t25.*t202.*t203.*(1.0./2.0);
t363 = t360+t361+t362;
t364 = t113.*t201;
t411 = t189.*t202.*(1.0./2.0);
t365 = t364-t411;
t366 = t365.*u5;
t367 = t140.*t202.*(1.0./2.0);
t368 = t55.*t202.*t203.*(1.0./2.0);
t412 = t201.*t217;
t369 = t367+t368-t412;
t370 = t369.*u3;
t371 = t132.*t201;
t372 = t143.*t202.*(1.0./2.0);
t413 = t3.*t19.*t202.*t203.*(1.0./2.0);
t373 = t371+t372-t413;
t374 = t206.*u6;
t410 = t363.*u4;
t414 = t373.*u2;
t415 = t216.*u1;
t375 = t366+t370+t374-t410-t414-t415;
t376 = t140.*t203.*(1.0./2.0);
t377 = t3.*t4.*t7.*t9.*(1.0./2.0);
t417 = t2.*t5.*t7.*(1.0./2.0);
t378 = t376+t377-t417;
t379 = t143.*t203.*(1.0./2.0);
t380 = t3.*t19.*(1.0./2.0);
t381 = t379+t380;
t382 = t381.*u2;
t383 = t29.*t110.*t203.*(1.0./2.0);
t384 = t2.*t8.*t9.*(1.0./2.0);
t419 = t3.*t14.*(1.0./2.0);
t385 = t383+t384-t419;
t386 = t385.*u4;
t387 = t222.*u1;
t388 = t189.*t203.*u5.*(1.0./2.0);
t418 = t378.*u3;
t389 = t382+t386+t387+t388-t418;
t390 = t222.*t389;
t416 = t216.*t375;
t391 = t359+t390-t416;
t432 = sin(q8);
t433 = cos(q8);
t434 = sin(q7);
t435 = cos(q7);
t436 = sin(q9);
t437 = cos(q9);
t438 = t435.*t436;
t439 = t433.*t434.*t437;
t440 = t438+t439;
t441 = t434.*t436;
t443 = t433.*t435.*t437;
t442 = t441-t443;
t444 = t434.*t437;
t445 = t433.*t435.*t436;
t446 = t444+t445;
t447 = t435.*t437;
t449 = t433.*t434.*t436;
t448 = t447-t449;
t450 = t213.*t442;
t451 = t231.*t440;
t452 = t206.*t432.*t437;
t453 = t450+t451+t452;
t454 = t216.*t446;
t455 = t222.*t448;
t484 = t209.*t432.*t436;
t456 = t454+t455-t484;
t457 = t216.*t442;
t458 = t222.*t440;
t459 = t209.*t432.*t437;
t460 = t457+t458+t459;
t461 = t213.*t446;
t462 = t231.*t448;
t479 = t206.*t432.*t436;
t463 = t461+t462-t479;
t464 = t209.*t433;
t465 = t216.*t432.*t435;
t467 = t222.*t432.*t434;
t466 = t464+t465-t467;
t468 = t206.*t433;
t469 = t213.*t432.*t435;
t474 = t231.*t432.*t434;
t470 = t468+t469-t474;
t471 = t225.*t433;
t472 = t234.*t432.*t435;
t653 = t307.*t432.*t434;
t473 = t471+t472-t653;
t476 = t234.*t446;
t477 = t307.*t448;
t526 = t225.*t432.*t436;
t478 = t476+t477-t526;
t480 = t234.*t442;
t481 = t307.*t440;
t482 = t225.*t432.*t437;
t483 = t480+t481+t482;
t485 = t234.*t433;
t735 = t225.*t432.*t435;
t486 = t485-t735;
t487 = t276.*t433;
t488 = t245.*t432.*t435;
t737 = t263.*t432.*t434;
t489 = t487+t488-t737;
t490 = t489.*u2;
t491 = t225.*t432;
t492 = t307.*t433.*t434;
t738 = t234.*t433.*t435;
t493 = t491+t492-t738;
t494 = t493.*u8;
t495 = t234.*t432.*t434;
t496 = t307.*t432.*t435;
t497 = t495+t496;
t498 = t497.*u7;
t499 = t238.*t432.*t435;
t500 = t259.*t432.*t434;
t739 = t268.*t433;
t501 = t499+t500-t739;
t502 = t280.*t433;
t503 = t249.*t432.*t435;
t504 = t127.*t203.*t432.*t434.*(1.0./2.0);
t505 = t502+t503+t504;
t506 = t505.*u5;
t507 = t272.*t433;
t508 = t255.*t432.*t434;
t509 = t241.*t432.*t435;
t510 = t507+t508+t509;
t511 = t510.*u3;
t736 = t486.*u6;
t740 = t501.*u4;
t512 = t490+t494+t498+t506+t511-t736-t740;
t513 = t470.*t512;
t514 = t225.*t433.*t437;
t515 = t234.*t432.*t435.*t437;
t741 = t307.*t432.*t434.*t437;
t516 = t514+t515-t741;
t517 = t241.*t442;
t518 = t272.*t432.*t437;
t743 = t255.*t440;
t519 = t517+t518-t743;
t520 = t519.*u3;
t521 = t245.*t442;
t522 = t263.*t440;
t523 = t276.*t432.*t437;
t524 = t521+t522+t523;
t525 = t524.*u2;
t527 = t234.*t440;
t745 = t307.*t442;
t528 = t527-t745;
t529 = t225.*t442;
t747 = t234.*t432.*t437;
t530 = t529-t747;
t531 = t530.*u6;
t532 = t259.*t440;
t533 = t268.*t432.*t437;
t748 = t238.*t442;
t534 = t532+t533-t748;
t535 = t534.*u4;
t536 = t249.*t442;
t537 = t280.*t432.*t437;
t749 = t127.*t203.*t440.*(1.0./2.0);
t538 = t536+t537-t749;
t539 = t538.*u5;
t742 = t516.*u8;
t744 = t478.*u9;
t746 = t528.*u7;
t540 = t520+t525+t531+t535+t539-t742-t744-t746;
t541 = t453.*t540;
t542 = t238.*t446;
t543 = t268.*t432.*t436;
t750 = t259.*t448;
t544 = t542+t543-t750;
t545 = t280.*t432.*t436;
t546 = t127.*t203.*t448.*(1.0./2.0);
t752 = t249.*t446;
t547 = t545+t546-t752;
t548 = t255.*t448;
t549 = t272.*t432.*t436;
t754 = t241.*t446;
t550 = t548+t549-t754;
t551 = t225.*t433.*t436;
t552 = t234.*t432.*t435.*t436;
t756 = t307.*t432.*t434.*t436;
t553 = t551+t552-t756;
t554 = t553.*u8;
t555 = t245.*t446;
t556 = t263.*t448;
t757 = t276.*t432.*t436;
t557 = t555+t556-t757;
t558 = t557.*u2;
t559 = t483.*u9;
t560 = t234.*t448;
t758 = t307.*t446;
t561 = t560-t758;
t562 = t225.*t446;
t563 = t234.*t432.*t436;
t564 = t562+t563;
t565 = t564.*u6;
t751 = t544.*u4;
t753 = t547.*u5;
t755 = t550.*u3;
t759 = t561.*u7;
t566 = t463.*(t554+t558+t559+t565-t751-t753-t755-t759);
t567 = t513+t541+t566;
t568 = t373.*t446;
t569 = t381.*t448;
t760 = t346.*t432.*t436;
t570 = t568+t569-t760;
t571 = t570.*u2;
t572 = t378.*t448;
t573 = t369.*t446;
t574 = t356.*t432.*t436;
t575 = t572+t573+t574;
t576 = t456.*u1;
t577 = t213.*t448;
t763 = t231.*t446;
t578 = t577-t763;
t579 = t578.*u7;
t580 = t206.*t446;
t581 = t213.*t432.*t436;
t582 = t580+t581;
t583 = t385.*t448;
t584 = t363.*t446;
t585 = t352.*t432.*t436;
t586 = t583+t584+t585;
t587 = t586.*u4;
t588 = t189.*t203.*t448.*(1.0./2.0);
t589 = t343.*t432.*t436;
t765 = t365.*t446;
t590 = t588+t589-t765;
t591 = t590.*u5;
t592 = t206.*t433.*t436;
t593 = t213.*t432.*t435.*t436;
t766 = t231.*t432.*t434.*t436;
t594 = t592+t593-t766;
t761 = t453.*u9;
t762 = t575.*u3;
t764 = t582.*u6;
t767 = t594.*u8;
t595 = t571+t576+t579+t587+t591-t761-t762-t764-t767;
t596 = t456.*t595;
t597 = t373.*t442;
t598 = t381.*t440;
t599 = t346.*t432.*t437;
t600 = t597+t598+t599;
t601 = t600.*u2;
t602 = t378.*t440;
t603 = t369.*t442;
t768 = t356.*t432.*t437;
t604 = t602+t603-t768;
t605 = t460.*u1;
t606 = t463.*u9;
t607 = t213.*t440;
t770 = t231.*t442;
t608 = t607-t770;
t609 = t608.*u7;
t610 = t206.*t442;
t771 = t213.*t432.*t437;
t611 = t610-t771;
t612 = t385.*t440;
t613 = t363.*t442;
t773 = t352.*t432.*t437;
t614 = t612+t613-t773;
t615 = t614.*u4;
t616 = t206.*t433.*t437;
t617 = t213.*t432.*t435.*t437;
t774 = t231.*t432.*t434.*t437;
t618 = t616+t617-t774;
t619 = t618.*u8;
t620 = t365.*t442;
t621 = t343.*t432.*t437;
t775 = t189.*t203.*t440.*(1.0./2.0);
t622 = t620+t621-t775;
t769 = t604.*u3;
t772 = t611.*u6;
t776 = t622.*u5;
t623 = t601+t605+t606+t609+t615+t619-t769-t772-t776;
t624 = t460.*t623;
t625 = t213.*t433;
t777 = t206.*t432.*t435;
t626 = t625-t777;
t627 = t626.*u6;
t628 = t206.*t432;
t629 = t231.*t433.*t434;
t778 = t213.*t433.*t435;
t630 = t628+t629-t778;
t631 = t356.*t433;
t632 = t378.*t432.*t434;
t780 = t369.*t432.*t435;
t633 = t631+t632-t780;
t634 = t633.*u3;
t635 = t213.*t432.*t434;
t636 = t231.*t432.*t435;
t637 = t635+t636;
t638 = t466.*u1;
t639 = t352.*t433;
t640 = t385.*t432.*t434;
t782 = t363.*t432.*t435;
t641 = t639+t640-t782;
t642 = t343.*t433;
t643 = t365.*t432.*t435;
t644 = t189.*t203.*t432.*t434.*(1.0./2.0);
t645 = t642+t643+t644;
t646 = t346.*t433;
t647 = t373.*t432.*t435;
t785 = t381.*t432.*t434;
t648 = t646+t647-t785;
t649 = t648.*u2;
t779 = t630.*u8;
t781 = t637.*u7;
t783 = t641.*u4;
t784 = t645.*u5;
t650 = t627+t634+t638+t649-t779-t781-t783-t784;
t651 = t466.*t650;
t652 = t596+t624+t651;
t654 = t294.*t433;
t655 = t336.*t432.*t435;
t786 = t314.*t432.*t434;
t656 = t654+t655-t786;
t657 = t656.*u2;
t658 = t216.*t433;
t787 = t209.*t432.*t435;
t659 = t658-t787;
t660 = t659.*u6;
t661 = t209.*t432;
t662 = t222.*t433.*t434;
t788 = t216.*t433.*t435;
t663 = t661+t662-t788;
t664 = t332.*t432.*t435;
t665 = t310.*t432.*t434;
t790 = t301.*t433;
t666 = t664+t665-t790;
t667 = t666.*u3;
t668 = t216.*t432.*t434;
t669 = t222.*t432.*t435;
t670 = t668+t669;
t671 = t325.*t432.*t435;
t672 = t317.*t432.*t434;
t793 = t297.*t433;
t673 = t671+t672-t793;
t674 = t673.*u4;
t675 = t289.*t433;
t676 = t119.*t203.*t432.*t434.*(1.0./2.0);
t794 = t432.*t435.*t475;
t677 = t675+t676-t794;
t789 = t663.*u8;
t791 = t670.*u7;
t792 = t470.*u1;
t795 = t677.*u5;
t678 = t657+t660+t667+t674-t789-t791-t792-t795;
t679 = t473.*t678;
t680 = t314.*t448;
t681 = t336.*t446;
t796 = t294.*t432.*t436;
t682 = t680+t681-t796;
t683 = t682.*u2;
t684 = t332.*t446;
t685 = t301.*t432.*t436;
t798 = t310.*t448;
t686 = t684+t685-t798;
t687 = t686.*u3;
t688 = t216.*t448;
t800 = t222.*t446;
t689 = t688-t800;
t690 = t689.*u7;
t691 = t209.*t446;
t692 = t216.*t432.*t436;
t693 = t691+t692;
t694 = t325.*t446;
t695 = t297.*t432.*t436;
t802 = t317.*t448;
t696 = t694+t695-t802;
t697 = t696.*u4;
t698 = t446.*(t327-t401);
t699 = t119.*t203.*t448.*(1.0./2.0);
t700 = t289.*t432.*t436;
t701 = t698+t699+t700;
t702 = t701.*u5;
t703 = t209.*t433.*t436;
t704 = t216.*t432.*t435.*t436;
t803 = t222.*t432.*t434.*t436;
t705 = t703+t704-t803;
t797 = t460.*u9;
t799 = t463.*u1;
t801 = t693.*u6;
t804 = t705.*u8;
t706 = t683+t687+t690+t697+t702-t797-t799-t801-t804;
t707 = t478.*t706;
t708 = t314.*t440;
t709 = t336.*t442;
t710 = t294.*t432.*t437;
t711 = t708+t709+t710;
t712 = t711.*u2;
t713 = t310.*t440;
t714 = t301.*t432.*t437;
t805 = t332.*t442;
t715 = t713+t714-t805;
t716 = t456.*u9;
t717 = t216.*t440;
t808 = t222.*t442;
t718 = t717-t808;
t719 = t718.*u7;
t720 = t209.*t442;
t809 = t216.*t432.*t437;
t721 = t720-t809;
t722 = t317.*t440;
t723 = t297.*t432.*t437;
t811 = t325.*t442;
t724 = t722+t723-t811;
t725 = t209.*t433.*t437;
t726 = t216.*t432.*t435.*t437;
t813 = t222.*t432.*t434.*t437;
t727 = t725+t726-t813;
t728 = t727.*u8;
t729 = t442.*t475;
t730 = t119.*t203.*t440.*(1.0./2.0);
t814 = t289.*t432.*t437;
t731 = t729+t730-t814;
t732 = t731.*u5;
t806 = t715.*u3;
t807 = t453.*u1;
t810 = t721.*u6;
t812 = t724.*u4;
t733 = t483.*(t712+t716+t719+t728+t732-t806-t807-t810-t812);
t734 = t679+t707+t733;
t815 = cos(q10);
t816 = sin(q10);
t817 = t113.*t815.*(1.0./2.0);
t818 = t29.*t203.*t815.*(1.0./2.0);
t832 = t189.*t816;
t819 = t817+t818-t832;
t820 = t122.*t815.*(1.0./2.0);
t821 = t37.*t203.*t815.*(1.0./2.0);
t823 = t119.*t816;
t822 = t820+t821-t823;
t824 = t189.*t815;
t825 = t113.*t816.*(1.0./2.0);
t826 = t29.*t203.*t816.*(1.0./2.0);
t827 = t824+t825+t826;
t828 = t119.*t815;
t829 = t122.*t816.*(1.0./2.0);
t830 = t37.*t203.*t816.*(1.0./2.0);
t831 = t828+t829+t830;
t833 = -t220+t221+t286;
t834 = t116.*t815.*(1.0./2.0);
t835 = t19.*t203.*t815.*(1.0./2.0);
t837 = t127.*t816;
t836 = t834+t835-t837;
t838 = t228-t229+t230;
t839 = t127.*t815;
t840 = t116.*t816.*(1.0./2.0);
t841 = t19.*t203.*t816.*(1.0./2.0);
t842 = t839+t840+t841;
t843 = t19.*t109.*t815;
t844 = t19.*t110.*t816.*(1.0./2.0);
t997 = t12.*t203.*t816.*(1.0./2.0);
t845 = t843+t844-t997;
t846 = t845.*u4;
t847 = t161.*t816.*(1.0./2.0);
t848 = t6.*t7.*t9.*t203.*t816.*(1.0./2.0);
t998 = t151.*t815;
t849 = t847+t848-t998;
t850 = t849.*u3;
t851 = t155.*t815;
t852 = t164.*t816.*(1.0./2.0);
t853 = t17.*t203.*t816.*(1.0./2.0);
t854 = t851+t852+t853;
t855 = t836.*u10;
t856 = t116.*t815;
t1000 = t127.*t816.*(1.0./2.0);
t857 = t856-t1000;
t999 = t854.*u2;
t1001 = t857.*u5;
t858 = t846+t850+t855-t999-t1001;
t859 = t253-t254;
t860 = t859.*u3;
t861 = t257+t258+t424;
t862 = t861.*u4;
t863 = t260-t261+t262;
t864 = t863.*u2;
t865 = t264+t860+t862+t864;
t866 = t838.*t865;
t867 = t19.*t109.*t816;
t868 = t12.*t203.*t815.*(1.0./2.0);
t1003 = t19.*t110.*t815.*(1.0./2.0);
t869 = t867+t868-t1003;
t870 = t869.*u4;
t871 = t161.*t815.*(1.0./2.0);
t872 = t151.*t816;
t873 = t6.*t7.*t9.*t203.*t815.*(1.0./2.0);
t874 = t871+t872+t873;
t875 = t164.*t815.*(1.0./2.0);
t876 = t17.*t203.*t815.*(1.0./2.0);
t1005 = t155.*t816;
t877 = t875+t876-t1005;
t878 = t877.*u2;
t879 = t127.*t815.*(1.0./2.0);
t880 = t116.*t816;
t881 = t879+t880;
t882 = t842.*u10;
t1004 = t874.*u3;
t1006 = t881.*u5;
t883 = t870+t878+t882-t1004-t1006;
t884 = t819.*t883;
t1002 = t827.*t858;
t885 = t866+t884-t1002;
t886 = t119.*t815.*(1.0./2.0);
t887 = t122.*t816;
t888 = t886+t887;
t889 = t819.*u1;
t890 = t831.*u10;
t891 = t184.*t815.*(1.0./2.0);
t892 = t2.*t19.*t203.*t815.*(1.0./2.0);
t977 = t174.*t816;
t893 = t891+t892-t977;
t894 = t37.*t109.*t816;
t895 = t35.*t203.*t815.*(1.0./2.0);
t979 = t37.*t110.*t815.*(1.0./2.0);
t896 = t894+t895-t979;
t897 = t896.*u4;
t898 = t171.*t816;
t899 = t84.*t203.*t815.*(1.0./2.0);
t980 = t180.*t815.*(1.0./2.0);
t900 = t898+t899-t980;
t901 = t900.*u3;
t976 = t888.*u5;
t978 = t893.*u2;
t902 = t889+t890+t897+t901-t976-t978;
t903 = -t305+t306+t340;
t904 = t308+t309+t397;
t905 = t904.*u3;
t906 = t312-t313;
t907 = t906.*u2;
t908 = t315+t316+t399;
t909 = t908.*u4;
t910 = t838.*u1;
t911 = t400+t905+t907+t909+t910;
t912 = t903.*t911;
t913 = t37.*t109.*t815;
t914 = t37.*t110.*t816.*(1.0./2.0);
t982 = t35.*t203.*t816.*(1.0./2.0);
t915 = t913+t914-t982;
t916 = t915.*u4;
t917 = t119.*t816.*(1.0./2.0);
t983 = t122.*t815;
t1050 = t917-t983;
t918 = t1050.*u5;
t919 = t171.*t815;
t920 = t180.*t816.*(1.0./2.0);
t984 = t84.*t203.*t816.*(1.0./2.0);
t921 = t919+t920-t984;
t922 = t921.*u3;
t923 = t822.*u10;
t924 = t174.*t815;
t925 = t184.*t816.*(1.0./2.0);
t926 = t2.*t19.*t203.*t816.*(1.0./2.0);
t927 = t924+t925+t926;
t928 = t927.*u2;
t985 = t827.*u1;
t929 = t916+t918+t922+t923+t928-t985;
t930 = t842.*t929;
t981 = t836.*t902;
t931 = t912+t930-t981;
t932 = t189.*t815.*(1.0./2.0);
t933 = t113.*t816;
t934 = t932+t933;
t935 = t934.*u5;
t936 = t143.*t815.*(1.0./2.0);
t937 = t3.*t19.*t203.*t815.*(1.0./2.0);
t986 = t132.*t816;
t938 = t936+t937-t986;
t939 = t938.*u2;
t940 = t822.*u1;
t941 = t29.*t109.*t816;
t942 = t25.*t203.*t815.*(1.0./2.0);
t988 = t29.*t110.*t815.*(1.0./2.0);
t943 = t941+t942-t988;
t944 = t140.*t815.*(1.0./2.0);
t945 = t217.*t816;
t990 = t55.*t203.*t815.*(1.0./2.0);
t946 = t944+t945-t990;
t987 = t827.*u10;
t989 = t943.*u4;
t991 = t946.*u3;
t947 = t935+t939+t940-t987-t989-t991;
t948 = t822.*t947;
t949 = t29.*t109.*t815;
t950 = t29.*t110.*t816.*(1.0./2.0);
t992 = t25.*t203.*t816.*(1.0./2.0);
t951 = t949+t950-t992;
t952 = t951.*u4;
t953 = t113.*t815;
t993 = t189.*t816.*(1.0./2.0);
t954 = t953-t993;
t955 = t217.*t815;
t956 = t55.*t203.*t816.*(1.0./2.0);
t995 = t140.*t816.*(1.0./2.0);
t1031 = t955+t956-t995;
t957 = t1031.*u3;
t958 = t132.*t815;
t959 = t143.*t816.*(1.0./2.0);
t960 = t3.*t19.*t203.*t816.*(1.0./2.0);
t961 = t958+t959+t960;
t962 = t961.*u2;
t963 = t819.*u10;
t964 = t831.*u1;
t994 = t954.*u5;
t965 = t952+t957+t962+t963+t964-t994;
t966 = t831.*t965;
t967 = t376-t377+t417;
t968 = t379-t380;
t969 = t968.*u2;
t970 = t383-t384+t419;
t971 = t970.*u4;
t972 = t833.*u1;
t996 = t967.*u3;
t973 = t388+t969+t971+t972-t996;
t974 = t833.*t973;
t975 = t948+t966+t974;
t1007 = sin(q12);
t1008 = cos(q12);
t1009 = sin(q11);
t1010 = cos(q11);
t1011 = sin(q13);
t1012 = cos(q13);
t1013 = t1010.*t1011;
t1014 = t1008.*t1009.*t1012;
t1015 = t1013+t1014;
t1016 = t1009.*t1011;
t1018 = t1008.*t1010.*t1012;
t1017 = t1016-t1018;
t1019 = t1009.*t1012;
t1020 = t1008.*t1010.*t1011;
t1021 = t1019+t1020;
t1022 = t1010.*t1012;
t1024 = t1008.*t1009.*t1011;
t1023 = t1022-t1024;
t1025 = t827.*t1017;
t1026 = t838.*t1015;
t1058 = t819.*t1007.*t1012;
t1027 = t1025+t1026-t1058;
t1028 = t831.*t1021;
t1029 = t822.*t1007.*t1011;
t1059 = t833.*t1023;
t1030 = t1028+t1029-t1059;
t1032 = t833.*t1015;
t1033 = t822.*t1007.*t1012;
t1054 = t831.*t1017;
t1034 = t1032+t1033-t1054;
t1035 = t827.*t1021;
t1036 = t838.*t1023;
t1037 = t819.*t1007.*t1011;
t1038 = t1035+t1036+t1037;
t1039 = t833.*t1007.*t1009;
t1040 = t831.*t1007.*t1010;
t1042 = t822.*t1008;
t1041 = t1039+t1040-t1042;
t1043 = t819.*t1008;
t1044 = t838.*t1007.*t1009;
t1049 = t827.*t1007.*t1010;
t1045 = t1043+t1044-t1049;
t1046 = t903.*t1007.*t1009;
t1047 = t842.*t1007.*t1010;
t1232 = t836.*t1008;
t1048 = t1046+t1047-t1232;
t1051 = t842.*t1021;
t1052 = t836.*t1007.*t1011;
t1098 = t903.*t1023;
t1053 = t1051+t1052-t1098;
t1055 = t903.*t1015;
t1056 = t836.*t1007.*t1012;
t1135 = t842.*t1017;
t1057 = t1055+t1056-t1135;
t1060 = t842.*t1008;
t1061 = t836.*t1007.*t1010;
t1062 = t1060+t1061;
t1063 = t1062.*u10;
t1064 = t877.*t1008;
t1065 = t863.*t1007.*t1009;
t1320 = t854.*t1007.*t1010;
t1066 = t1064+t1065-t1320;
t1067 = t1066.*u2;
t1068 = t836.*t1007;
t1069 = t903.*t1008.*t1009;
t1070 = t842.*t1008.*t1010;
t1071 = t1068+t1069+t1070;
t1072 = t1071.*u12;
t1073 = t842.*t1007.*t1009;
t1321 = t903.*t1007.*t1010;
t1074 = t1073-t1321;
t1075 = t869.*t1008;
t1076 = t845.*t1007.*t1010;
t1077 = t861.*t1007.*t1009;
t1078 = t1075+t1076+t1077;
t1079 = t1078.*u4;
t1080 = t881.*t1008;
t1081 = t857.*t1007.*t1010;
t1323 = t127.*t203.*t1007.*t1009.*(1.0./2.0);
t1082 = t1080+t1081-t1323;
t1083 = t859.*t1007.*t1009;
t1084 = t849.*t1007.*t1010;
t1325 = t874.*t1008;
t1085 = t1083+t1084-t1325;
t1086 = t1085.*u3;
t1322 = t1074.*u11;
t1324 = t1082.*u5;
t1087 = t1063+t1067+t1072+t1079+t1086-t1322-t1324;
t1088 = t903.*t1007.*t1009.*t1012;
t1089 = t842.*t1007.*t1010.*t1012;
t1327 = t836.*t1008.*t1012;
t1090 = t1088+t1089-t1327;
t1091 = t1090.*u12;
t1092 = t859.*t1015;
t1093 = t874.*t1007.*t1012;
t1328 = t849.*t1017;
t1094 = t1092+t1093-t1328;
t1095 = t854.*t1017;
t1096 = t863.*t1015;
t1330 = t877.*t1007.*t1012;
t1097 = t1095+t1096-t1330;
t1099 = t1053.*u13;
t1100 = t842.*t1015;
t1101 = t903.*t1017;
t1102 = t1100+t1101;
t1103 = t1102.*u11;
t1104 = t836.*t1017;
t1105 = t842.*t1007.*t1012;
t1106 = t1104+t1105;
t1107 = t1106.*u10;
t1108 = t845.*t1017;
t1109 = t869.*t1007.*t1012;
t1332 = t861.*t1015;
t1110 = t1108+t1109-t1332;
t1111 = t1110.*u4;
t1112 = t857.*t1017;
t1113 = t127.*t203.*t1015.*(1.0./2.0);
t1114 = t881.*t1007.*t1012;
t1115 = t1112+t1113+t1114;
t1329 = t1094.*u3;
t1331 = t1097.*u2;
t1333 = t1115.*u5;
t1116 = t1091+t1099+t1103+t1107+t1111-t1329-t1331-t1333;
t1117 = t1027.*t1116;
t1118 = t861.*t1023;
t1119 = t869.*t1007.*t1011;
t1334 = t845.*t1021;
t1120 = t1118+t1119-t1334;
t1121 = t857.*t1021;
t1122 = t127.*t203.*t1023.*(1.0./2.0);
t1336 = t881.*t1007.*t1011;
t1123 = t1121+t1122-t1336;
t1124 = t849.*t1021;
t1125 = t874.*t1007.*t1011;
t1338 = t859.*t1023;
t1126 = t1124+t1125-t1338;
t1127 = t1126.*u3;
t1128 = t903.*t1007.*t1009.*t1011;
t1129 = t842.*t1007.*t1010.*t1011;
t1339 = t836.*t1008.*t1011;
t1130 = t1128+t1129-t1339;
t1131 = t854.*t1021;
t1132 = t863.*t1023;
t1133 = t877.*t1007.*t1011;
t1134 = t1131+t1132+t1133;
t1136 = t1057.*u13;
t1137 = t842.*t1023;
t1138 = t903.*t1021;
t1139 = t1137+t1138;
t1140 = t1139.*u11;
t1141 = t836.*t1021;
t1342 = t842.*t1007.*t1011;
t1142 = t1141-t1342;
t1143 = t1142.*u10;
t1335 = t1120.*u4;
t1337 = t1123.*u5;
t1340 = t1130.*u12;
t1341 = t1134.*u2;
t1144 = t1127+t1136+t1140+t1143-t1335-t1337-t1340-t1341;
t1145 = t1038.*t1144;
t1326 = t1045.*t1087;
t1146 = t1117+t1145-t1326;
t1147 = t961.*t1021;
t1148 = t938.*t1007.*t1011;
t1343 = t968.*t1023;
t1149 = t1147+t1148-t1343;
t1150 = t1149.*u2;
t1151 = t967.*t1023;
t1152 = t1021.*(t955+t956-t995);
t1345 = t946.*t1007.*t1011;
t1153 = t1151+t1152-t1345;
t1154 = t1153.*u3;
t1155 = t1030.*u1;
t1156 = t827.*t1023;
t1346 = t838.*t1021;
t1157 = t1156-t1346;
t1158 = t1157.*u11;
t1159 = t819.*t1021;
t1347 = t827.*t1007.*t1011;
t1160 = t1159-t1347;
t1161 = t1160.*u10;
t1162 = t970.*t1023;
t1163 = t943.*t1007.*t1011;
t1348 = t951.*t1021;
t1164 = t1162+t1163-t1348;
t1165 = t954.*t1021;
t1166 = t189.*t203.*t1023.*(1.0./2.0);
t1350 = t934.*t1007.*t1011;
t1167 = t1165+t1166-t1350;
t1168 = t819.*t1008.*t1011;
t1169 = t838.*t1007.*t1009.*t1011;
t1352 = t827.*t1007.*t1010.*t1011;
t1170 = t1168+t1169-t1352;
t1171 = t1170.*u12;
t1344 = t1027.*u13;
t1349 = t1164.*u4;
t1351 = t1167.*u5;
t1172 = t1150+t1154+t1155+t1158+t1161+t1171-t1344-t1349-t1351;
t1173 = t1030.*t1172;
t1174 = t968.*t1015;
t1175 = t938.*t1007.*t1012;
t1353 = t961.*t1017;
t1176 = t1174+t1175-t1353;
t1177 = t967.*t1015;
t1178 = t1017.*t1031;
t1179 = t946.*t1007.*t1012;
t1180 = t1177+t1178+t1179;
t1181 = t1180.*u3;
t1182 = t1038.*u13;
t1183 = t827.*t1015;
t1356 = t838.*t1017;
t1184 = t1183-t1356;
t1185 = t1184.*u11;
t1186 = t819.*t1017;
t1187 = t827.*t1007.*t1012;
t1188 = t1186+t1187;
t1189 = t1188.*u10;
t1190 = t951.*t1017;
t1191 = t943.*t1007.*t1012;
t1357 = t970.*t1015;
t1192 = t1190+t1191-t1357;
t1193 = t1192.*u4;
t1194 = t819.*t1008.*t1012;
t1195 = t838.*t1007.*t1009.*t1012;
t1358 = t827.*t1007.*t1010.*t1012;
t1196 = t1194+t1195-t1358;
t1197 = t954.*t1017;
t1198 = t934.*t1007.*t1012;
t1199 = t189.*t203.*t1015.*(1.0./2.0);
t1200 = t1197+t1198+t1199;
t1354 = t1176.*u2;
t1355 = t1034.*u1;
t1359 = t1196.*u12;
t1360 = t1200.*u5;
t1201 = t1181+t1182+t1185+t1189+t1193-t1354-t1355-t1359-t1360;
t1202 = t827.*t1008;
t1203 = t819.*t1007.*t1010;
t1204 = t1202+t1203;
t1205 = t1204.*u10;
t1206 = t819.*t1007;
t1207 = t827.*t1008.*t1010;
t1362 = t838.*t1008.*t1009;
t1208 = t1206+t1207-t1362;
t1209 = t1208.*u12;
t1210 = t946.*t1008;
t1211 = t1007.*t1010.*t1031;
t1363 = t967.*t1007.*t1009;
t1212 = u3.*(t1210+t1211-t1363);
t1213 = t827.*t1007.*t1009;
t1214 = t838.*t1007.*t1010;
t1215 = t1213+t1214;
t1216 = t1041.*u1;
t1217 = t943.*t1008;
t1218 = t951.*t1007.*t1010;
t1219 = t970.*t1007.*t1009;
t1220 = t1217+t1218+t1219;
t1221 = t1220.*u4;
t1222 = t934.*t1008;
t1223 = t954.*t1007.*t1010;
t1365 = t189.*t203.*t1007.*t1009.*(1.0./2.0);
t1224 = t1222+t1223-t1365;
t1225 = t961.*t1007.*t1010;
t1226 = t968.*t1007.*t1009;
t1367 = t938.*t1008;
t1227 = t1225+t1226-t1367;
t1228 = t1227.*u2;
t1364 = t1215.*u11;
t1366 = t1224.*u5;
t1229 = t1205+t1209+t1212+t1216+t1221+t1228-t1364-t1366;
t1230 = t1041.*t1229;
t1361 = t1034.*t1201;
t1231 = t1173+t1230-t1361;
t1233 = t927.*t1007.*t1010;
t1234 = t906.*t1007.*t1009;
t1368 = t893.*t1008;
t1235 = t1233+t1234-t1368;
t1236 = t1235.*u2;
t1237 = t831.*t1008;
t1238 = t822.*t1007.*t1010;
t1239 = t1237+t1238;
t1240 = t1239.*u10;
t1241 = t822.*t1007;
t1242 = t833.*t1008.*t1009;
t1243 = t831.*t1008.*t1010;
t1244 = t1241+t1242+t1243;
t1245 = t1244.*u12;
t1246 = t900.*t1008;
t1247 = t921.*t1007.*t1010;
t1248 = t904.*t1007.*t1009;
t1249 = t1246+t1247+t1248;
t1250 = t1249.*u3;
t1251 = t831.*t1007.*t1009;
t1369 = t833.*t1007.*t1010;
t1252 = t1251-t1369;
t1253 = t1045.*u1;
t1254 = t896.*t1008;
t1255 = t915.*t1007.*t1010;
t1256 = t908.*t1007.*t1009;
t1257 = t1254+t1255+t1256;
t1258 = t1257.*u4;
t1259 = t1007.*t1010.*t1050;
t1260 = t119.*t203.*t1007.*t1009.*(1.0./2.0);
t1371 = t888.*t1008;
t1261 = u5.*(t1259+t1260-t1371);
t1370 = t1252.*u11;
t1262 = t1236+t1240+t1245+t1250+t1253+t1258+t1261-t1370;
t1263 = t1048.*t1262;
t1264 = t927.*t1021;
t1265 = t893.*t1007.*t1011;
t1372 = t906.*t1023;
t1266 = t1264+t1265-t1372;
t1267 = t1266.*u2;
t1268 = t1034.*u13;
t1269 = t904.*t1023;
t1270 = t900.*t1007.*t1011;
t1373 = t921.*t1021;
t1271 = t1269+t1270-t1373;
t1272 = t831.*t1023;
t1273 = t833.*t1021;
t1274 = t1272+t1273;
t1275 = t1274.*u11;
t1276 = t822.*t1021;
t1376 = t831.*t1007.*t1011;
t1277 = t1276-t1376;
t1278 = t1277.*u10;
t1279 = t908.*t1023;
t1280 = t896.*t1007.*t1011;
t1377 = t915.*t1021;
t1281 = t1279+t1280-t1377;
t1282 = t1021.*t1050;
t1283 = t888.*t1007.*t1011;
t1379 = t119.*t203.*t1023.*(1.0./2.0);
t1284 = t1282+t1283-t1379;
t1285 = t1284.*u5;
t1286 = t831.*t1007.*t1010.*t1011;
t1287 = t833.*t1007.*t1009.*t1011;
t1380 = t822.*t1008.*t1011;
t1288 = t1286+t1287-t1380;
t1374 = t1271.*u3;
t1375 = t1038.*u1;
t1378 = t1281.*u4;
t1381 = t1288.*u12;
t1289 = t1267+t1268+t1275+t1278+t1285-t1374-t1375-t1378-t1381;
t1290 = t1053.*t1289;
t1291 = t906.*t1015;
t1292 = t893.*t1007.*t1012;
t1382 = t927.*t1017;
t1293 = t1291+t1292-t1382;
t1294 = t921.*t1017;
t1295 = t900.*t1007.*t1012;
t1384 = t904.*t1015;
t1296 = t1294+t1295-t1384;
t1297 = t1296.*u3;
t1298 = t1030.*u13;
t1299 = t831.*t1015;
t1300 = t833.*t1017;
t1301 = t1299+t1300;
t1302 = t1301.*u11;
t1303 = t822.*t1017;
t1304 = t831.*t1007.*t1012;
t1305 = t1303+t1304;
t1306 = t1305.*u10;
t1307 = t915.*t1017;
t1308 = t896.*t1007.*t1012;
t1386 = t908.*t1015;
t1309 = t1307+t1308-t1386;
t1310 = t1309.*u4;
t1311 = t831.*t1007.*t1010.*t1012;
t1312 = t833.*t1007.*t1009.*t1012;
t1387 = t822.*t1008.*t1012;
t1313 = t1311+t1312-t1387;
t1314 = t1313.*u12;
t1315 = t888.*t1007.*t1012;
t1316 = t119.*t203.*t1015.*(1.0./2.0);
t1388 = t1017.*t1050;
t1317 = t1315+t1316-t1388;
t1383 = t1293.*u2;
t1385 = t1027.*u1;
t1389 = t1317.*u5;
t1318 = t1297+t1298+t1302+t1306+t1310+t1314-t1383-t1385-t1389;
t1390 = t1057.*t1318;
t1319 = t1263+t1290-t1390;
wc = reshape([0.0,t32.*t73+t21.*t95+t6.*t9.*t61,t113.*(-t106+t188+t127.*(t172+t175+t176+t124.*u1+u5.*(t120-t35.*t110)))+t116.*t148+t122.*(-t72+t159+t169),-t209.*t285-t206.*t339+t225.*t391,-t466.*t567+t473.*t652+t470.*t734,t822.*t885-t819.*t931-t836.*t975,t1041.*t1146+t1048.*t1231-t1045.*t1319,t2.^2.*u1+t3.^2.*u1,t19.*t61-t37.*t73+t29.*t95,t19.*t148+t29.*(-t106+t188+t127.*(t172+t175+t176-t122.*u5+u1.*(t123-t125)))+t37.*(-t72+t159+t169),t222.*t285+t231.*t339-t307.*t391,t456.*t567-t478.*t652-t463.*t734,-t833.*t885-t838.*t931+t903.*t975,-t1030.*t1146-t1053.*t1231-t1038.*t1319,0.0,-t12.*t61+t35.*t73-t25.*t95,-t189.*(-t106+t188+t127.*(t172+t175+t176-t122.*u5-t189.*u1))+t119.*(t72-t169+t158.*t189)-t127.*(t60+t136+t122.*(t137+t144+t145-t196+t189.*u5)),t216.*t285+t213.*t339-t234.*t391,-t460.*t567+t483.*t652+t453.*t734,t831.*t885-t827.*t931-t842.*t975,-t1034.*t1146-t1057.*t1231+t1027.*(t1263+t1290-t1390)],[7, 3]);