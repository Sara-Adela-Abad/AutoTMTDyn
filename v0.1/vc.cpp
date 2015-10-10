  vc[0][0] = -u1*sin(q1);
  vc[0][2] = -u1*cos(q1);
  vc[1][0] = -u1*(sin(q1+q2)+sin(q1)*2.0)-u2*sin(q1+q2);
  vc[1][2] = -u1*(cos(q1+q2)+cos(q1)*2.0)-u2*cos(q1+q2);
  vc[2][0] = -u2*(sin(q1+q2+q3)+sin(q1+q2)*2.0)-u1*(sin(q1+q2+q3)+sin(q1+q2)*2.0+sin(q1)*2.0)-u3*sin(q1+q2+q3);
  vc[2][2] = -u2*(cos(q1+q2+q3)+cos(q1+q2)*2.0)-u1*(cos(q1+q2+q3)+cos(q1+q2)*2.0+cos(q1)*2.0)-u3*cos(q1+q2+q3);