  Dd[0][0] = -u1*cos(q1);
  Dd[2][0] = u1*sin(q1);
  Dd[3][0] = -u1*(cos(q1+q2)+cos(q1)*2.0)-u2*cos(q1+q2);
  Dd[3][1] = -cos(q1+q2)*(u1+u2);
  Dd[5][0] = u1*(sin(q1+q2)+sin(q1)*2.0)+u2*sin(q1+q2);
  Dd[5][1] = sin(q1+q2)*(u1+u2);
  Dd[6][0] = -u2*(cos(q1+q2+q3)+cos(q1+q2)*2.0)-u1*(cos(q1+q2+q3)+cos(q1+q2)*2.0+cos(q1)*2.0)-u3*cos(q1+q2+q3);
  Dd[6][1] = -u1*(cos(q1+q2+q3)+cos(q1+q2)*2.0)-u2*(cos(q1+q2+q3)+cos(q1+q2)*2.0)-u3*cos(q1+q2+q3);
  Dd[6][2] = -cos(q1+q2+q3)*(u1+u2+u3);
  Dd[8][0] = u2*(sin(q1+q2+q3)+sin(q1+q2)*2.0)+u1*(sin(q1+q2+q3)+sin(q1+q2)*2.0+sin(q1)*2.0)+u3*sin(q1+q2+q3);
  Dd[8][1] = u1*(sin(q1+q2+q3)+sin(q1+q2)*2.0)+u2*(sin(q1+q2+q3)+sin(q1+q2)*2.0)+u3*sin(q1+q2+q3);
  Dd[8][2] = sin(q1+q2+q3)*(u1+u2+u3);