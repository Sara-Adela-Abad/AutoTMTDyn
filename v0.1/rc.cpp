  rc[0][0] = cos(q1);
  rc[0][2] = -sin(q1);
  rc[1][0] = cos(q1+q2)+cos(q1)*2.0;
  rc[1][2] = -sin(q1+q2)-sin(q1)*2.0;
  rc[2][0] = cos(q1+q2+q3)+cos(q1+q2)*2.0+cos(q1)*2.0;
  rc[2][2] = -sin(q1+q2+q3)-sin(q1+q2)*2.0-sin(q1)*2.0;