  rj[1][0] = cos(q1)*2.0;
  rj[1][2] = sin(q1)*-2.0;
  rj[2][0] = cos(q1+q2)*2.0+cos(q1)*2.0;
  rj[2][2] = sin(q1+q2)*-2.0-sin(q1)*2.0;
  rj[3][0] = cos(q1+q2+q3)+cos(q1+q2)*2.0+cos(q1)*2.0;
  rj[3][2] = -sin(q1+q2+q3)-sin(q1+q2)*2.0-sin(q1)*2.0;