  T[0][0] = -sin(q1);
  T[2][0] = -cos(q1);
  T[3][0] = -sin(q1+q2)-sin(q1)*2.0;
  T[3][1] = -sin(q1+q2);
  T[5][0] = -cos(q1+q2)-cos(q1)*2.0;
  T[5][1] = -cos(q1+q2);
  T[6][0] = -sin(q1+q2+q3)-sin(q1+q2)*2.0-sin(q1)*2.0;
  T[6][1] = -sin(q1+q2+q3)-sin(q1+q2)*2.0;
  T[6][2] = -sin(q1+q2+q3);
  T[8][0] = -cos(q1+q2+q3)-cos(q1+q2)*2.0-cos(q1)*2.0;
  T[8][1] = -cos(q1+q2+q3)-cos(q1+q2)*2.0;
  T[8][2] = -cos(q1+q2+q3);
  T[10][0] = 1.0;
  T[13][0] = 1.0;
  T[13][1] = 1.0;
  T[16][0] = 1.0;
  T[16][1] = 1.0;
  T[16][2] = 1.0;