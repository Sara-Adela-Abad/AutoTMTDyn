  rj[0][0] = 1.0;
  rj[1][0] = cos(q2)*sin(q3)+1.0;
  rj[1][1] = -cos(q3)*sin(q1)+cos(q1)*sin(q2)*sin(q3);
  rj[1][2] = cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3);
  rj[2][0] = cos(q2)*sin(q3)+sqrt(3.0)*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))*(1.0/2.0)+cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-cos(q5)*cos(q6)*sin(q2)*(3.0/2.0)+1.0;
  rj[2][1] = sqrt(3.0)*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))*(-1.0/2.0)-cos(q3)*sin(q1)+(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)+cos(q1)*sin(q2)*sin(q3)+cos(q1)*cos(q2)*cos(q5)*cos(q6)*(3.0/2.0);
  rj[2][2] = cos(q1)*cos(q3)+sqrt(3.0)*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))*(1.0/2.0)+(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+sin(q1)*sin(q2)*sin(q3)+cos(q2)*cos(q5)*cos(q6)*sin(q1)*(3.0/2.0);
  rj[3][0] = cos(q2)*sin(q3)-sqrt(3.0)*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))*(1.0/2.0)+cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-cos(q5)*cos(q6)*sin(q2)*(3.0/2.0)+1.0;
  rj[3][1] = sqrt(3.0)*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))*(1.0/2.0)-cos(q3)*sin(q1)+(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)+cos(q1)*sin(q2)*sin(q3)+cos(q1)*cos(q2)*cos(q5)*cos(q6)*(3.0/2.0);
  rj[3][2] = cos(q1)*cos(q3)-sqrt(3.0)*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))*(1.0/2.0)+(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+sin(q1)*sin(q2)*sin(q3)+cos(q2)*cos(q5)*cos(q6)*sin(q1)*(3.0/2.0);
  rj[3][3] = 2.0;
  rj[4][0] = cos(q2)*sin(q3)-(sin(q7)*sin(q9)-cos(q7)*cos(q9)*sin(q8))*(cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))-cos(q5)*cos(q6)*sin(q2))+sqrt(3.0)*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))*(1.0/2.0)+(cos(q7)*sin(q9)+cos(q9)*sin(q7)*sin(q8))*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))-cos(q8)*cos(q9)*(cos(q2)*cos(q3)*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+cos(q2)*sin(q3)*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+cos(q5)*sin(q2)*sin(q6))+cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-cos(q5)*cos(q6)*sin(q2)*(3.0/2.0)+1.0;
  rj[4][1] = sqrt(3.0)*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))*(-1.0/2.0)-cos(q3)*sin(q1)-(cos(q7)*sin(q9)+cos(q9)*sin(q7)*sin(q8))*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))+(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(sin(q7)*sin(q9)-cos(q7)*cos(q9)*sin(q8))*((sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))+cos(q1)*cos(q2)*cos(q5)*cos(q6))+cos(q8)*cos(q9)*(-(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+cos(q1)*cos(q2)*cos(q5)*sin(q6))+cos(q1)*sin(q2)*sin(q3)+cos(q1)*cos(q2)*cos(q5)*cos(q6)*(3.0/2.0);
  rj[4][2] = cos(q1)*cos(q3)+sqrt(3.0)*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))*(1.0/2.0)+(cos(q7)*sin(q9)+cos(q9)*sin(q7)*sin(q8))*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))+(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(sin(q7)*sin(q9)-cos(q7)*cos(q9)*sin(q8))*((cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))+cos(q2)*cos(q5)*cos(q6)*sin(q1))+sin(q1)*sin(q2)*sin(q3)+cos(q8)*cos(q9)*(-(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+cos(q2)*cos(q5)*sin(q1)*sin(q6))+cos(q2)*cos(q5)*cos(q6)*sin(q1)*(3.0/2.0);
  rj[4][3] = 3.0;
  rj[5][0] = cos(q2)*sin(q3)-(sin(q10)*sin(q12)-cos(q10)*cos(q12)*sin(q11))*(cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))-cos(q5)*cos(q6)*sin(q2))-sqrt(3.0)*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))*(1.0/2.0)+(cos(q10)*sin(q12)+cos(q12)*sin(q10)*sin(q11))*(sin(q2)*sin(q5)-cos(q2)*cos(q3)*cos(q4)*cos(q5)+cos(q2)*cos(q5)*sin(q3)*sin(q4))-cos(q11)*cos(q12)*(cos(q2)*cos(q3)*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+cos(q2)*sin(q3)*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+cos(q5)*sin(q2)*sin(q6))+cos(q2)*cos(q3)*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)+cos(q2)*sin(q3)*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-cos(q5)*cos(q6)*sin(q2)*(3.0/2.0)+1.0;
  rj[5][1] = sqrt(3.0)*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))*(1.0/2.0)-cos(q3)*sin(q1)-(cos(q10)*sin(q12)+cos(q12)*sin(q10)*sin(q11))*(cos(q4)*cos(q5)*(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))+cos(q5)*sin(q4)*(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))+cos(q1)*cos(q2)*sin(q5))+(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(sin(q10)*sin(q12)-cos(q10)*cos(q12)*sin(q11))*((sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))-(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))+cos(q1)*cos(q2)*cos(q5)*cos(q6))+cos(q11)*cos(q12)*(-(sin(q1)*sin(q3)+cos(q1)*cos(q3)*sin(q2))*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+(cos(q3)*sin(q1)-cos(q1)*sin(q2)*sin(q3))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+cos(q1)*cos(q2)*cos(q5)*sin(q6))+cos(q1)*sin(q2)*sin(q3)+cos(q1)*cos(q2)*cos(q5)*cos(q6)*(3.0/2.0);
  rj[5][2] = cos(q1)*cos(q3)-sqrt(3.0)*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))*(1.0/2.0)+(cos(q10)*sin(q12)+cos(q12)*sin(q10)*sin(q11))*(cos(q4)*cos(q5)*(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))+cos(q5)*sin(q4)*(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))-cos(q2)*sin(q1)*sin(q5))+(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))*(3.0/2.0)-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))*(3.0/2.0)-(sin(q10)*sin(q12)-cos(q10)*cos(q12)*sin(q11))*((cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*sin(q6)+cos(q6)*sin(q4)*sin(q5))-(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(sin(q4)*sin(q6)-cos(q4)*cos(q6)*sin(q5))+cos(q2)*cos(q5)*cos(q6)*sin(q1))+sin(q1)*sin(q2)*sin(q3)+cos(q11)*cos(q12)*(-(cos(q1)*cos(q3)+sin(q1)*sin(q2)*sin(q3))*(cos(q4)*cos(q6)-sin(q4)*sin(q5)*sin(q6))+(cos(q1)*sin(q3)-cos(q3)*sin(q1)*sin(q2))*(cos(q6)*sin(q4)+cos(q4)*sin(q5)*sin(q6))+cos(q2)*cos(q5)*sin(q1)*sin(q6))+cos(q2)*cos(q5)*cos(q6)*sin(q1)*(3.0/2.0);
  rj[5][3] = 4.0;