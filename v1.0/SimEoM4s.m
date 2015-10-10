%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM4s ( Tcnt , q , u , dt )

format longE

global qt;
global ut;

qt = q; ut = u;

[ tmp nq ] = size ( q );
[ tmp , tmp , ncnv ] = size ( Tcnt );
m = nq + ncnv * 3; % # of states and constraint vectors (Lagrangian multiplyers)
par = [ nq ncnv*3 ];
z0 = 1e-3 * ones ( 1 , 2 * m );
% z0( 1, 1 : nq ) = pi/3 * ones ( 1 , nq );
% z0( 1, 1 : nq ) = 0.26*pi * [-1, 1e-3, 2, 1e-3, -1, 1, 1e-3, -2, 1e-3, 1, 1e-3, -2, 1e-3, 1, 1];
z0( 1, 1 : nq ) = 1/6*pi * [-1, 1, 2, -1, -1, 1, 1, -2, -1, 1, 1, -2, -1, 1, 1];
t0 = 0;

% Standard ODE solver:
options = odeset ();%'abstol',1*1e-6,'reltol',1*1e-6);
tspan = linspace ( t0 , t0 + dt , 500);
[ t , z , tfinal ] = ode15s ( @EOM , tspan , z0 , options , par );

% Simple Runge-Kutta
% st = 1e-3; z = z0; t = t0;
% for tc = t0 : st : dt
%     t = [ t , tc ];
%     z = [ z ; ( EOM ( tc , z(end,:)' , par ) )' * st ];
% end


function dz = EOM ( t , z , par )
t

global qt;
global ut;

nq = par(1); ncn = par (2);
ns = nq + ncn;
qu = [ qt , ut ];
zq = [ z( 1 : nq , 1 ) ; z( ns + 1 : ns + nq , 1 ) ]'; % only states
u = z( ns + 1 : end );

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); q7 = z(7);
q8 = z(8); q9 = z(9); q10 = z(10); q11 = z(11); q12 = z(12); q13 = z(13);
q14 = z(14); q15 = z(15);
u1 = z( ns + 1 ); u2 = z( ns + 2 ); u3 = z( ns + 3 ); u4 = z( ns + 4 );
u5 = z( ns + 5 ); u6 = z( ns + 6 ); u7 = z( ns + 7 ); u8 = z( ns + 8 );
u9 = z( ns + 9 ); u10 = z( ns + 10 ); u11 = z( ns + 11 ); u12 = z( ns + 12 );
u13 = z( ns + 13 ); u14 = z( ns + 14 ); u15 = z( ns + 15 );

% AnimEOM4 ( t , z' , rjf , qt , ut , nq , ncn );
% pause(1e-2);

Mq = M;
Tq = T(q1,q2,q3,q4,q5);
Dq = Dd(q1,q2,q3,q4,q5,u1,u2,u3,u4,u5);
fgq = fg;
fjq = fj;
Tefq = Tef;
Tcnq = Tcn(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15);
Dcnq = Dcn(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15);

A = zeros ( ns , ns );
A ( 1 : nq , : ) = [ Tq.'*Mq*Tq -Tcnq(:,:,1).' -Tcnq(:,:,2).' ];
A ( nq + 1 : nq + 3, 1 : nq ) = Tcnq(:,:,1);
A ( nq + 4 : nq + 6 , 1 : nq ) = Tcnq(:,:,2);

B = zeros ( ns , 1 );
B ( 1 : nq , 1 ) = Tq.' * ( fgq - Mq * Dq * u ( 1 : nq ) ) + fjq;
B ( nq + 1 : nq + 3 , 1 ) = - Dcnq(:,:,1) * u ( 1 : nq );
B ( nq + 4 : nq + 6 , 1 ) = - Dcnq(:,:,2) * u ( 1 : nq );

% % matrix inverse
% dzt = A \ B;
% % eig(A)

% SVD decomposition:
[ U , S , V ] = svd ( A ); % A=U*S*V' & A^-1=V*S^-1*U'
Si = inv ( S );
for i = 1 : ns
    if Si(i,i) > 1e2; Si(i,i) = 0; end
end
dzt = V * Si * U' * B;

% variation vector
dz = [ u ; dzt ];

