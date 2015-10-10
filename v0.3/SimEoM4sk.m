%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM4sk ( Tcnt , q , u , dt )

format longE

global qt;
global ut;

qt = q; ut = u;

[ tmp nq ] = size ( q );
[ tmp , tmp , ncnv ] = size ( Tcnt );
par = [ nq ncnv*3 ];
z0 = 0 * ones ( 1 , nq );
z0( 1, 1 : nq ) = pi/6 * [-1, 1e-3, 2, -1e-3, -1, 1, -1e-3, -2, 1e-3, 1, 1 , -1e-3, -2, 1e-3, 1];
t0 = 0;

% Standard ODE solver:
options = odeset ();%'abstol',1*1e-6,'reltol',1*1e-6);
tspan = linspace ( t0 , t0 + dt , 500);
[ t , z , tfinal ] = ode113 ( @EOM , tspan , z0 , options , par );

% Simple Runge-Kutta
% st = 1e-3; z = z0; t = t0;
% for tc = t0 : st : dt
%     t = [ t , tc ];
%     z = [ z ; ( EOM ( tc , z(end,:)' , par ) )' * st ];
% end


function dz = EOM ( t , z , par )
% t

nq = par(1); ncn = par (2);

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); q7 = z(7);
q8 = z(8); q9 = z(9); q10 = z(10); q11 = z(11); q12 = z(12); q13 = z(13);
q14 = z(14); q15 = z(15);

% AnimEOM4 ( t , z' , rjf , qt , ut , nq , ncn );
% pause(1e-2);

Tq = T(q1,q2,q3,q4,q5);
Tcnq = Tcn(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15);

% [-1, 1e-3, 2, -1e-3, -1, 1, -1e-3, -2, 1e-3, 1, 1 , -1e-3, -2, 1e-3, 1]
A = zeros ( nq , nq );
A ( 1 : 3 , : ) = Tcnq(:,:,1);
A ( 4 : 6 , : ) = Tcnq(:,:,2);
A(7,1) = 1; A(7,5) = 1;
A(8,2) = 1; A(8,4) = -1;
A(9,6) = 1; A(9,10) = 1;
A(10,7) = 1; A(10,9) = -1;
A(11,11) = 1; A(11,15) = 1;
A(12,12) = 1; A(12,14) = -1;
A ( 13 : 15 , : ) = Tq ( 1 : 3 , : );

w = 15;
B = zeros ( nq , 1 );
B ( 13 : 15 , 1 ) = 0.07 * w * [ cos(w*t) ; -sin(w*t) ; 2*cos(2*w*t) ];

% matrix inverse
dz = A \ B;
% eig(A)

% % SVD decomposition:
% [ U , S , V ] = svd ( A ); % A=U*S*V' & A^-1=V*S^-1*U'
% Si = inv ( S );
% for i = 1 : ns
%     if Si(i,i) > 1e2; Si(i,i) = 0; end
% end
% dz = V * Si * U' * B;

