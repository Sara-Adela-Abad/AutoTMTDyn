% EOM Simulation:

% properties:
lc1 = 1; lc2 = 1; lc3 = 1; % link COM pos.
m1 = 1; m2 = 1; m3 = 1; % mass
I1 = 1e-3; I2 = 1e-3; I3 = 1e-3; % rotational inertia
l1 = 2; l2 = 2; % link length
k1 = 1e-1; k2 = 1e-1; k3 = 1e-1; % joint spring coeff.
g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ lc1 0 0 ; lc2 0 0; lc3 0 0 ];
m = [ m1 , m2 , m3 ];
I = zeros ( 3 , 3 , 3 );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 );
j = zeros ( 1 , 5 , 3 );
j(:,:,1) = [ 2 inf 0 0 0 ]; j(:,:,2) = [ 2 inf l1 0 0 ]; j(:,:,3) = [ 2 inf l2 0 0 ];
jkd = zeros ( 3 , 2 , 3 );
jkd(1,:,1) = [ k1 0 ]; jkd(1,:,2) = [ k2 0 ]; jkd(1,:,3) = [ k3 0 ];

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 7 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );

