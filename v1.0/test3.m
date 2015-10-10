% EOM Simulation:

%% properties:
lc1 = 1; lc2 = 1; lc3 = 1; lc4 = 1; % link COM pos.
lcn1 = 2; lcn2 = 2;
m1 = 1; m2 = 1; m3 = 1; m4 = 1; % mass
I1 = 1e-3; I2 = 1e-3; I3 = 1e-2; I4 = 1e-3; % rotational inertia
l1 = 1; l2 = 1; l3 = 1;
rb1 = [ 1 0 0 ]; rb2 = [ cos(pi/3) sin(pi/3) 0 ]; rb3 = [ cos(pi/3) -sin(pi/3) 0 ]; % joint pos. vec. in base from base centre
rp2 = [ -sin(pi/3) 1+cos(pi/3) 0 ]; rp3 = [ sin(pi/3) 1+cos(pi/3) 0 ]; % joint pos. vec. in platform from joint II
g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ 0 0 lc1 0 0 ; -lc2 0 0 0 0 ; 0 0 -lc3 0 0 ; 0 0 -lc4 0 2 ; 0 0 -lcn1 2 3 ; 0 0 -lcn2 2 4 ];
m = [ m1 , m2 , m3 , m4 ];
I = sym ( zeros ( 3 , 3 , 4 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 ); I(:,:,4) = I4 * eye ( 3 );

% j = sym ( zeros ( 4 , 5 , 6 ) );
% j(1,:,1) = [ 1 inf rb1 ]; j(2,:,1) = [ 3 inf 0 0 0 ]; j(3,:,1) = [ 2 inf 0 0 0 ]; % x-z-y
% j(1,:,2) = [ 2 inf 0 0 l1 ]; j(2,:,2) = [ 3 inf 0 0 0 ]; j(3,:,2) = [ 1 inf 0 0 0 ]; % y-z-x
% j(1,:,3) = [ 3 inf rp2 ]; j(2,:,3) = [ 1 inf 0 0 0 ]; j(3,:,3) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,4) = [ 3 inf rp3 ]; j(2,:,4) = [ 1 inf 0 0 0 ]; j(3,:,4) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,5) = [ 2 inf 0 0 -l2 ]; j(2,:,5) = [ 3 inf 0 0 0 ]; j(3,:,5) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,5) = [ 0 0 rb2 ];
% j(1,:,6) = [ 2 inf 0 0 -l3 ]; j(2,:,6) = [ 3 inf 0 0 0 ]; j(3,:,6) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,6) = [ 0 0 rb3 ];

j = sym ( zeros ( 3 , 5 , 6 ) );
j(1,:,1) = [ 1 inf rb1 ]; j(2,:,1) = [ 3 inf 0 0 0 ]; j(3,:,1) = [ 2 inf 0 0 0 ]; % x-z-y
j(1,:,2) = [ 2 inf 0 0 l1 ]; j(2,:,2) = [ 3 inf 0 0 0 ]; j(3,:,2) = [ 1 inf 0 0 0 ]; % y-z-x
j(1,:,3) = [ 3 inf rp2 ]; j(2,:,3) = [ 1 inf 0 0 0 ]; j(3,:,3) = [ 2 inf 0 0 0 ]; % z-x-y
j(1,:,4) = [ 3 inf rp3 ]; j(2,:,4) = [ 1 inf 0 0 0 ]; j(3,:,4) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,5) fixed end in 2nd and 3rd links frames
% j(1,:,6) 

jkd = sym (zeros ( 3 , 2 , 18 ) );

%% Derivation:

% matlabpool open local 4

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM3 ( M , T , Dd , fg , fj , Tef , Tcn , Dcn , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );

matlabpool close
