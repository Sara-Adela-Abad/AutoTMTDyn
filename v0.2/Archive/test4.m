% EOM Simulation:

%% properties:
lc = 0.5; % link COM pos.
l = 1; % Link lenght
ml = 1; % mass
Il = 1e-3; % rotational inertia
rb01 = [ 1 0 0 ]; rb02 = [ -sin(pi/6) cos(pi/6) 0 ]; rb03 = [ -sin(pi/6) -cos(pi/6) 0 ]; % joint pos. vec. in base from base centre
lcp = 1; rp12 = [ -1-sin(pi/6) cos(pi/6) 0 ]; rp13 = [ -1-sin(pi/6) -cos(pi/6) 0 ]; % joint pos. vec. in platform from joint II
g = [ 0 , 0 , -1 ]; % gravity vector
% pi = 3.1415;

%inputs:
lc = [ lc 0 0 0 0 ; -lc 0 0 0 0 ; -lcp 0 0 0 0 ; lc 0 0 0 0 ; -lc 0 0 0 0 ; lc 0 0 0 3 ; -lc 0 0 0 0 ; 0 0 0 2 0 ; 0 0 0 2 5 ];
m = ml  * ones(1,7);
I = sym ( zeros ( 3 , 3 , 7 ) );
I(:,:,1) = Il * eye ( 3 ); I(:,:,2) = Il * eye ( 3 ); I(:,:,3) = 2 * Il * eye ( 3 ); I(:,:,4) = Il * eye ( 3 ); I(:,:,5) = Il * eye ( 3 ); I(:,:,6) = Il * eye ( 3 ); I(:,:,7) = Il * eye ( 3 );

% j = sym ( zeros ( 4 , 5 , 6 ) );
% j(1,:,1) = [ 1 inf rb1 ]; j(2,:,1) = [ 3 inf 0 0 0 ]; j(3,:,1) = [ 2 inf 0 0 0 ]; % x-z-y
% j(1,:,2) = [ 2 inf 0 0 l1 ]; j(2,:,2) = [ 3 inf 0 0 0 ]; j(3,:,2) = [ 1 inf 0 0 0 ]; % y-z-x
% j(1,:,3) = [ 3 inf rp2 ]; j(2,:,3) = [ 1 inf 0 0 0 ]; j(3,:,3) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,4) = [ 3 inf rp3 ]; j(2,:,4) = [ 1 inf 0 0 0 ]; j(3,:,4) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,5) = [ 2 inf 0 0 -l2 ]; j(2,:,5) = [ 3 inf 0 0 0 ]; j(3,:,5) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,5) = [ 0 0 rb2 ];
% j(1,:,6) = [ 2 inf 0 0 -l3 ]; j(2,:,6) = [ 3 inf 0 0 0 ]; j(3,:,6) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,6) = [ 0 0 rb3 ];

j = sym ( zeros ( 3 , 5 , 9 ) );
j(1,:,1) = [ 2 inf rb01 ]; % y
j(1,:,2) = [ 1 inf l 0 0 ]; j(2,:,2) = [ 2 inf 0 0 0 ]; j(3,:,2) = [ 1 inf 0 0 0 ]; % x-y-x
j(1,:,3) = [ 2 inf -l 0 0 ]; % y
j(1,:,4) = [ 3 2*pi/3 rp12 ]; j(2,:,4) = [ 2 inf 0 0 0 ]; % z-y
j(1,:,5) = [ 1 inf l 0 0 ]; j(2,:,5) = [ 2 inf 0 0 0 ]; j(3,:,5) = [ 1 inf 0 0 0 ]; % x-y-x

j(1,:,6) = [ 3 -2*pi/3 rp13 ]; j(2,:,6) = [ 2 inf 0 0 0 ]; % z-y
j(1,:,7) = [ 1 inf l 0 0 ]; j(2,:,7) = [ 2 inf 0 0 0 ]; j(3,:,7) = [ 1 inf 0 0 0 ]; % x-y-x

j(1,:,8) = [ 2 inf -l 0 0 ]; % y
j(1,:,9) = [ 2 inf -l 0 0 ]; % y


jkd = sym (zeros ( 3 , 2 , 15 ) );

%% Derivation:

% matlabpool open local 4

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );
 
% numerical simulation
[ t , z , tfinal ] = SimEoM4 ( M , T , Dd , fg , fj , Tef , Tcn , Dcn , qf , uf , 1 , rj );

% plot:
figure;
plot ( t , z );

[ tmp nq ] = size ( q );
[ tmp , tmp , ncnv ] = size ( Tcn );
ncn = ncnv * 3;
ns = nq + ncn;

sk = 1; % step skip
p = 1e-3; % pause time
qu = [ q , u ];

[ n , m ] = size ( rj ); % number of links
ss = length ( t ); % simulation steps
sa = floor ( ( ss - 1 ) / sk ) + 2; % animation steps
rcp = zeros ( n , 4 , sa );

parfor i = 1 : ss
    rcp(:,:,i) = subs (  rc(3,:) , qu , [ z( i , 1 : nq ) ; z( i , ns + 1 : ns + nq ) ] );
end
rcp(:,:,end) = subs (  rc(3,:) , qu , [ z( end , 1 : nq ) ; z( end , ns + 1 : ns + nq ) ] );
rcp = double ( rcp );

figure;
plot( rcp )
pause;

% animation
AnimEOM4 ( t , z , rj , qf , uf , nq , ncn );

matlabpool close
