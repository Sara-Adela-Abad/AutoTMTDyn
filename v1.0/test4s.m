% EOM Simulation:
clc
close all
pause(1e-2);

%% properties:
l = 1; % Link lenght
ml = 1; % mass
Il = 1e-3; % rotational inertia
rb01 = [ 1 0 0 ]; rb02 = [ -l*sin(pi/6) l*cos(pi/6) 0 ]; rb03 = [ -l*sin(pi/6) -l*cos(pi/6) 0 ]; % joint pos. vec. in base from base centre
lcp = 1; rp12 = [ -l-l*sin(pi/6) l*cos(pi/6) 0 ]; rp13 = [ -l-l*sin(pi/6) -l*cos(pi/6) 0 ]; % joint pos. vec. in platform from joint II
g = [ 0 , 0 , 1 ]; % gravity vector
% pi = 3.1415;

%inputs:
lc = [ -lcp 0 0 0 0 ; 0 0 0 2 0 ; 0 0 0 2 1 ];
m = ml  * ones(1,7);
I = sym ( zeros ( 3 , 3 , 1 ) );
I(:,:,1) = Il * eye ( 3 );

% j = sym ( zeros ( 4 , 5 , 6 ) );
% j(1,:,1) = [ 1 inf rb1 ]; j(2,:,1) = [ 3 inf 0 0 0 ]; j(3,:,1) = [ 2 inf 0 0 0 ]; % x-z-y
% j(1,:,2) = [ 2 inf 0 0 l1 ]; j(2,:,2) = [ 3 inf 0 0 0 ]; j(3,:,2) = [ 1 inf 0 0 0 ]; % y-z-x
% j(1,:,3) = [ 3 inf rp2 ]; j(2,:,3) = [ 1 inf 0 0 0 ]; j(3,:,3) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,4) = [ 3 inf rp3 ]; j(2,:,4) = [ 1 inf 0 0 0 ]; j(3,:,4) = [ 2 inf 0 0 0 ]; % z-x-y
% j(1,:,5) = [ 2 inf 0 0 -l2 ]; j(2,:,5) = [ 3 inf 0 0 0 ]; j(3,:,5) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,5) = [ 0 0 rb2 ];
% j(1,:,6) = [ 2 inf 0 0 -l3 ]; j(2,:,6) = [ 3 inf 0 0 0 ]; j(3,:,6) = [ 1 inf 0 0 0 ]; % y-z-x
% j(4,:,6) = [ 0 0 rb3 ];

j = sym ( zeros ( 6 , 5 , 3 ) );
j(1,:,1) = [ 2 inf rb01 ]; j(2,:,1) = [ 1 inf l 0 0 ]; j(3,:,1) = [ 2 inf 0 0 0 ]; j(4,:,1) = [ 1 inf 0 0 0 ]; j(5,:,1) = [ 2 inf -l 0 0 ]; % t-y-t-x-y-x-t-y
j(1,:,2) = [ 3 2*pi/3 rp12 ]; j(2,:,2) = [ 2 inf 0 0 0 ]; j(3,:,2) = [ 1 inf l 0 0 ]; j(4,:,2) = [ 2 inf 0 0 0 ]; j(5,:,2) = [ 1 inf 0 0 0 ]; j(6,:,2) = [ 2 inf -l 0 0 ]; % t-z-y-t-x-y-x-t-y
j(1,:,3) = [ 3 -2*pi/3 rp13 ]; j(2,:,3) = [ 2 inf 0 0 0 ]; j(3,:,3) = [ 1 inf l 0 0 ]; j(4,:,3) = [ 2 inf 0 0 0 ]; j(5,:,3) = [ 1 inf 0 0 0 ]; j(6,:,3) = [ 2 inf -l 0 0 ]; % t-z-y-t-x-y-x-t-y

jkd = sym (zeros ( 3 , 2 , 15 ) );

%% Derivation:

% matlabpool open local 2

% EOM:
[ M1 , T1 , Dd1 , fg1 , fj1 , rj1 , rc1 , vc1 , wc1 , ref1 , rcn1 ,  Tef1 , Tcn1 , Dcn1 , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );
 
% numerical simulation
% [ t , z , tfinal ] = SimEoM4 ( M1 , T1 , Dd1 , fg1 , fj1 , Tef1 , Tcn1 , Dcn1 , qf , uf , 5 , rj1 );
[ t , z , tfinal ] = SimEoM4s ( Tcn1 , qf , uf , 1 );

% plot:
[ tmp nq ] = size ( qf );
[ tmp , tmp , ncnv ] = size ( Tcn1 );
ncn = ncnv * 3;
ns = nq + ncn;

figure;
plot ( t , z(:,1:nq) );
pause(1e-1)

figure;
plot ( t , z(:,ns+1:end) );
pause(1e-1)


sk = 1; % step skip
p = 1e-3; % pause time
qu = [ qf , uf ];

[ n , m ] = size ( rj1 ); % number of links
ss = length ( t ); % simulation steps
sa = floor ( ( ss - 1 ) / sk ) + 2; % animation steps
rcp = zeros ( sa , 3 );

for i = 1 : ss
    rcp(i,:) = subs (  rc1(1,:) , qu , [ z( i , 1 : nq ) ; z( i , ns + 1 : ns + nq ) ] );
end
rcp(end,:) = subs (  rc1(1,:) , qu , [ z( end , 1 : nq ) ; z( end , ns + 1 : ns + nq ) ] );
rcp = double ( rcp );

figure;
plot( t , rcp(1:end-1,:) )
pause(1e-1);

% % animation
% AnimEOM4 ( t , z , rj , qf , uf , nq , ncn );

% matlabpool close
