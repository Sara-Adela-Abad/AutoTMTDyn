% EOM Simulation:
clc
clear all
close all

% properties:
m1 = 1;
springCoeff = 0;
springInit = 0;
dampViscous = 0;
dampCoulomb = 10*0.1;
r = 20e-3;
phi_s = pi/2;
lt = [ r 0 -r ]; % [ r*sin(phi_s) 0 -r*(1-cos(phi_s)];

g = [ 0 , 0 , -9.81*0]; % gravity vector

%inputs:
lc = ones(6, 1)*[ lt/2 0 0 ];
lc (1,:) = 2*[ lt/2 0 0 ];
lc (4,:) = 2*[ lt/2 0 0 ];

m = m1*ones(1, 6);
m(1) = 2*m1;
m(4) = 2*m1;

I = sym ( zeros ( 3 , 3 , 6 ) ); % number of matrices in the cube (here 7) equal to number of bodies
for i = 1:6
	I(:,:,i) = 1e-3 * eye ( 3 );
end
	I(:,:,1) = 2*1e-3 * eye ( 3 );
	I(:,:,4) = 2*1e-3 * eye ( 3 );

j = sym ( zeros ( 2 , 5 , 6 ) );
j(1,:,1) = [ 1 inf 0 0 r ];
for i = 2:6
	j(:,:,i) = [ 2 phi_s lt ;
        1 inf 0 0 0 ];
end
	j(:,:,2) = [ 2 2*phi_s 0 0 -2*r ;
        1 inf 0 0 0 ];
	j(:,:,5) = [ 2 2*phi_s 0 0 -2*r ;
        1 inf 0 0 0 ];

jkd = sym (zeros ( 3 , 2 , 6 ) );	
jkd(1,:,i) = [ springCoeff, springInit ];
jkd(2,1,i) = dampViscous;
jkd(3,1,i) = dampCoulomb; 
for i = 2:6
    jkd (4,1,i) = 1;
end

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 3 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );

