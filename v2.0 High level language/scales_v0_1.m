% EOM Simulation:
clc
clear all
close all

% properties:
m1 = 1e-1;
springCoeff = 0;
springInit = 0;
dampViscous = 0;
dampCoulomb = 0*1e-3;
r = 20e-3;
phi_s = pi/2;
lt = [ r 0 -r ]; % [ r*sin(phi_s) 0 -r*(1-cos(phi_s)];

g = [ 0 , 9.81*1, 0]; % gravity vector

%inputs:
lc = ones (8, 1)*[ lt/2 0 0 ];

m = m1*ones (1, 8);

I = sym ( zeros ( 3 , 3 , 8 ) ); % number of matrices in the cube (here 7) equal to number of bodies
for i = 1:8
	I(:,:,i) = 1e-3 * eye ( 3 );
end

j = sym ( zeros ( 2 , 5 , 8 ) );
j(1,:,1) = [ 1 inf 0 0 r ];
for i = 2:8
	j(:,:,i) = [ 2 phi_s lt ;
        1 inf 0 0 0 ];
end

jkd = sym (zeros ( 4 , 2 , 8 ) ); 
jkd(1,:,i) = [ springCoeff, springInit ];
jkd(2,1,i) = dampViscous;
jkd(3,1,i) = dampCoulomb;
for i = 2:8
    jkd (4,1,i) = 1;
end

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 0.15 );
plot ( t , z );
pause;

% animation
figure;
AnimEOM ( t , z , rj , qf , uf );

