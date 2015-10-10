% EOM Simulation:
clc
clear all
close all

% properties:
m1 = 1;
springCoeff = 0;
springInit = 0;
dampViscous = 1e-2;
dampCoulomb = 0;
 
g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ 0 0 0.5 0 0 ;
    0 0 0.5 0 0 ;
    0 0 -0.5 0 0 ; 0 0 -0.5 0 0 ;
    0 0 0.5 0 2 ;
    0 0 -0.5 0 0 ;
    0 0 -0.5 0 5 ;
    0 0 0 1 1 ;
    0 0 -1.0 1 4 ;
    0 0 -1.0 2 4 ];

m = [ m1 , m1 , m1 , m1 , m1 , m1 , m1 ];

I = sym ( zeros ( 3 , 3 , 7 ) );
for i = 1:7
	I(:,:,i) = 1e-3 * eye ( 3 );
end

j = sym ( zeros ( 3 , 5 , 7 ) );
j(1,:,1) = [ 2 inf 0 0 0 ];
j(1,:,2) = [ 2 inf 0 0 1.0 ];
j(1,:,3) = [ 2 inf 1.0 0 0 ];
j(1,:,4) = [ 2 inf -1.0 0 0 ];
j(:,:,5) = [ 2 inf 0 0 1.0 ;
			 1 inf 0 0 0 ;
			 3 inf 0 0 0 ];
j(1,:,6) = [ 1 inf 0 0 1.0 ];
j(1,:,7) = [ 2 inf 0 0 1.0 ];

jkd = sym (zeros ( 3 , 2 , 9 ) );
for i = 1:9
	jkd(1,:,i) = [ springCoeff, springInit ];
    jkd(2,1,i) = dampViscous;
    jkd(3,1,i) = dampCoulomb;
end

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );

