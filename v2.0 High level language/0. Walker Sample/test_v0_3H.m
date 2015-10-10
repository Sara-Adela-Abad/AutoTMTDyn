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
lc = [ 0 0 0.5 0 0 ; // number of rows equal to number of important points (number of required relative frames): number of the joints + number of constraints locations + number of external force locations
    0 0 0.5 0 0 ;
    0 0 -0.5 0 0 ; 0 0 -0.5 0 0 ;
    0 0 0.5 0 2 ;
    0 0 -0.5 0 0 ;
    0 0 -0.5 0 5 ;
    0 0 0 1 1 ;
    0 0 -1.0 1 4 ;
    0 0 0 2 4 ];

m = [ m1 , m1 , m1 , m1 , m1 , m1 , m1 ]; // number of elements equal to number of bodies

I = sym ( zeros ( 3 , 3 , 7 ) ); // number of matrices in the cube (here 7) equal to number of bodies
for i = 1:7
	I(:,:,i) = 1e-3 * eye ( 3 );
end

j = sym ( zeros ( 3 , 5 , 10 ) ); // number of matrices in the cube (here 10) equal to number of lc rows (number of requiered relative frames)
j(1,:,1) = [ 2 inf 0 0 0 ];
j(1,:,2) = [ 2 inf 0 0 1.0 ];
j(1,:,3) = [ 2 inf 0 1.0 1.0 ];
j(1,:,4) = [ 2 inf 0 0 -1.0 ];
j(:,:,5) = [ 2 inf 0 0.5 1.0 ;
			 1 inf 0 0 0 ;
			 3 inf 0 0 0 ];
j(1,:,6) = [ 1 inf 0 -0.25 1.0 ];
j(1,:,7) = [ 1 inf 0 0.25 1.0 ];
j(1,:,8) = [ 0 0 0 0 0 ]; // external load, not neccesary (it was zero)
j(1,:,9) = [ 0 0 0 0 0 ]; // external load
j(1,:,10) = [ 2 inf 0 0 -1.0 ]; // constraint joint

jkd = sym (zeros ( 3 , 2 , 10 ) ); // number of matrices in the cube (here 10) equal to number of states (degrees of freedom)
for i = 1:10
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

