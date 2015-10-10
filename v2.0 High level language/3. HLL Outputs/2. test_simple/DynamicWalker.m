% EOM Simulation:
clc
clear all
close all

% Gravity vector
g = [0.0, 0.0, -9.81];

% Inputs

% Locations
lc = [
	% Position data from body link1 for a joint j1 relative to base
	0.0 0.0 0.5 0 0;
	% Position data from body link2 for a joint j2 relative to j1
	0.0 0.0 0.5 0 1;
	% Position data from body link3 for a joint j3 relative to j2
	0.0 0.0 0.5 0 2;
];

% Mass values
m = [
	% link1
	1.0,
	% link2
	1.0,
	% link3
	1.0];

% Inertia values
I = sym (zeros (3, 3, 3));
% Inertia for body link1
I (:, :, 1) = [
	1.0e-3 0.0 0.0;
	0.0 1.0e-3 0.0;
	0.0 0.0 1.0e-3];

% Inertia for body link2
I (:, :, 2) = [
	1.0e-3 0.0 0.0;
	0.0 1.0e-3 0.0;
	0.0 0.0 1.0e-3];

% Inertia for body link3
I (:, :, 3) = [
	1.0e-3 0.0 0.0;
	0.0 1.0e-3 0.0;
	0.0 0.0 1.0e-3];

% Joint specifications
j = sym (zeros (2, 5, 3));
% Joint rotations for joint j1
j (1, :, 1) = [
	2 inf 0 0 0];
% Joint rotations for joint j2
j (1:2, :, 2) = [
	0 0 0.0 0.0 1.0;
	1 inf 0 0 0];
% Joint rotations for joint j3
j (1:2, :, 3) = [
	0 0 0.0 1.0 1.0;
	1 inf 0 0 0];

% Stiffness data
jkd = sym (zeros (4, 2, 3));
% Stiffness values for joint j1
jkd (1, :, 1) = [ 5.0 5.0 ];
jkd (2, 1, 1) = 5.0;
jkd (3, 1, 1) = 5.0;
jkd (4, 1, 1) = 0;

% Stiffness values for joint j2
jkd (1, :, 2) = [ 5.0 5.0 ];
jkd (2, 1, 2) = 5.0;
jkd (3, 1, 2) = 5.0;
jkd (4, 1, 2) = 0;

% Stiffness values for joint j3
jkd (1, :, 3) = [ 5.0 5.0 ];
jkd (2, 1, 3) = 5.0;
jkd (3, 1, 3) = 5.0;
jkd (4, 1, 3) = 0;

% Run program -- Should this really be generated?

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
	TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );
