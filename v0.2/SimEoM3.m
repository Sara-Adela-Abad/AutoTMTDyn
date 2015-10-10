%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM3 ( M , T , D , fg , fj , Tef , Tcn , Dcn , q , u , dt )

global Mt;
global Tt;
global Dt;
global fgt;
global fjt;
global qt;
global ut;
global Teft;
global Tcnt;
global Dcnt;

Mt = M; Tt = T; Dt = D; fgt = fg; fjt = fj; Teft = Tef; Tcnt = Tcn; Dcnt = Dcn ; qt = q; ut = u;

[ tmp nq ] = size ( q );
[ tmp , tmp , ncnv ] = size ( Tcn );
m = nq + ncnv * 3; % # of states and constraint vectors (Lagrangian multiplyers)
par = [ nq ncnv*3 ];
z0 = 1e-6 * ones ( 1 , 2 * m );
t0 = 0;

% Standard ODE solver:
options = odeset (); %,'abstol',1*1e-6,'reltol',1*1e-6);
tspan = linspace ( t0 , t0 + dt , 500);
[ t , z , tfinal ] = ode113 ( @EOM , tspan , z0 , options , par );


function dz = EOM ( t , z , par )
t

global Mt;
global Tt;
global Dt;
global fgt;
global fjt;
global qt;
global ut;
global Teft;
global Tcnt;
global Dcnt;

Mt2 = Mt; Tt2 = Tt; Dt2 = Dt; fgt2 = fgt; fjt2 = fjt; Teft2 = Teft; Tcnt2 = Tcnt; Dcnt2 = Dcnt;

nq = par(1); ncn = par (2);
ns = nq + ncn;
qu = [ qt , ut ];
zq = [ z( 1 : nq , 1 ) ; z( ns + 1 : ns + nq , 1 ) ]'; % only states
u = z( ns + 1 : end );

value = {0*Mt; 0*Tt; 0*Dt; 0*fgt; 0*fjt; 0*Teft; 0*Tcnt; 0*Dcnt};
out = struct('f',value);

parfor i = 1 : 8

    switch i
        case 1
            Mt3 = subs (  Mt2 , qu , zq );
%             out(i).f = ( Mt3 );
            out(i).f = double ( Mt3 );
            
        case 2
            Tqt3 = subs ( Tt2 , qu , zq );
%             out(i).f = ( Tqt3 );
            out(i).f = double ( Tqt3 );
            
        case 3
            Dt3 = subs ( Dt2 , qu , zq );
%             out(i).f = ( Dt3 );
            out(i).f = double ( Dt3 );
            
        case 4
            fgt3 = subs ( fgt2 , qu , zq );
%             out(i).f = ( fgt3 );
            out(i).f = double ( fgt3 );
            
        case 5
            fjt3 = subs ( fjt2 , qu , zq );
%             out(i).f = ( fjt3 );
            out(i).f = double ( fjt3 );
            
        case 6
            Teft3 = subs ( Teft2 , qu , zq );
%             out(i).f = ( Teft3 );
            out(i).f = double ( Teft3 );
            
        case 7
            Tcnt3 = subs ( Tcnt2 , qu , zq );
%             out(i).f = ( Tcnt3 );
            out(i).f = double ( Tcnt3 );
            
        case 8
            Dcnt3 = subs ( Dcnt2 , qu , zq );
%             out(i).f = ( Dcnt3 );
            out(i).f = double ( Dcnt3 );
    end
    
end

t
% M = double ( out(1).f );
% Tq = double ( out(2).f );
% D = double ( out(3).f );
% fg = double ( out(4).f );
% fj = double ( out(5).f );
% Tef = double ( out(6).f );
% Tcn = double ( out(7).f );
% Dcn = double ( out(8).f );

M = out(1).f;
Tq = out(2).f;
D = out(3).f;
fg = out(4).f;
fj = out(5).f;
Tef = out(6).f;
Tcn = out(7).f;
Dcn = out(8).f;

A = zeros ( ns , ns );
A ( 1 : nq , : ) = [ Tq.'*M*Tq -Tcn(:,:,1).' -Tcn(:,:,2).' ];
A ( nq + 1 : nq + 3, 1 : nq ) = Tcn(:,:,1);
A ( nq + 4 : nq + 6 , 1 : nq ) = Tcn(:,:,2);

B = zeros ( ns , 1 );
B ( 1 : nq , 1 ) = Tq.' * ( fg - M * D * u ( 1 : nq ) ) + fj;
B ( nq + 1 : nq + 3 , 1 ) = - Dcn(:,:,1) * u ( 1 : nq );
B ( nq + 4 : nq + 6 , 1 ) = - Dcn(:,:,2) * u ( 1 : nq );

dzt = A \ B;
dz = [ u ; dzt ];

