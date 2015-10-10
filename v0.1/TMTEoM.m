%% TMT EOM Derivator:
% ===========================================
% [ M , T , Dd , fg , fj , rj , rc , vc , wc , qf , uf ] = ...
%    TMTEoM ( lc , m , I , j , jkd , g );
%
% ===========================================
% Author:
%   S.M.Hadi Sadati
%   PhD student - King's College London
%   2014
%
% ===========================================
% Help:
%
% Inputs:
% ( n ) Number of links (every thing that has mass and/or rotational inertia
% ( nq ) Number of generalized coordinates
% lc : nx3 vector of link COM positions in link frame
%     Joint frame origin attached to joint
% m : nx1 vector of link masses
% I : 3x3xn cube of link inertia matrices w.r.t. link frame
% j : rx5xn cube of each joint (translation and rotation leading to each joint)
%     translation from previous joint position:
%       ( r , 1:3 , n ) = linear transformation value along previous link frame x-y-z axes 
%     and rotation of new link frame w.r.t. previous link frame
%       ( r , 4 , n ) = 1 : 3 axis indicator number for Euler post-multiplication rotations about sequensive temporary frames' x-y-z axes ),
%     and the rotation values:
%       ( r , 5 , n ) = value.
%     PLACE "inf" FOR "VALUE" IF IT IS ONE OF THE GENERALIZED COORDINATES!!!
% jkd : 3x2xnq cube generalized coordinates 3 pair of 
%     [ spring coeff.s ; viscous damping coeff.s ; external input ]
%     and their initial pos.s ( Only spring coeff. needs it! )
% g : 1x3 gravity acceleration vector
% 
% Outputs:
% M : Mass matrix
% T : Transformation matrix
% D : Damping/Stiffness matrix
% fg : Gravity force virtual work
% fj : Joint stiffness/damping virtual work acting directly on generalized coordinates
% rj : Joint absolute position vector in base frame
% rc : Joint COM absolute position vector in base frame
% vc : Joint COM absolute linear velocity vector in base frame
% wc : Joint COM absolute rotational velocity vector in link frame
% qf : Generalized coordinates
% uf : Generalized coordinates derivatives
% 
% ============================================
% Example:
%
% syms lc1 lc2 lc3 m1 m2 m3 I1 I2 I3 l1 l2 k1 k2 k3 gx gy gz
% lc = [ lc1 0 0 ; lc2 0 0; lc3 0 0 ];
% m = [ m1 , m2 , m3 ];
% I = sym ( zeros ( 3 , 3 , 3 ) );
% I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 );
% j = sym ( zeros ( 1 , 5 , 3 ) );
% j(:,:,1) = [ 2 inf 0 0 0 ]; j(:,:,2) = [ 2 inf l1 0 0 ]; j(:,:,3) = [ 2 inf l2 0 0 ];
% jkd = sym (zeros ( 3 , 2 , 3 ) );
% jkd(1,:,1) = [ k1 0 ]; jkd(1,:,2) = [ k2 0 ]; jkd(1,:,3) = [ k3 0 ];
% g = [ gx , gy , gz ];
%
% ============================================



%% Main:
function [ M , T , Dd , fg , fj , rj , rc , vc , wc , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g )


% >> Initialization:

j = sym ( j );
lc = sym ( lc );
m = sym ( m );
I = sym ( I );
j = sym ( j );
jkd = sym ( jkd );
g = sym ( g );

[ n , tmp ] = size ( lc ); % number of links
[ rl , tmp , tmp ] = size ( j ); % maximum number of rotation in each joint
nq = 0; % number of generalized coordinates
for i = 1 : n
    for i1 = 1 : rl
        for i2 = 2 : 5
        if j( i1 , i2 , i ) == inf
            nq = nq + 1;
        end; end; end; end

M = sym ( zeros ( 6 * n , 6 * n ) );
T = sym ( zeros ( 6 * n , nq ) );
Dd = sym ( zeros ( 6 * n , nq ) );
fg = sym ( zeros ( 6 * n , 1 ) );
fj = sym ( zeros ( nq , 1 ) );
rj = sym ( zeros ( n , 3 ) );
rc = sym ( zeros ( n , 3 ) );
vc = sym ( zeros ( n , 3 ) );
wcb = sym ( zeros ( n , 3 ) );
wc = sym ( zeros ( n , 3 ) );

syms q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20
syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
q = [ q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 ];
u = [ u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 ];
iq = 0;
qf = sym ( [] );
uf = sym ( [] );

R = sym ( zeros ( 3 , 3 , n ) );
TR = sym ( zeros ( 4 , 4 , n ) );
Rt = sym ( zeros ( 3 , 3 , n ) );
TRt = sym ( zeros ( 4 , 4 , n ) );
Wc = sym ( zeros ( 3 , 3 , n ) );
fgt = sym ( [] );


% >> Derivation:
    
for i = 1 : n
    
    im = 3 * i;
    iI = 3 * n + 3 * i;
    
    M ( im-2 : im , im-2 : im ) = m(i) * eye ( 3 );
    M ( iI-2 : iI , iI-2 : iI ) = I(:,:,i);
    
    R(:,:,i) = sym ( eye ( 3 ) );
    TR(:,:,i) = sym ( eye ( 4 ) );
    
    for i1 = 1 : rl
        
        if j(i1,:,i) == 0
            break; end
        
        for i2 = 2 : 5            
            if j(i1,i2,i) == inf
                
                iq = iq + 1;
                
                j(i1,i2,i) = q(iq);
                qf = [ qf , q(iq) ];
                uf = [ uf , u(iq) ];
                
                fj(iq,1) = ...
                    - jkd(1,1,iq) * ( q(iq) - jkd(1,2,iq) ) ... % spring
                    - jkd(2,1,iq) * u(iq) ... % viscous
                    + jkd(3,1,iq); % external input
                                
            end; end
        
        [ tmp1 , tmp2 ] = TRm ( j(i1,1:5,i) );
        R(:,:,i) = R(:,:,i) * tmp1;
        TR(:,:,i) = TR(:,:,i) * tmp2;
        
                
    end
    
    if i > 1
        Rt(:,:,i) = Rt(:,:,i-1) *  R(:,:,i);
        TRt(:,:,i) = TRt(:,:,i-1) *  TR(:,:,i);
    else
        Rt(:,:,i) = R(:,:,i);
        TRt(:,:,i) = TR(:,:,i);
    end
    
%     R(:,:,i) = fclose ( R(:,:,i) );
%     TR(:,:,i) = simple ( TR(:,:,i) );
%     Rt(:,:,i) = simple ( Rt(:,:,i) );
%     TRt(:,:,i) = simple ( TRt(:,:,i) );
    
    tmp3 = sym ( [] );
    for i1 = 1 : 3
        tmp4 = jacobian ( Rt(:,i1,i) , qf ) * uf.';
        tmp3 = [ tmp3 , tmp4 ];
    end
    Wc(:,:,i) = tmp3 * Rt(:,:,i).'; % angular velocity tensor
    wcb(i,:) = [ Wc(3,2,i) , Wc(1,3,i) , Wc(2,1,i) ]; % w in base frame
    wc(i,:) = ( Rt(:,:,i).' * wcb(i,:).' ).'; % w in link frame
    
%     wc(i,:) = simple ( wc(i,:) );
    
    for i1 = 1 : 3
        for i2 = 1 : iq
            [ tmp5 , tmp6 ] = coeffs ( wc(i,i1) , uf(i2) );
            if isempty ( tmp5 )
                tmp5 = sym ( 0 );
            end
            T ( iI - 3 + i1 , i2 ) = tmp5(1);
%             T ( iI - 3 + i1 , i2 ) = simple ( T ( iI - 3 + i1 , i2 ) );
        end; end
    
    Dd ( iI - 2 : iI , 1 : iq ) = jacobian ( wc(i,:).' , qf );
%     D ( iI - 2 : iI , 1 : iq ) = simple ( D ( iI - 2 : iI , 1 : iq ) );
    
    rj(i,:) = TRt(1:3,4,i).';
    tmp7 = TRt(:,:,i) * [ lc(i,:) , 1 ].';
    rc(i,:) = tmp7(1:3).';
    tmp8 = jacobian ( rc(i,:).' , qf );
    vc(i,:) = ( tmp8 * uf.' ).';
    
%     rc(i,:) = simple ( rc(i,:) );
%     vc(i,:) = simple ( vc(i,:) );    
    
    T ( im - 2 : im , 1 : iq ) = tmp8;
    Dd ( im - 2 : im , 1 : iq ) = jacobian ( vc(i,:).' , qf );
     
%     T ( im - 2 : im , 1 : iq ) = simple ( T ( im - 2 : im , 1 : iq ) );
%     D ( im - 2 : im , 1 : iq ) = simple ( D ( im - 2 : im , 1 : iq ) );
    
    fgt = [ fgt , g ];
    
end

rj( i + 1 ,:) = rc(i,:);
fg ( 1 : 3 * n , 1 ) = M ( 1 : 3 * n , 1 : 3 * n ) * fgt.';
    
M = simple( M );
T = simple( T );
Dd = simple( Dd );
fg = simple( fg );
fj = simple( fj );
rj = simple( rj );
rc = simple( rc );
vc = simple( vc );
wc = simple( wc );

matlabFunction ( M , 'file' , 'M.m' ); % save as matlab function
matlabFunction ( T , 'file' , 'T.m' );
matlabFunction ( Dd , 'file' , 'Dd.m' );
matlabFunction ( fg , 'file' , 'fg.m' );
matlabFunction ( fj , 'file' , 'fj.m' );
matlabFunction ( rj , 'file' , 'rj.m' );
matlabFunction ( rc , 'file' , 'rc.m' );
matlabFunction ( vc , 'file' , 'vc.m' );
matlabFunction ( wc , 'file' , 'wc.m' );

% ccode ( M , 'file' , 'M.cpp' ); % save as c file
% ccode ( T , 'file' , 'T.cpp' );
% ccode ( Dd , 'file' , 'Dd.cpp' );
% ccode ( fg , 'file' , 'fg.cpp' );
% ccode ( fj , 'file' , 'fj.cpp' );
% ccode ( rj , 'file' , 'rj.cpp' );
% ccode ( rc , 'file' , 'rc.cpp' );
% ccode ( vc , 'file' , 'vc.cpp' );
% ccode ( wc , 'file' , 'wc.cpp' );

M_c = ccode ( M );  % save as c format text file
T_c = ccode ( T );
Dd_c = ccode ( Dd );
fg_c = ccode ( fg );
fj_c = ccode ( fj );
rj_c = ccode ( rj );
rc_c = ccode ( rc );
vc_c = ccode ( vc );
wc_c = ccode ( wc );

M_o = fopen ( 'M.cpp' , 'w' );
T_o = fopen ( 'T.cpp' , 'w' );
Dd_o = fopen ( 'Dd.cpp' , 'w' );
fg_o = fopen ( 'fg.cpp' , 'w' );
fj_o = fopen ( 'fj.cpp' , 'w' );
rj_o = fopen ( 'rj.cpp' , 'w' );
rc_o = fopen ( 'rc.cpp' , 'w' );
vc_o = fopen ( 'vc.cpp' , 'w' );
wc_o = fopen ( 'wc.cpp' , 'w' );

fwrite ( M_o , M_c );
fwrite ( T_o , T_c );
fwrite ( Dd_o , Dd_c );
fwrite ( fg_o , fg_c );
fwrite ( fj_o , fj_c );
fwrite ( rj_o , rj_c );
fwrite ( rc_o , rc_c );
fwrite ( vc_o , vc_c );
fwrite ( wc_o , wc_c );

fclose ( M_o );
fclose ( T_o );
fclose ( Dd_o );
fclose ( fg_o );
fclose ( fj_o );
fclose ( rj_o );
fclose ( rc_o );
fclose ( vc_o );
fclose ( wc_o );


%% Complementary Functions:

function [ R , TR ] = TRm ( r ) % base rotation matices

i = r(1);
x = r(2);

switch i
    
    case 1
        R = [ ...
            1       0        0 ;
            0  cos(x)  -sin(x);
            0  sin(x)   cos(x)];
        
    case 2
        R = [ ...
             cos(x)  0  sin(x);
                  0  1       0;
            -sin(x)  0  cos(x)];
        
    case 3
        R = [ ...
            cos(x)  -sin(x)  0;
            sin(x)   cos(x)  0;
                 0        0  1];
    otherwise
        R = sym ( eye ( 3 ) );
        
end

TR = sym ( [ [ R , r(3:5).' ] ; 0 0 0 1 ] );
 

%% Notes:
% Each TR contains translation and then a rotation

