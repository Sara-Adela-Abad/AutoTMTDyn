% Symbolic Derivation:

% properties:
syms lc1 lc2 lc3 m1 m2 m3 I1 I2 I3 l1 l2 k1 k2 k3 gx gy gz

%inputs:
lc = [ lc1 0 0 ; lc2 0 0; lc3 0 0 ];
m = [ m1 , m2 , m3 ];
I = sym ( zeros ( 3 , 3 , 3 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 );
j = sym ( zeros ( 1 , 5 , 3 ) );
j(:,:,1) = [ 2 inf 0 0 0 ]; j(:,:,2) = [ 2 inf l1 0 0 ]; j(:,:,3) = [ 2 inf l2 0 0 ];
jkd = sym (zeros ( 3 , 2 , 3 ) );
jkd(1,:,1) = [ k1 0 ]; jkd(1,:,2) = [ k2 0 ]; jkd(1,:,3) = [ k3 0 ];
g = [ gx , gy , gz ];

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

