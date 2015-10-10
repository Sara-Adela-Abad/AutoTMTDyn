% EOM Simulation:

clc
clear all
close all
pause( 1e-2 )


%% initialization:

sym l_sg f_f

l_c1 = 18e-2 / 2 ; l_c2 = 8e-2 / 2 ; l_cx3 = 10e-2 / 2 ; l_cz3 = 9.3e-2 / 2 ; l_c4 = 3.41e-2 / 2 ; l_c5 = 5.78e-2 / 2 ; l_cy6 = 0 ; l_cz6 = 1.3e-2 / 2 ;

m__b = 1e-2 ;
m_1 = m_b ; m_2 = m_b ; m_3 = m_b ; m_4 = m_b ; m_5 = m_b ; m_6 = m_b ;

I_b = 1e-6 ;
I_1 = I_b ; I_2 = I_b ; I_3 = I_b ; I_4 = I_b ; I_5 = I_b ; I_6 = I_b ;

q_b = - 75 * pi / 180 ;
l_1 = 18e-2 ; l_2 = 8e-2 ; l_x3 = 10e-2 ; l_z3 = 9.3e-2 ; l_4 = 3.41e-2 ; l_5 = 5.78e-2 ; l_y6 = 0 ; l_z6 = 1.3e-2 ;

cv_b = 1e-4 ;
k_s1 = 363.63 ; k_s5 = 300;
k_g = 1e3 ; c_vg = 1 ;


%% inputs:

g = [ 0 , 0 , -9.81 ]; % gravity vector

lc = [ 0 0 -l_c1 0 0 ; 0 0 -l_c2 0 0 ; l_cx3 0 -l_cz3 0 4 ; 0 0 -l_c4 0 0 ; 0 0 -l_c5 0 0 ; 0 l_cy6 -l_cz6 0 0 ;
	0 0 -l_se4 1 4 ; l_sfx4 0 -l_sfz4 1 5 ;
	0 -l_sfy2 -l_sfz23 1 6 ; 0 -l_sey23 -l_sez23 1 7 ;
	0 l_sfy3 -l_sfz23 1 6 ; 0 l_sey23 -l_sez23 1 7 ] ;

m = [ m_1 , m_2 , m_3 , m_4 , m_5 , m_6 ] ;

I = sym( zeros( 3, 3, 6 ) ) ;
I( : , : , 1 ) = I_1 * eye( 3 ) ; I( : , : , 2 ) = I_2 * eye( 3 ) ; I( : , : , 3 ) = I_3 * eye( 3 ) ; I( : , : , 4 ) = I_4 * eye( 3 ); I( : , : , 5 ) = I_5 * eye( 3 ) ; I( : , : , 6 ) = I_6 * eye( 3 ) ;

j = sym( zeros( 1 , 5 , 10 ) ) ; i = 1 ;
j( 1 , : , i ) = [ 0 0 0 0 0 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 1 q_b 0 0 -l_1 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 0 0 0 0 -l_2 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 0 0 0 0 inf ] ; i = i + 1 ;
j( 1 , : , i ) = [ 3 inf 0 l_x3 -l_z3 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 2 inf 0 0 -l_4 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 1 inf 0 0 -l_5 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 0 0 0 l_y6 -l_z6 ] ; i = i + 1 ;
j( 1 , : , i ) = [ 0 0 0 inf inf ] ; i = i + 1 ;

jkd = sym (zeros ( 3 , 2 , 6 ) ) ;
jkd( 1 , 1 , 1 ) = k_s1 ; jkd( 2 , 1 , 1 ) = c_vb ;
jkd( 2 , 1 , 2 ) = c_vb ;
jkd( 1 , 1 , 3 ) = k_s5 ; jkd( 2 , 1 , 3 ) = c_vb ;
jkd( 2 , 1 , 4 ) = c_vb ;
jkd( 3 , 1 , 5 ) = f_f ;
jkd( 1 , : , 6 ) = [ k_g l_sg ] ; jkd( 2 , 1 , 6 ) = c_vg ;


%% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

	
%% numerical simulation
[ t , z , tfinal ] = SimEoM5_S_v0_2( M , T , Dd , fg , fj , ref , Tef , qf , uf , rj , 1 );
plot ( t , z );
pause;


%% animation
AnimEOM( t , z , rj , qf , uf );