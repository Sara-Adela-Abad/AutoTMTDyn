%% Animation:
% ====================================================

function AnimEOM ( t , z , rj , q , u )

% Data:

sk = 1; % step skip
p = 1e-3; % pause time
qu = [ q , u ];

[ n , m ] = size ( rj );
ss = length ( t ); % simulation steps
sa = floor ( ( ss - 1 ) / sk ) + 2; % animation steps
rjt = zeros ( n , 3 , sa );

for i = 1 : sk : ss
    rjt(:,:,i) = subs (  rj , qu , z(i,:) );
end
rjt(:,:,end) = subs (  rj , qu , z(end,:) );
rjt = double ( rjt );


% Windows set:

clf

window_xmin = min ( min ( rjt(:,1,:) ) );
window_xmax = max ( max ( rjt(:,1,:) ) );
window_ymin = min ( min ( rjt(:,2,:) ) );
window_ymax = max ( max ( rjt(:,2,:) ) );
window_zmin = min ( min ( rjt(:,3,:) ) );
window_zmax = max ( max ( rjt(:,3,:) ) );

ad = max( abs ( [ window_xmin window_xmax window_ymin window_ymax window_zmin window_zmax ] ) );

axis('equal')
axis ( [ window_xmin window_xmax window_ymin window_ymax window_zmin window_zmax ] ...
    + [ -ad ad -ad ad -ad ad ] );

set(gcf,'Color',[1,1,1])


% Animation:

clr = { 'blue' , 'red' , 'green' , 'magenta' , 'cyan' , 'yellow' , 'black' };

for i = 1 : n - 1
    
    xs = [ rjt(i,1,1) , rjt(i+1,1,1) ];
    ys = [ rjt(i,2,1) , rjt(i+1,2,1) ];
    zs = [ rjt(i,3,1) , rjt(i+1,3,1) ];
    
    ln(i) = line ( xs , ys , zs , ...
        'linewidth' , 2 , 'erase' , 'xor' ,'color' , clr{i} );

end

for i1 = 1 : sk : ss
    
    title( t(i1) );
    
    for i = 1 : n - 1
        
        xs = [ rjt(i,1,i1) , rjt(i+1,1,i1) ];
        ys = [ rjt(i,2,i1) , rjt(i+1,2,i1) ];
        zs = [ rjt(i,3,i1) , rjt(i+1,3,i1) ];
                
        set( ln(i) , 'XData' , xs , 'YData' , ys , 'ZData' , zs );

    end
            
    drawnow 
    pause( p );

end