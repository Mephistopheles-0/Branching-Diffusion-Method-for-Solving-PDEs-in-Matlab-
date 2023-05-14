function Branching_Matlab ()
% Parameters for the model
  T = 0.3; t0 = 0; x0 = 0; d = 100; m = d ;
  mu = zeros (d ,1); sigma = eye ( d )* sqrt (2);
  a = [0 2 0 -1] ;
  g = @ ( x ) 1./(1+0.2* norm ( x )^2)*1/2;

% Parameters for the algorithm
rng ( 'default' ); M = 10^7; beta = 1; p = [0 0.5 0 0.5] ;

% Branching method
tic ;
[mn, sd] = MC_BM (mu, sigma, beta, p, a, t0, x0, T, g, M);
runtime = toc ;

% Output
disp([' Terminal condition : u(T,x0) = ', num2str(g(x0)), '; ']);
disp([' Branching method : u(0,x0) ~ ', num2str(mn), '; ']);
disp([ ' Estimated standard deviation : ' num2str(sd) '; ' ]);
disp([ ' Estimated L2 - appr . error = ' num2str(sd/sqrt(M)) '; ' ]);
disp([ ' Elapsed runtime = ' num2str(runtime) '; ' ]);
end

function [ mn , sd ] = MC_BM ( mu , sigma , beta , p , a , t0 , x0 , T , g , M )
  mn = 0; sd = 0;
  for m =1: M
    result = BM_Eval ( mu , sigma , beta , p , a , t0 , x0 , T , g );
    mn = mn + result ;
    sd = sd + result ^2;
  end
  mn = mn / M ; sd = sqrt ( ( sd - mn ^2/ M )/ M );
end

function result = BM_Eval ( mu , sigma , beta , p , a , t0 , x0 , T , g )
  bp = BP ( mu , sigma , beta , p , t0 , x0 , T );
  result = 1;
  for k =1: size ( bp {1} ,2)
    result = result * g ( bp {1}(: , k ) );
  end
  if norm (a - p ) > 0
    for k =1: length ( a )
      if p ( k ) > 0
         result = result * ( a ( k )/ p ( k ) )^( bp {2}( k ) );
      elseif a ( k ) ~= 0
         error ( 'a ( k ) zero but p ( k ) non - zero ' );
      end
    end
  end
end

function bp = BP ( mu , sigma , beta , p , t0 , x0 , T )
  bp = cell (2 ,1);
  bp {2} = p *0;
  tau = exprnd (1/ beta );
  new_t0 = min ( tau + t0 , T );
  delta_t = new_t0 - t0 ;
  m = size ( sigma ,2);
  new_x0 = x0 + mu * delta_t + sigma * sqrt ( delta_t )* randn (m ,1);
  if tau >= T - t0
    bp {1} = new_x0 ;
  else
    [ tmp , nonlinearity ] = max ( mnrnd (1 , p ));
    bp {2}( nonlinearity ) = bp {2}( nonlinearity ) + 1;
    for k =1: nonlinearity -1
      tmp = BP ( mu , sigma , beta , p , new_t0 , new_x0 , T );
      bp{1} = horzcat(bp{1}, tmp{1});
      bp {2} = bp {2} + tmp {2};
    end
  end
end





