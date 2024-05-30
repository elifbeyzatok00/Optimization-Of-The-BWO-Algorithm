function [bestSolution, bestFitness, iteration]=besd(fhd, dimension, maxIteration, fNumber)

config;

N=30;
D=dimension;
low=lbArray;
up=ubArray;
maxcycle=ceil(maxIteration/N);
% Initialization
K = 5 ;   %  LINE 2

S=zeros(K*N,D);
for i=1:K*N
    S(i,:) = rand( 1 , D ) .* (( up - low ) + low) ;   
end
fitS = testFunction(S', fhd, fNumber);
[gmin, ind] = min( fitS ) ;                         %  LINE 4
gbest =S ( ind, : ) ;                                      %  LINE 4

% Iterative search phase
for epk = 1 : maxcycle    %  LINE 5
   
   %  Generation of j0, j1, and j2
   while 1,  j1 = randperm( K * N ) ;     j2 = randperm( K*N ) ;    if sum( j1 == j2 ) == 0; break ; end,    end      %  LINE 7
   j1 = j1( 1 : N ) ;       j2 = j2( 1 : N )  ;        %  LINE 8
   j0 = j1 ;                                                          %  LINE 9
   
   %  Setting sub-pattern matrix P and fitP
   P = S( j1 , : ) ;   fitP = fitS( j1 ) ;     %  LINE 11
   
   % Generation of Bijective Vectors ; dv1
   dv1 = S( j2 , : ) ;         %  LINE 13
   
   % Top-N-Best pattern vectors
   [~ , index ] = sort( fitS , 'ascend' ) ;    H = S( index , : ) ;    %     LINE 16
   tbest = P ;   %  pre-memory
   for i = 1 : N ,  tbest( i  ,  :  ) = H( ceil( rand^3 * K * N ) , :  ) ; end        %  LINE 16
   
   % Generation of Bezier Mutation Vectors ; dv2        %  LINES 18 - 19
   while 1
      j1 = randperm( N ) ;  j2 = randperm( N ) ;
      if sum( j1 == 1:N ) == 0  &&  sum( j2 == 1:N ) == 0  && sum( j1 == j2 ) == 0; break ; end
   end
   dv2 = P ;  %  pre-memory
   for i =1 : N
      X = [ P(i , : )  ;  P( j1(i)  , : ) ;  P( j2(i)  ,  : ) ;  tbest(i , : )  ] ;
      X = X( randperm( 4 ) , : )  ;
      B = bernsteinMatrix( 3 , rand) ;
      dv2( i , : ) = B * X ;
   end
   
   % Evolutionary Step Size
   a = randn( N , 1 ) ;  b = 1 + rand( N , 1 ) ;   c = randn( N , 1 ) .^ randi( 7 , N , 1 ) ;   F = a .* b .^ c ;      %  LINE 21
   % F = random( 'stable' , 1 , 0 , 1 , 0 , N , 1 )  ;      %   Users can use this line instead of line 59 for CAUCHY-SCALING OF BeSD ; cauchy(mu=0,std=1)
   % F = random('stable', 0.50 , 1 , 1 , 0 , N , 1 )  ;   %  Users can use this line instead of line 59 for LEVY-SCALING OF BeSD
   
   
   % Generation of Crossover Control Matrix ; map
   [ map1 , map2 ] = genmap( N , D ) ;          %  LINES 23 - 26
   if rand < rand, map = map1 ; else map = map2 ; end    %  LINE 27
   
   % Generation of Trial Pattern Vectors ; T
   w1 = randn(N,1) ;  w2 = randn(N,1) ;                          %  LINE 29

   for iii = 1 : N
       T(iii,:) = P(iii,:) + map(iii,:) .* F(iii) .* ( w1(iii) .* ( dv2(iii,:) - P(iii,:) ) + w2(iii) .* ( dv1(iii,:) - P(iii,:) ) ) ;
   end
   %  LINE 30
   
   % Boundary Control Mechanism   ;    %  LINE 32
   for i = 1 : N
      for index = 1 : D
         if T(i,index) < low, T(i,index) = rand * ( up(index) - low(index) ) + low(index) ; end
         if T(i,index) > up,  T(i,index) = rand * ( up(index) - low(index) ) + low(index) ; end
      end
   end
   
   %  Update the sub-pattern matrix, P and fitP  ;      %  LINES 34 - 35
   fitT = testFunction(T', fhd, fNumber);
   ind = fitT < fitP ;
   P( ind , :  ) = T( ind , : ) ;
   fitP( ind )  = fitT( ind ) ;
   
   %  Update the Global Solution, gbest
   [BestVal, index ] = min( fitP ) ;            %  LINE 37
   BestP = P( index , : ) ;
   if BestVal < gmin ,  gmin = BestVal ;   gbest = BestP ;   end        %   LINE 38
   
   % Update the pattern matrix, S and fitS
   S( j0 , : ) = P ;    fitS( j0 ) = fitP ;   %  LINE 40
   
   % Report the objective function value of global solution, gmin     ;      %  LINE 42 - 43
   
end  % for epk
bestSolution=gbest;
bestFitness=gmin;
iteration=1;
end

%   Generation of Crossover Control Matrix; map    ;  LINES 23 - 27
function [ map1 , map2 ] = genmap( N , D )
map1 = zeros( N, D ) ;  % LINE 23
map2 = zeros( N, D ) ;   % LINE 23
for j = 1 :  N
   h = randperm( D ) ; w = rand .^ randi( 7 )  ;           %  LINE 25
   map1( j , h( 1 : ceil( w * D ) ) ) = 1;                           %  LINE 25
   h = randperm( D ) ; w = 1 - rand .^ randi( 7 ) ;       %  LINE 26
   map2( j , h( 1 : ceil( w * D ) ) ) = 1 ;                           %  LINE 26
end
end