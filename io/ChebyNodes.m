function Y = ChebyNodes(a, b, N)
% Computes the Chebyshev zero points (Y) in the domian -- (a, b)
% N = # of nodes

if nargin < 3
    N = 20;
    a = -1.0;
    b = 1.0;
end

index = linspace(1.0, N, N)';

% Scaling from (a, b) --> (-1, 1)
Y = (a + b)/2.0 + (b - a)/2.0*cos((2.0*index - 1.0)*pi/2.0/N);

% % Just checking:
% plot(Y, zeros(k,1), 'ro')

end           