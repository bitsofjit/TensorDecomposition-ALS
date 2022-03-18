%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [val,der] = sphericalharmonic(L,M,rvec)
% Purpose
% =======
% Get the value of the spherical harmonic Y_LM and its derivatives at xg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val,der] = sphericalharmonic(L,M,rvec)

if (L < 0 || L > 3) error('sphericalharmonic error: invalid L'); end
if (abs(M) > L) error('sphericalharmonic error: invalid M'); end

r = norm(rvec);
x = rvec(1); 
y = rvec(2); 
z = rvec(3); 

switch L
  case 0
    % Y_0^0
    val = 0.5/sqrt(pi); der = [0. 0. 0.];

  case 1
    c = sqrt(3./pi);
    r3 = 2.*r*r*r;
    switch M
      case -1
        % Y_1^{-1}
        val = 0.5*c*y/r;
        der = [ - c*x*y/r3  c*(r*r - y*y)/r3  -c*y*z/r3 ];
      case 0
        % Y_1^{0}  
        val = 0.5*c*z/r;
        der = [ - c*x*z/r3  -c*y*z/r3  c*(r*r - z*z)/r3 ];
      case 1
        % Y_1^{1}  
        val = 0.5*c*x/r;
        der = [ c*(r*r - x*x)/r3  -c*x*y/r3  -c*x*z/r3 ];
    end

  case 2
    c1 = sqrt(15./pi);
    c2 = sqrt(5./pi);
    r2 = r*r;
    switch M
      case -2
        val = c1*x*y/(2.*r2);
        der = [ c1*(r*r - 2*x*x)*y/(2.*r2*r2) c1*(r*r - 2*y*y)*x/(2.*r2*r2) ...
                -c1*x*y*z/(r2*r2) ];
      case -1
        val = c1*y*z/(2.*r2);
        der = [ -c1*x*y*z/(r2*r2) c1*(r*r - 2*y*y)*z/(2.*r2*r2) ...
                c1*(r*r - 2*z*z)*y/(2.*r2*r2) ];
      case 0
        val = c2*(-x*x - y*y + 2.*z*z)/(4.*r2);
        der = [ -3.*c2*x*z*z/(2.*r2*r2) -3.*c2*y*z*z/(2.*r2*r2) ...
                3.*c2*(r*r - z*z)*z/(2.*r2*r2) ];
      case 1
        val = c1*x*z/(2.*r2);
        der = [ c1*(r*r - 2.*x*x)*z/(2.*r2*r2) -c1*x*y*z/(r2*r2) ...
                c1*(r*r - 2.*z*z)*x/(2.*r2*r2) ];
      case 2
        val = c1*(x*x - y*y)/(4.*r2);
        der = [ c1*(r*r - x*x + y*y)*x/(2.*r2*r2) ...
                -c1*(r*r + x*x - y*y)*y/(2.*r2*r2) -c1*(x*x - y*y)*z/(2.*r2*r2) ];
    end

  case 3
    r3 = r*r*r;
    r5 = r*r*r*r*r;
    switch M
      case -3
        c = sqrt(35./(2.*pi));
        val = -c*y*(-3.*x*x + y*y)/(4.*r3);
        der = [ 3.*c*x*y*(2.*r*r - 3.*x*x + y*y)/(4.*r5) ...
                3.*c*(-3.*x*x*y*y + y*y*y*y + r*r*(x*x - y*y))/(4.*r5) ...
                3.*c*y*z*(-3.*x*x + y*y)/(4.*r5) ];
      case -2
        c = sqrt(105./pi);
        val = c*x*y*z/(2.*r3);
        der = [ c*(r*r - 3.*x*x)*y*z/(2.*r5) ...
                c*(r*r - 3.*y*y)*x*z/(2.*r5) ...
                c*(r*r - 3.*z*z)*x*y/(2.*r5) ];
      case -1
        c = sqrt(21./(2.*pi));
        val = -c*y*(x*x + y*y - 4.*z*z)/(4.*r3);
        der = [ c*(r*r - 15.*z*z)*x*y/(4.*r5) ...
                -c*(15.*y*y*z*z + r*r*(x*x - 4.*z*z))/(4.*r5) ...
                c*(11.*r*r - 15.*z*z)*y*z/(4.*r5) ];
      case 0
        c = sqrt(7./pi);
        val = c*z*(-3.*x*x - 3.*y*y + 2.*z*z)/(4.*r3);
        der = [ 3.*c*(r*r - 5.*z*z)*x*z/(4.*r5) ...
                3.*c*(r*r - 5.*z*z)*y*z/(4.*r5) ...
                -3.*c*(r*r*r*r - 6.*r*r*z*z + 5.*z*z*z*z)/(4.*r5) ];
      case 1
        c = sqrt(21./(2.*pi));
        val = -c*x*(x*x + y*y - 4.*z*z)/(4.*r3);
        der = [ -c*(15.*x*x*z*z + r*r*(y*y - 4.*z*z))/(4.*r5) ...
                c*(r*r - 15.*z*z)*x*y/(4.*r5) ...
                c*(11.*r*r - 15.*z*z)*x*z/(4.*r5) ];
      case 2
        c = sqrt(105./pi);
        val = c*z*(x*x - y*y)/(4.*r3);
        der = [ c*(2.*r*r - 3.*x*x + 3.*y*y)*x*z/(4.*r5) ...
                -c*(2.*r*r + 3.*x*x - 3.*y*y)*y*z/(4.*r5) ...
                c*(x*x - y*y)*(r*r - 3*z*z)/(4.*r5) ];
      case 3
        c = sqrt(35./(2.*pi));
        val = c*x*(x*x - 3.*y*y)/(4.*r3);
        der = [ 3.*c*(-x*x*x*x + 3.*x*x*y*y + r*r*(x*x - y*y))/(4.*r5) ...
                -3.*c*(2.*r*r + x*x - 3.*y*y)*x*y/(4.*r5) ...
                -3.*c*(x*x - 3.*y*y)*x*z/(4.*r5) ];
    end
end
