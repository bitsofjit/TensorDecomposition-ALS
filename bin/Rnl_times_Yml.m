function [enr_val, Yml] = Rnl_times_Yml(l, m, x0, x, alpha)
% R_{nl}(r)*Y^m_l(\theta,\phi) => Product of radial solution of
% Schrodinger equation and the spherical harmonics. The radial
% solution is taken to be the model pseduso-potential which
% mimics a Gaussian.
% 
% The Spherical harmonic (the quantum mechanics version)
% incorporating Condon?Shortley phase
% Works only for 3-dimensions
% Ref: specialFunctions.f90 of pufee package
% =====
% NOTE:
% =====
% real spherical harmonic function
% Real orthonormal set defined in terms of the complex spherical
% harmonics Y_l^m:
%            Ylm,                                   m = 0
%            val =   (-1)**m sqrt(2) Re(Ylm),       m > 0
%                    (-1)**|m| sqrt(2) Im(Yl|m|),   m < 0

if nargin < 5,
    l = 5;
    m = 3;
    x0 = [0.0; 0.0; 0.0];
    x = [0.5; 0.5; 0.5];
    alpha = 0.5;
end

if (abs(m) > l),
    error('m value outside range! Choose -l <= m <= l')
end

M = 0:l;

R = norm(x);    theta = acos(x(3)/R);   phi = atan(x(2)/x(1));

enr_val = zeros(l+1,1);

for i = 0:l,
    
    emm = M(i+1);
    
    if emm == 0,
        
        norm_const = (-1)^emm*sqrt((2*l + 1)/4.0/pi);
        
        enr_val(i+1) = norm_const;
        
    elseif emm < 0,
        
        emm = abs(emm);
        
        norm_const = (-1)^emm*sqrt(2.)*sqrt((2*l + 1)*factorial(l - emm)/ ...
                     4.0/pi/factorial(l + emm));
        
        enr_val(i+1) = norm_const*sin(emm*phi);
        
    elseif emm > 0,
        
        norm_const = (-1)^emm*sqrt(2.)*sqrt((2*l + 1)*factorial(l - emm)/ ...
                     4.0/pi/factorial(l + emm));
        
        enr_val(i+1) = norm_const*cos(emm*phi);
        
    end
    
end

% P = legendre(N,X) computes the associated Legendre functions 
% of degree N and order M = 0, 1, ..., N, evaluated for each element
% of X.  N must be a scalar integer and X must contain real values
% between -1 <= X <= 1. 

Yml = enr_val.*legendre(l,cos(theta));

if x0 == zeros(3,1);
    arg = x.^2;
else 
    arg = (x - x0).^2;
end

val_vec = exp(-alpha*sum(arg))*Yml;

m_out = M(m+1);

enr_val = val_vec(m_out+1);

Yml = Yml(m_out+1);

end