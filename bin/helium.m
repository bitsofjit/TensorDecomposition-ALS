%% Problem set up
clear;
nx = 100;
N = nx^3;
M = 10;
x = linspace(-M,M,nx);
h = x(2)-x(1);
[X, Y, Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);
e = ones(nx,1)/h^2;
L1 = spdiags([e -2*e e], -1:1, nx, nx);
I = speye(nx);
MLap = -1*kron(kron(L1,I),I) - kron(kron(I,L1),I) - kron(speye(nx^2),L1);

%% set up SCF iteration
max_its = 100;
Vext = -2./R(:);
Vtot = Vext;
Vh = Vtot;
Psi = Vtot;
for k = 1:max_its
    
    %eigenvalue problem solution
    opts.issym = 1;
    opts.isreal = 1;
    opts.v0 = Psi*h^(3/2);
    [Psi, E] = eigs(MLap/2+spdiags(Vtot,0,N,N),1,'SA',opts);
    % normalize for the integral
    Psi = Psi / h^(3/2);
    % 2 electrons
    
    %compute the parts of Vtot that depend on the density
    rho = 2*Psi.^2;
    Vxc = (-(3/pi)^(1/3)) *rho.^(1/3);
    [Vh, flag, rr, it, rv] = pcg(MLap,4*pi*rho,1e-6,200,[],[],Vh);
    Vtot = Vext + Vh + Vxc;
    
    % compute energies 
    K =  (h^3)*Psi'*(MLap*Psi);
    Eext = sum(rho.*Vext)*(h^3);
    Eh = sum(rho.*Vh)*(h^3)/2;
    Exc = (-3/4)*((3/pi)^(1/3))*sum(rho.^(4/3))*(h^3);
    Etot =  K + Eext + Eh + Exc;
    
    str1 = sprintf('Kinetic       | External       | Potential     | Exchange       | Total');
    str2 = sprintf('%d  | %d  | %d  | %d  | %d',K,Eext,Eh,Exc,Etot);
    disp(str1)
    disp(str2)
end
