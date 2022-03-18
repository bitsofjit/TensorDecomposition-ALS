function [X,Y,Z, eVecs, eVals] = hydrogen3dFD(n)
% 3d finite difference hydrogen atom
% INPUT
% n = number of grid points in each direction
% OUTPUT
% X,Y,Z the meshgrids used for plotHydrogen()
% eVecs = the eigenvectors of H
% eVals = the eigenvalues of H

    % define grid
    L = 20;
    h = 2*L/(n+1);
    x = (-L+h:h:L-h)';
    y = x;
    z = x;
    
    % construct FD matrix
    I = speye(length(x));
    e = ones(length(x),1);
    D2 = spdiags([e -2*e e]/h^2, [-1 0 1], length(x), length(x));
    Delta = kron(kron(D2,I),I) + kron(kron(I,D2),I) + kron(kron(I,I),D2);
    V = externalPotential(x,y,z);
    H = -(1/2)*Delta - V;
    
    [eVecs,D] = eigs(H,12,'sa');
    eVals = diag(D);

    [X,Y,Z] = meshgrid(x,y,z);
end

function V = externalPotential(x,y,z)
    nx = length(x);
    ny = length(y);
    nz = length(z);
    X = kron(kron(x,ones(ny,1)),ones(nz,1));
    Y = kron(kron(ones(nx,1),y),ones(nz,1));
    Z = kron(kron(ones(nx,1),ones(ny,1)),z);
    absX = sqrt(X.^2+Y.^2+Z.^2);
    V = spdiags(1./(absX),[0],nx*ny*nz,nx*ny*nz);
end