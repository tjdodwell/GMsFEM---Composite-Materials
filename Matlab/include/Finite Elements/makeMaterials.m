function [isotropic,composite] = makeMaterials(E_R,nu_R,E1,E2,E3,nu_21,nu_31,nu_32,G_12,G_13,G_23)

% Isotropic Resin
lambda = E_R*nu_R/((1+nu_R)*(1-2*nu_R)); mu = E_R/(2*(1+nu_R));
isotropic = zeros(6);
isotropic(1:3,1:3) = lambda*ones(3);
isotropic = isotropic + diag(mu*[2,2,2,1,1,1]);

% Orthotropic Composite
S = zeros(6);
S(1,1) = 1/E1; S(1,2) = -nu_21/E2; S(1,3) = -nu_31/E3;
S(2,1) = S(1,2);    S(2,2) = 1/E2; S(2,3) = -nu_32/E3;
S(3,1) = S(1,3);    S(3,2) = S(2,3);    S(3,3) = 1/E3;
S(4,4) = 1/G_23;
S(5,5) = 1/G_13;
S(6,6) = 1/G_12;

composite = inv(S);


end