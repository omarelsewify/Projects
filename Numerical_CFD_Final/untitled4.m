clear, clc, close all
syms rho u v P rhoL uL vL PL rhoR uR vR PR gam nux nuy 

W = [rho; rho*u; rho*v; P/(gam-1) + rho*(u^2 + v^2)/2];
WL = [rhoL; rhoL*uL; rhoL*vL; PL/(gam-1) + rhoL*(uL^2 + vL^2)/2];
WR = [rhoR; rhoR*uR; rhoR*vR; PR/(gam-1) + rhoR*(uR^2 + vR^2)/2];

V = [rho; u; v; P];

FxL = [rhoL*uL; rhoL*(uL^2) + PL; rhoL*uL*vL; (PL/(gam-1) + rhoL*(uL^2 + vL^2)/2 + PL)*uL];
FyL = [rhoL*vL; rhoL*(uL*vL); rhoL*vL^2 + PL; (PL/(gam-1) + rhoL*(uL^2 + vL^2)/2 + PL)*vL];

FxR = [rhoR*uR; rhoR*(uR^2) + PR; rhoR*uR*vR; (PR/(gam-1) + rhoR*(uR^2 + vR^2)/2 + PR)*uR];
FyR = [rhoR*vR; rhoR*(uR*vR); rhoR*vR^2 + PR; (PR/(gam-1) + rhoR*(uR^2 + vR^2)/2 + PR)*vR];

FnuL = [FxL,FyL]*[nux;nuy];
FnuR = [FxR,FyR]*[nux;nuy];

LHS = (FnuL - FnuR);
RHS = WL-WR;

% dWdV = jacobian(W,V)
% dFdV = jacobian(F,V)
% % 
% dFdW = simplify(inv(dWdV)*dFdV)

