

function [yk_pred, Jac] = observation_model(xk, d, xl)
ang = xk(3);
cc = cos(ang);
ss = sin(ang);
% xl is the stack of X and Y (intermediate variables)
dx = xl - xk(1:2) - d*[cc; ss];
dx_sqrd_norm = dx' * dx;
dx_norm = sqrt(dx_sqrd_norm);
yk_pred = [dx_norm; atan2(dx(2), dx(1)) - ang];
yk_pred(2) = wrapToPi(yk_pred(2));

if (nargout > 1)
    Jac = zeros(2,3);

    Jac(1,1) = -dx(1) / dx_norm;
    Jac(1,2) = -dx(2) / dx_norm;
    Jac(1,3) = d*(dx(1)*ss-dx(2)*cc) / (dx_norm);

    Jac(2,1) = dx(2) / dx_norm * dx_norm ;
    Jac(2,2) = -dx(1) / dx_norm * dx_norm;
    Jac(2,3) = -d * ((dx(1) * cc + dx(2) * ss) / dx_norm * dx_norm) -1;
end

