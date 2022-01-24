% Motion model

function [xk, Jac, Jac_wrt_w] = motion_model(x_prev, u_input, T_sampling)
assert(numel(x_prev)==3);
ang = x_prev(3);

cc = cos(ang);
ss = sin(ang);
xk_ = x_prev(:) + T_sampling * [cc, 0; ss, 0; 0, 1];
xk= xk_ * u_input;

if (nargout > 1 )
    Jac = eye(3);
    Jac(1,3) = - T_sampling * ss * u_input(1);
    Jac(2,3) = T_sampling * cc * u_input(1);

    Jac_wrt_w = T_sampling * [cc , 0; ss , 0; 0, 1];
end
