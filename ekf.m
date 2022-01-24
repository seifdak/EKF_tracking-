load dataset2.mat

% Call motion model
T = t(2)-t(1); % Sampling periods of the 10 Hz laser rangefinder scans
idx = 12609;
u_input = [v(idx); om(idx)];
x = [x_true(idx-1); y_true(idx-1); th_true(idx-1)];
% Call motion model with dataset pararms
xk = motion_model(x, u_input, T);
% Compute numerical jacobian
F_Jac = fdjac(x, 'motion_model', u_input, T);
% Implement Analytcal Jacobians
[xk, Jac, Jac_wrt_w] = motion_model(x, u_input, T);
%Compare both Jacobian calc. methods
motion_fro_norm = norm(Jac - F_Jac, 'fro') / norm(Jac,'fro');

%Call observation model
idx_1 = 17;
y_measured = [r(idx,idx_1); b(idx,idx_1)];
xk = [x_true(idx); y_true(idx); th_true(idx)];

y_pred = observation_model(xk, d, l(idx_1,:)'); % Innovation

G_Jac = fdjac(xk, 'observation_model', d, l(idx_1,:)');

[y_pred, G_Jac_ex] = observation_model( ...
    xk, ...
    d, ...
    l(idx_1,:)');
obs_fro_norm= norm(G_Jac_ex - G_Jac, "fro") / norm(G_Jac_ex,'fro');

% Initialize EKF
rmax = 1;
% x0 = [x_true(1); y_true(1); th_true(1)];
% 4.b start with a poor initial condition
 x0 = [1, 1, 0.1]';
P0 = diag([1 1 0.1]);

%covariances
num_landmarks = 17;
covar_speeds = diag([v_var, om_var]);
var_num_valid_meas = repmat([r_var, b_var], 1 , num_landmarks);
covar_obs_noise_all = diag(var_num_valid_meas);  % R cov matrix

% State and cov EKF
state_cur = x0;
cov_cur = P0;
num_steps = numel(t);
state_history = zeros(3, num_steps);
state_cov_history = zeros(3, 3, num_steps);

% RUN EKF
tic
b_eval_Jac_at_gt = false;
for k= 1:num_steps
    %Prediction mean
    u_input =[v(k); om(k)];
    if k == 1
        x = state_cur;
    else
        x = [x_true(k-1); y_true(k-1); th_true(k-1)];
    end
    [state_pred, Jac2, Jac_wrt_w] = motion_model(x, u_input, T);
    %state_history(:,k) = state_pred;
    %Prediction covariance
    cov_mot_noise_all = Jac_wrt_w * covar_speeds * Jac_wrt_w.';
    cov_pred = Jac2 * state_cov_history(:,:,k) * Jac2.' + cov_mot_noise_all;
    %state_cov_history(:,:,k) = cov_pred;
    %correction step
    % Call observation_model(state_pred, d, l) with paramatrized max range
    for i = 1:num_landmarks
        if  r(k,i) + d < rmax && r(k,i) ~= 0 && b(k,i) ~= 0
            % Kalman Gain
            if b_eval_Jac_at_gt == false
                [state_corr, G_Jac_ana] = observation_model([x_true(k), y_true(k), th_true(k)]', d, l(i,:)');
        else
            [state_corr, G_Jac_ana] = observation_model(state_history(:,k), d, l(i,:)');
            end
            K_gain = cov_pred * G_Jac_ana.' / inv(G_Jac_ana * cov_pred * G_Jac_ana.' + covar_obs_noise_all(1:2,1:2));
            y_measured = [r(k, i); b(k,i)];
            %compute posterior
            state_post =  state_pred + K_gain * (y_measured - state_corr); %x_hat = x_check + K_gain * innovation
            % correction covariance
            cov_corr = (1 - K_gain * G_Jac_ana) * cov_pred;
            state_cov_history(:,:,k) = cov_corr;

        else
            fprintf('landmark %d was skipped at timestemp %d.\n',i,k)
        end
            % Save states in state_history matrix (3x12609)
            state_history(:,k) = state_post;
    end 
end
toc

% Extract and calculate the uncertainty envelopes for x, y and th
upper_bound_x = 3 * sqrt(diag(state_cov_history(1:1,:,num_steps)));
lower_bounx_x = - 3 * sqrt(diag(state_cov_history(1:1,:,num_steps)));
upper_bound_y = 3 * sqrt(diag(state_cov_history(2:2,:,num_steps)));
lower_bound_y = - 3 * sqrt(diag(state_cov_history(2:2,:,num_steps)));
upper_bound_th = 3 * sqrt(diag(state_cov_history(3:3,:,num_steps)));
lower_bound_th = - 3 * sqrt(diag(state_cov_history(3:3,:,num_steps)));

% Plot figures for estimation erros with user defined parameter (rmax)
figure, plot(t, state_history(1:1,1:num_steps) - x_true.')
xlabel('time [s]')
ylabel('x-position error')

figure, plot(t, state_history(2:2,1:num_steps) - y_true.')
xlabel('time [s]')
ylabel('y-position error ')

figure, plot(t, state_history(3:3,1:num_steps) - th_true.')
xlabel('time [s]')
ylabel('heading-position error ')