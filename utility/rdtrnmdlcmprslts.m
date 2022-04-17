%% Compute the results

for idx = 1:num_st
    % Compute the absolute error between linear and nonlinear models
    e_abs(:,idx) = ...
        abs(x_lin(:,x_lin_st(idx)) - x_nlin(:,x_nlin_st(idx)));
    % Compute the relative estimation error
    % for the point of maximum absolute actual value of the real signal
    [x_nlin_max(idx),max_idx(idx)] = max(x_nlin(:,x_nlin_st(idx)));
    e_rel(idx) = max(e_abs(max_idx(idx),idx)./x_nlin_max(idx));
    
    % Compute the normalized cross-correlation
    [r(:,idx),lag(:,idx)] = ...
        xcorr(x_lin(:,x_lin_st(idx)),x_nlin(:,x_nlin_st(idx)),...
        100,...
        'normalized');
    [r_max(idx),r_max_idx(idx)] = max(r(:,idx));
    lag_max(idx)                = lag(r_max_idx(idx))*t_s;
    r_0(idx)                    = r(lag(:,idx) == 0,idx);
end