%% Plot the results

figure;
for idx = 1:num_st
    subplot(num_st,1,idx);
    hold on;
    plot(t_lin,x_lin(:,x_lin_st(idx)));
    plot(t_nlin,x_nlin(:,x_nlin_st(idx)));
    hold off;
end