function[] = fn_plot_eGap(subData)
% plot energy gap values across all G
figure;

hold on, plot(subData.simParams.G_range,subData.E_max,'k');
hold on, plot(subData.simParams.G_range,subData.E_mean,'k:');

legend('E max', 'E mean', 'location', 'best');
% ylim([0.1 0.6]);
xlabel('G'); ylabel('Se difference between adjacent attractors');
title(['Energy gap profle across G']);

end