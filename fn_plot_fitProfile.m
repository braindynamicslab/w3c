function [] = fn_plot_fitProfile(subData)

    G_select1 = find(subData.fitFC_xAttr_r == max(subData.fitFC_xAttr_r));

    figure, subplot(3,2,1), imagesc(subData.FC);
    title('Real FC');
    
    hold on, subplot(3,2,2), imagesc(subData.SC);
    title('SC');
    
    hold on, subplot(3,2,3), imagesc(subData.xAttrFC{G_select1(1)});
    title(['xAttr Conn, G = ', num2str(subData.simParams.G_range(G_select1(1)))]);
    
    hold on, subplot(3,2,4),plot(subData.simParams.G_range,subData.fitFC_xAttr_r,'k');
    hold on, subplot(3,2,4),plot(subData.simParams.G_range,subData.fitFC_xAttrIntra_r,'r');
    hold on, subplot(3,2,4),plot(subData.simParams.G_range,subData.fitFC_xAttrInter_r,'b');
    legend('Full', 'Intra', 'Inter', 'location', 'best');
    ylim([0.1 0.6]);
    xlabel('G'); ylabel('R');
    title(['xAttr Fit']);
    
    hold on, subplot(3,2,6),plot(subData.simParams.G_range,subData.fitFC_xAttr_r_condSC,'k');
    hold on, subplot(3,2,6),plot(subData.simParams.G_range,subData.fitFC_xAttrIntra_r_condSC,'r');
    hold on, subplot(3,2,6),plot(subData.simParams.G_range,subData.fitFC_xAttrInter_r_condSC,'b');
    legend('Full', 'Intra', 'Inter', 'location', 'best');
    ylim([0.1 0.6]);
    xlabel('G'); ylabel('R');
    title(['xAttr Fit (condSC)']);
        
    set(gcf, 'position', [500 90 900 700]);
    
    suptitle([subData.subID]);
end