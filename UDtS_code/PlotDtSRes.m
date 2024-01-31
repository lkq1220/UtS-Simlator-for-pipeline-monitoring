%plot data aggregation
close all;
clear;
clc;

figure
load("1e5DtSRes.mat");
% figure
set(gcf, 'Units', 'centimeters'); 
afFigurePosition = [2 7 14 11];
set(gcf, 'Position', afFigurePosition,'PaperSize',[18 8],'PaperPositionMode','auto'); 
plot(Distance/1000,PSNR_UDtS,'-','linewidth',2.5,'SeriesIndex',1);
hold on 
plot(Distance/1000,PSIR_UDtS,'-','linewidth',2.5,'SeriesIndex',2);
hold on 
plot(Distance/1000,PSF11_UDtS,'-','linewidth',2.5,'SeriesIndex',3);
hold on 
grid on
ylabel('Success probability', 'Interpreter', 'Latex','fontsize',14);
xlabel('Distance from city center to satellite (km)','Interpreter','Latex','FontSize', 14);
axis([Distance(1)/1000 Distance(end)/1000 0 1]);
legend('$P_{SNR}^{DtS}$','$P_{SIR}^{DtS}$','$P_{S}^{DtS}$','Interpreter', 'Latex','fontsize',15,'Location','northeast');
set(gca,'fontsize',15);
