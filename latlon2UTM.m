%%%MATLAB����ʵ�־�γ��ת����ƽ������꣺
function [new_position,Ctrs,D,Idx] = latlon2UTM(data)
    M_PI=3.14159265358979323846;
    L = 6381372 * M_PI * 2; %�����ܳ�  
    W = L; % ƽ��չ����x������ܳ�  
    H = L / 2; % y��Լ�����ܳ�һ��  
    mill = 2.3; % ����ͶӰ�е�һ����������Χ��Լ������2.3֮�� 
    position = load(data);
%     geoplot(position(:,2),position(:,1),'r.')
%     hold on
%     geobasemap streets
%     opts = statset('Display', 'final');
%     [Idxplot, Ctrsplot] = kmeans(position, 1, 'Replicates', 10, 'Options', opts);
% 
%     geoplot(Ctrsplot(:,2),Ctrsplot(:,1),'kx','MarkerSize',14, 'LineWidth',2);
%     legend('Cluster 1','Centroids', 'Location','NE');
%     hold on;

    %%lon=120.7015202;%����
    %%lat=36.37423;%γ��
    n=size(position,1);
    new_position=zeros(n,2);

    for i =1:n
        lon=position(i,1);
        lat=position(i,2);
        x = lon * M_PI / 180; % �����ȴӶ���ת��Ϊ����  
        y = lat * M_PI / 180; %��γ�ȴӶ���ת��Ϊ����  
        y1 = -1.25 * log(tan(0.25 * M_PI + 0.4 * y)); % ����ͶӰ��ת��  
        % ����תΪʵ�ʾ���  
        dikaerX = (W / 2) + (W / (2 * M_PI)) * x ; %�ѿ�������x
        dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;%�ѿ�������y
        new_position(i,1)=dikaerX;
        new_position(i,2)=dikaerY;
    end
    
    opts = statset('Display', 'final');
    [Idx, Ctrs] = kmeans(new_position, 1, 'Replicates', 10, 'Options', opts);
    D=sqrt(abs((new_position(Idx==1,1)-Ctrs(1,1)).^2)+abs((new_position(Idx==1,2)-Ctrs(1,2)).^2));
%     figure
%     plot(new_position(Idx==1,1),new_position(Idx==1,2),'*r');
%     hold on
%     plot(Ctrs(:,1),Ctrs(:,2),'kx','MarkerSize',14, 'LineWidth',2);
%     legend('Cluster 1','Centroids', 'Location','NE');
%     hold on;
%     ylabel('UTM North (m)','Interpreter','Latex','FontSize', 12);
%     xlabel('UTM East (m)','Interpreter','Latex','FontSize', 12);
end