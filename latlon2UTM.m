%%%MATLAB程序实现经纬度转换成平面尔坐标：
function [new_position,Ctrs,D,Idx] = latlon2UTM(data)
    M_PI=3.14159265358979323846;
    L = 6381372 * M_PI * 2; %地球周长  
    W = L; % 平面展开后，x轴等于周长  
    H = L / 2; % y轴约等于周长一半  
    mill = 2.3; % 米勒投影中的一个常数，范围大约在正负2.3之间 
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

    %%lon=120.7015202;%经度
    %%lat=36.37423;%纬度
    n=size(position,1);
    new_position=zeros(n,2);

    for i =1:n
        lon=position(i,1);
        lat=position(i,2);
        x = lon * M_PI / 180; % 将经度从度数转换为弧度  
        y = lat * M_PI / 180; %将纬度从度数转换为弧度  
        y1 = -1.25 * log(tan(0.25 * M_PI + 0.4 * y)); % 米勒投影的转换  
        % 弧度转为实际距离  
        dikaerX = (W / 2) + (W / (2 * M_PI)) * x ; %笛卡尔坐标x
        dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;%笛卡尔坐标y
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