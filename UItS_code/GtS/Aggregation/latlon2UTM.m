function [new_position,Idx,Ctrs,D_center,Ctrs_center] = latlon2UTM(data,GW_num)
    rng(1);
    M_PI=3.14159265358979323846;
    L = 6381372 * M_PI * 2; 
    W = L;  
    H = L / 2;
    mill = 2.3; 
    position = load(data);
    n=size(position,1);
    new_position=zeros(n,2);

    for i =1:n
        lon=position(i,1);
        lat=position(i,2);
        x = lon * M_PI / 180;  
        y = lat * M_PI / 180; 
        y1 = -1.25 * log(tan(0.25 * M_PI + 0.4 * y));
        dikaerX = (W / 2) + (W / (2 * M_PI)) * x ;
        dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;
        new_position(i,1)=dikaerX;
        new_position(i,2)=dikaerY;
    end
    GW_number=GW_num;
    opts = statset('Display', 'final');
    [Idx, Ctrs] = kmeans(new_position, GW_number, 'Replicates', 10, 'Options', opts);
    [Idx_center, Ctrs_center] = kmeans(new_position, 1, 'Replicates', 10, 'Options', opts);

    D_center=sqrt(abs((new_position(Idx_center==1,1)-Ctrs_center(1,1)).^2)+abs((new_position(Idx_center==1,2)-Ctrs_center(1,2)).^2));
end