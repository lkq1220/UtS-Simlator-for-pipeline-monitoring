function [new_position,Ctrs,D,Idx] = latlon2UTM(data)
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

        dikaerX = (W / 2) + (W / (2 * M_PI)) * x;
        dikaerY = (H / 2) - (H / (2 * mill)) * y1;
        new_position(i,1)=dikaerX;
        new_position(i,2)=dikaerY;
    end
    
    opts = statset('Display', 'final');
    [Idx, Ctrs] = kmeans(new_position, 1, 'Replicates', 10, 'Options', opts);
    D=sqrt(abs((new_position(Idx==1,1)-Ctrs(1,1)).^2)+abs((new_position(Idx==1,2)-Ctrs(1,2)).^2));
end