function [Satellite_Link_Farms,Satellite_Link_Center_found] = Simulations_Distance_Points(Earth_Radius,Orbit_Height,Grid,Center_found)

Elevation_Simulation = [90:-1:10];
Slant_Range=zeros(1,length(Elevation_Simulation));
Base =zeros(1,length(Elevation_Simulation));
    
Location_X_Y=Grid;
    
for WhichAngle=1:length(Elevation_Simulation)  
    Eq_pt_1 = cosd(Elevation_Simulation(WhichAngle)).*cosd(Elevation_Simulation(WhichAngle));  
    Eq_pt_2 = ((Orbit_Height + Earth_Radius)./Earth_Radius)^2;
    
    Slant_Range(WhichAngle) = Earth_Radius.*(sqrt(Eq_pt_2- Eq_pt_1) - sind(Elevation_Simulation(WhichAngle)));  % Slant_Range(E): distance from user to satellite
        
        if(Elevation_Simulation(WhichAngle)==90)
            Base (WhichAngle) = abs(sqrt(Slant_Range(WhichAngle)^2-Orbit_Height^2));
        else
            Base (WhichAngle) = sqrt(Slant_Range(WhichAngle)^2-Orbit_Height^2);
        end
    
%% simulation points (distance): Satellite mobility
    
    Satellite_Link_Farms(WhichAngle,:,:)=[(Location_X_Y(:,1)+ Base(WhichAngle)), Location_X_Y(: ,2)];
    Satellite_Link_Center_found(WhichAngle,:,:)=[(Center_found(1,1)+ Base(WhichAngle)), Center_found(1 ,2)];
end

end