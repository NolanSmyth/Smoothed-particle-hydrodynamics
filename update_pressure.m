function [ density, pressure ] = update_pressure(n, poss, pm, md0,h,k )
%This function calculates the pressure for each particle through an SPH
%method
%Inputs:
%n: number of particles
%poss: position array of particles
%pm: initial particle mass
%md0: "rest density" corresponding to 0 pressure force when fluid is at
%rest
%h: kernel size/smoothing length
%k: Gas stiffness constant
%Outputs:
%ad: mean density over all interacting particles
%pm_new: the updated particle mass

for i = 1:n
    d = 0;
    for j = 1:n
        r = poss(i,:) - poss(j,:);
        if((r*r') < h^2)
            d = d + pm*Equations(r,h,1); %calculate the density for each particle
        end
    end
    density(i) = d; %save this value in an array for later use
    pressure(i) = k*(d-md0); %calculate the pressure force based on density and rest density
end

end

