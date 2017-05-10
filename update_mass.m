function [pm_new ] = update_mass(n, poss, pm, md0, h)
%This function adjusts the initial mass of the particles. 
%Inputs:
%n: number of particles
%poss: position array of particles
%pm: initial particle mass
%h: kernel size/smoothing length
%md0: "rest density" corresponding to 0 pressure force when fluid is at
%rest
%Outputs:
%pm_new: the updated particle mass

%This could be optimized more by using a nearest neighbor algorithm, but
%we couldn't get it working... :(

d = 0; %density summation variable
for i = 1:n
    for j = 1:n
        r = poss(i,:) - poss(j,:); %get position vector
        if((r*r') < h^2) %condition for computing density
            d = d + Equations(r,h,1); %compute density using kernel
        end
    end
end
ad = d/n; %mean density
pm_new = (ad*md0)/(ad*ad); %update particle mass based on rest density and calculated mean density


end

