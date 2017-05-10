function [ poss, velocity ] = boundary_bucket( poss, velocity, cr, i, dt )
%This function evaluates the boundary conditions for the flow into a bucket
%simulation. If the conditions are met, it adjusts the position and
%velocity to simulate a colission.
%Inputs:
%poss: position of the particles
%velocity: velocity of the particles
%cr: Scaling factor for strength of interactions
%i:current index to look at
%dt:time step size
%Outputs:
%poss: updated position of the particles
%velocity: updated velocity of the particles

if poss(i,2) <= -1 && poss(i,1)>=1.5 && poss(i,1) <= 2.5 %bucket edge cases
        x = poss(i,:);
        cp = [poss(i,1)+0.01 poss(i,2)+0.01];
        d = sqrt(dot(cp-x,cp-x));
        normal = [0 1];
       
        ui = velocity(i,:);
        poss(i,:) =cp + d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
    if poss(i,2) >= -1 && poss(i,2) <= 0 && abs(poss(i,1)-1.5) <=0.05
        x = poss(i,:);
        cp = [poss(i,1)+0.01 poss(i,2)+0.01];
        d = sqrt(dot(cp-x,cp-x));
        normal = [1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp +  d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
    if poss(i,2) >= -1 && poss(i,2) <= 0 && abs(poss(i,1)-2.55) <=0.1
        x = poss(i,:);
        cp = [poss(i,1)+0.01 poss(i,2)+0.01];
        d = sqrt(dot(cp-x,cp-x));
        normal = [-1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp +  d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
    if(poss(i,1) < 0)
        x = poss(i,:);
        cp = [0 poss(i,2)];
        d = sqrt(dot(cp-x,cp-x));
        normal = [1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp + d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
end

