function [ poss, velocity ] = boundary_block( poss, velocity, cr,i,dt)
%This function evaluates the boundary conditions for the flow past a block
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

if(poss(i,1) > 4)
        x = poss(i,:);
        cp = [4 poss(i,2)];
        d = sqrt(dot(cp-x,cp-x));
        normal = [-1 0]; %Direction of corrected velocity
       
        ui = velocity(i,:);
        poss(i,:) =cp +  d*normal; %Update velocity to simulate collision with wall
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
    
    if abs(poss(i,1)-1.3) <= 0.05 && (abs(poss(i,2) -0.5) < 0.05) %flow around square boundary conditions
        x = poss(i,:);
        cp = [poss(i,1)+.01, poss(i,2)+.01];
        d = sqrt(dot(cp-x,cp-x));
        normal = [-1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp +  d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end


end

