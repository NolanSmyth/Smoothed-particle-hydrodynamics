%Computational Project - Flow past block
% Initialize stuff
clear all
close all

b_size = 0.2; %Size of block to flow past
b_pos = [1.25, 0.45, .1, .1]; %Position of block's lower left corner
n = 10; %Size of particle grid is (n+1) x (n+1)
xbegin = 0.1; %Starting x value
ybegin = 0; %Starting y value
[X,Y] = meshgrid(xbegin:1/n:(xbegin+1),ybegin:1/n:(ybegin+1)); %Creat meshgrid for plotting
N = size(meshgrid(0:1/n:1,0:1/n:1)); %Size of initial grid
Nx = N(2);
Ny = N(1);
n = Nx*Ny; %New size of array
X = reshape(X,[n 1]); %Make X and Y into column vectors
Y = reshape(Y,[n 1]);

poss = zeros(n,2); %initialize position of every particle
poss(:,1) = X;
poss(:,2) = Y;

F_tot = zeros(n,2); %initialize acceleration matrix
density = zeros(n,1);
pressure = zeros(n,1);

velocity = [2 .* ones(n,1), zeros(n,1)]; %initialize velocity matrix with flow

% Constants/integration limits
pm = 1; %Particle mass (changes later)

x = 20;%Average number of particles in kernel
A = 1*1;%Area (m^2)
h = sqrt((A*x)/(n*pi));%Kernel size (m)

dt = 0.005; %integration time step size (may need to be adjusted)
k = 30; %Gas stiffness constant

g = -9.8;   %gravity constant m/s^2
md0 = 988;  %density (roughly that of water for now) kg/m^3
surfaceTension = 0.0728; %surface tension value recommended by Johansson N/m
surfaceLimit = sqrt(md0/x);
cr = 0.001; %Strength of interaction with boarders scaling factor
viscosityConstant = 3.5; %Recommended by Johansson (Pa*s)
runTime = 3;


% Video recording
mov = VideoWriter('Simulation_Block_test','Uncompressed AVI');
mov.FrameRate = round(1/dt);
open(mov);

%Write particles
figure
plot(poss(:,1),poss(:,2),'b.') %plot the initial position of the particles
hold on
rectangle('Position', b_pos, 'EdgeColor', 'r'); %plot block
hold off
axis tight
set(gca,'nextplot','replacechildren');
xlim([0 4])
ylim([-1.5 2.5])

% Adjust mass after desierd density
[pm] = update_mass(n,poss,pm,md0,h);

% Time loop
for t = 0:dt:runTime
    
    % Density and pressure
    
    [density, pressure] = update_pressure(n,poss,pm,md0,h,k);
    
    % Forces
    for i = 1:n
        fp = 0; %Pressure force
        fv = 0; %Viscosity force
       
        fg = 0; %no gravity for this case (top down view)
        
        for j = 1:n
            r = poss(i,:) - poss(j,:);
            if((r*r')< h^2)
                if r ~= 0
                    %Calculate pressure force
                    fp = fp +(((pressure(i)/(density(i)^2)))+((pressure(j)/(density(j)^2))))*pm*Equations(r,h,2);
                    %Calculate viscosity force
                    u = velocity(j,:) - velocity(i,:); %Find difference in velocities
                    fv = fv + u*(pm/density(j))*Equations(r,h,3);
                end
            end
        end
       
        
        fp = -density(i)*fp; %Correct pressure force based on ith particle's density
        fv = (viscosityConstant)*fv; %Correct viscocisty force based on viscosity constant
        
        fa = 0; %other force if needed
        F_tot(i,:) = (fp + fv + [0 fg] + fa); %Sum of forces
    end
    % Time integration
    for i = 1:n
        a = F_tot(i,:)/density(i); %Calculate acceleration (F/m)
        % Use Modified Verlet method for simple integration
        velocity(i,:) = velocity(i,:) + a*dt;
        poss(i,:) = poss(i,:) + velocity(i,:)*dt + 0.5*a*dt^2;
        %end
        
        % Boundary conditions (different for each scenerio)
        [poss, velocity] = boundary_block(poss, velocity, cr,i,dt);
        
        
    end
    
    % Save movie
    currFrame = getframe;
    writeVideo(mov,currFrame);
    
    plot(poss(:,1),poss(:,2),'b.')
    xlim([0 4])
    ylim([-1.5 2.5])
    hold on
    rectangle('Position', b_pos, 'EdgeColor', 'r'); %plot block
    hold off
    xlim([0 4])
    ylim([-1.5 2.5])
end
close(mov)