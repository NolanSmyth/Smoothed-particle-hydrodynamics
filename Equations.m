function [ Equation_ret ] = Equations(r,h,n)
%This function contains the equations needed to solve the Navier Stokes
%equations for a smooth particle hydrodynamics (SPH) simulation. 

%Inputs:
%r: 2 dimensional position vector
%h: characteristic width of the kernel
%n: Input determining which equation to return

Equation_ret = 0;

if(n == 1)
    %Kernel smoother used to approximate dirac delta function when
    %interpolating properties of "fluid" particles.
    Equation_ret = (315/(64*pi*h^9))*(h^2 -(r*r'))^3;
end

if(n == 2)
    %Gradient of the Kernel. Used to find the presure when solving navier
    %stokes.
    Equation_ret = -(45/(pi*h^6))*(r/(sqrt(r*r')))*(h-sqrt((r*r')))^2;
end

if(n == 3)
    %Laplacian of the Kernel. Used to find the viscocity when solving navier stokes. 
    Equation_ret = (45/(pi*h^6))*(h-(sqrt((r*r'))));  
end

end