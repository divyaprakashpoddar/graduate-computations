% Simulates the motion of an elliptical cylinder under the influence of a body force in a channel using Lattice Boltzmann Method
%
% Author: Divyaprakash (2020AMZ8461), PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 14 January 2021

clear; clc; close all;
iId = 0;

% Simulation parameter values
a       = 8;                            % Semi-Major Axis [Lattice Units]
b       = 26;                           % Semi-Minor Axis [Lattice Units]
nx      = 1000;                         % No. of Lattice Units (x-Direction)
ny      = 1000;                          % No. of Lattice Units (y-Direction)
T       = 30000;                        % Simulation Time


% Get grid                                            
ci      = 1*nx/10;                      % Cylinder center (x)
cj      = 0.2*ny;                       % Cylinder center (y)
theta   = 45; 
% isn: [0--> Fluid, 1--> Solid], icn: [0--> fluid, 1--> Fictitious Fluid]
[isn, icn, xb, yb, nb, X] = defineEllipse(a,b,nx,ny,ci,cj,theta);
%writematrix(isn,'isn_0.txt');
% Visualize the Mesh and the boundaries
% hold on
% spy(isn','b.')
% spy((1-isn)','r.')
% plot(xb,yb,'k.','MarkerSize',10)

% LBM Parameters
na      = 9;                                    % D2Q9
rho0    = 1;                                    % Initial Density
tau     = 0.55;                                  % Relaxation Time
cs      = 1/sqrt(3);                            % Speed of Sound
nu      = cs^2*(tau - 1/2);                     % Kinematic Viscosity
ex      = [0 1 0 -1 0 1 -1 -1 1];
ey      = [0 0 1 0 -1 1 1 -1 -1];
wt      = [4/9 1/9*ones(1,4) 1/36*ones(1,4)];   % Weights

% Wall Boundary
uWall   = 0;
uWT     = -uWall;           % wall Velocity
uWB     = uWall;            % wall Velocity
    

% Define and preallocate the variables
f           = zeros(na,nx,ny);          % Distribution function
ft          = zeros(na,nx,ny);          % Post Collision distribution
feq         = zeros(1,na);              % Equilibrium distribution at a node
UX          = 0.0*ones(nx,ny);         % x-velocity
UY          = 0.0*ones(nx,ny);         % y-velocity
    
% Initialize the distribution function
ucx = 0; ucy = 0; % Velocity of fluid inside the cylinder
for i = 1:nx
    for j = 1:ny
        for a = 1:na
            u2 = ucx*ucx + ucy*ucy;
            term1 = ucx*ex(a) + ucy*ey(a);
            term2 = term1*term1;
            f(a,i,j) = wt(a)*rho0*(1 + 3*term1 + 4.5*term2 - 1.5*u2);
        end
    end
end

% Force calculation
Fbx = zeros(1,nb);  % X-direction force at each boundary node
Fby = zeros(1,nb);  % Y-direction force at each boundary node
FRX = zeros(1,T);   % X-direction resultant force
FRY = zeros(1,T);   % Y-direction resultant force
FXTemp = 0;
FYTemp = 0;

% Torque calculation
TqR = zeros(1,T); %Total Torque
TqTemp = 0;

% Cylinder motion
UCX = zeros(1,T);        % Translation velocity of cylinder
UCY = zeros(1,T);
omega = zeros(1,T);    % Rotational velocity of cylinder
UBX = zeros(1,nb);   % Vector of all the boundary node velocities
UBY = zeros(1,nb);   % Vector of all the boundary node velocities
CX = zeros(1,T);
CY = zeros(1,T);
CD = zeros(1,T);
CL = zeros(1,T);
Theta = zeros(1,T);

% Mass and (Polar) Moment of Inertia
rhoC    = 2*rho0;           % Density of the cylinder
M       = rhoC*(pi*a*b);
I       = pi/4*(a*b^3 + a^3*b);

% Applied Force
g = 1e-4;
FA = M*g;

% Start the time loop
for t = 1:T
    
    if mod(t,100) == 0
        % Print data every 100 time steps
        clc
        fprintf('Time Step: %d (max: %d)\n',t,T)
        fprintf('Maximum Velocity: %.4f\n',max(max(max(UX)),max(max(UY))))
        fprintf('Average Density: %.4f\n',rho_avg)
        fprintf('Drag Force: %.8f\n',FRX(t-1))
        fprintf('Cylinder Velocity (x): %.8f\n',UCX(t-1))
        fprintf('Angular Velocity: %.8f\n',omega(t-1))
        fprintf('cx: %.4f\n',ci)
        fprintf('cy: %.4f\n',cj)
    end
    
    % Initialize the forces
    Fbx = zeros(1,nb); % X-direction force at each boundary node
    Fby = zeros(1,nb); % Y-direction force at each boundary node
    
    % Average density
    rho_avg = 0;
    % Calculate rho and u
    for i = 1:nx
        for j = 1:ny
            % Initialize and calculate ux, uy and rho for a node
            ux = 0; uy = 0; rho = 0; 
            for a = 1:na
                rho = rho + f(a,i,j);
                ux = ux + f(a,i,j)*ex(a); % These velocities are for
                uy = uy + f(a,i,j)*ey(a); % a node.
            end
            ux = ux/rho;
            uy = uy/rho;
            u2 = ux*ux + uy*uy;

            % Store the velocities
            UX(i,j) = ux;
            UY(i,j) = uy;
            rho_avg = rho + rho_avg;
            % Calculate the equilibrium distribution
            for a = 1:na
                term1 = ux*ex(a) + uy*ey(a);
                term2 = term1*term1;
                feq(a) = wt(a)*rho*(1 + 3*term1 + 4.5*term2 - 1.5*u2);

                % Calculate the post collision distribution
                ft(a,i,j) = f(a,i,j) - (f(a,i,j) - feq(a))/tau;
            end
        end
    end

    % Average density
    rho_avg = rho_avg/nx/ny;
    
    % Streaming of post collision distribution (Fluid nodes)
    for i = 1:nx
        for j = 1:ny
            if (isn(i,j) == 0)
                for a = 1:na
                    ia = i + ex(a);
                    ja = j + ey(a);

                    % Periodic Boundary Condition on the left and right
                    % boundary
                    if ia < 1;  ia = nx; end
                    if ia > nx; ia = 1; end
                        
                    f(a,ia,ja) = ft(a,i,j);
                   
                end
            end
        end
    end
    
     % Streaming of post collision distribution (Fictitious Fluid nodes)
    rev = [1 4 5 2 3 8 9 6 7]; % Reversed directions of distributions
     for i = 1:nx
        for j = 1:ny
            if (icn(i,j) == 1)
                for a = 1:na
                    ia = i + ex(a);
                    ja = j + ey(a);

                    if (icn(ia,ja) == 0) % It means a boundary is encountered
                        
                        % Get the coordinate of the boundary node
                        XB = i + 1/2*ex(a);
                        YB = j + 1/2*ey(a);
                        bID = and(xb == XB, yb == YB);
                        
                        % Get the boundary node velocity
                        ubx = UBX(bID);
                        uby = UBY(bID);
                        
                        % Unknown: f(rev(a)), Unknown direction: rev(a)
                        f(rev(a),i,j) = ft(a,i,j) + (2*wt(rev(a))*rho0/cs^2)*dot([ex(rev(a)), ey(rev(a))], [ubx uby]);
                    else
                        f(a,ia,ja) = ft(a,i,j);
                    end
                end
            end
        end
    end

    % Apply Boundary Condition: Cylinder Walls
    rev = [1 4 5 2 3 8 9 6 7]; % Reversed directions of distributions
    for i = 3:nx-3
        for j = 3:ny-3
            if (isn(i,j) == 0)
                for a = 1:na
                    ia = i + ex(a);
                    ja = j + ey(a);
                    if (isn(ia,ja) == 1) % It means a boundary is encountered
                        % Get the coordinate of the boundary node
                        XB = i + 1/2*ex(a);
                        YB = j + 1/2*ey(a);
                        bID = and(xb == XB, yb == YB);
                        
                        % Get the boundary node velocity
                        ubx = UBX(bID);
                        uby = UBY(bID);
                        
                        % Calculate force at the boundary nodes
                        Fbx(bID) = Fbx(bID) + 2*(ft(a,i,j) - ft(rev(a),ia,ja)...
                                    - 2*wt(a)*rho0/cs^2*dot([ex(a), ey(a)], [ubx uby]))*ex(a);
                        Fby(bID) = Fby(bID) + 2*(ft(a,i,j) - ft(rev(a),ia,ja)...
                            - 2*wt(a)*rho0/cs^2*dot([ex(a), ey(a)], [ubx uby]))*ey(a);

                        
                        f(rev(a),i,j) = ft(a,i,j) + (2*wt(rev(a))*rho0/cs^2)*dot([ex(rev(a)), ey(rev(a))], [ubx uby]);
                    end
                end
            end
        end
    end
     
    % Moving Boundary Condition on the top wall
    %rev = [1 4 5 2 3 8 9 6 7]; % Reversed directions of distributions

    j = ny-1;
    for i = 1:nx
            for a = [8 5 9]
                    f(a,i,j) = ft(rev(a),i,j) + (2*wt(a)*rho0/cs^2)*dot([ex(a) ey(a)],[uWT 0]);
            end
    end
    
    % Moving Boundary Condition on the bottom wall

    j = 2;
    for i = 1:nx
            for a = [7 3 6]
                    f(a,i,j) = ft(rev(a),i,j) + (2*wt(a)*rho0/cs^2)*dot([ex(a) ey(a)],[uWB 0]);
            end
    end
        
    
    if mod(t,20) == 0
        % Plot velocity contour every 100 time steps
        UMag = sqrt(UX.^2 + UY.^2);
        iId = iId+1;
        contourf(flipud(UMag(:,2:end-1)))
        % colorbar
        % title(sprintf('t = %d',iId))
        axis equal
        axis off
        fName = sprintf('image-%d.png',iId);
        saveas(gca,fName)
    end
    
    if mod(t,100) == 0
        % Write data files every 1000 time steps
        fName = sprintf('U_%d.txt',t);
        writematrix(UX,fName)
        fName = sprintf('V_%d.txt',t);
        writematrix(UY,fName)
        
        %fName = sprintf('isn_%d.txt',t);
        %writematrix(isn,fName)
        
        %writematrix(UCX,'UCX.txt')
        %writematrix(UCY,'UCY.txt')
        %writematrix(FRY,'FRY.txt')
        %writematrix(FRX,'FRX.txt')
        %writematrix(TqR,'TqR.txt')
        %writematrix(omega,'omega.txt')
        %writematrix(CX,'CX.txt')
        %writematrix(CY,'CY.txt')
        %writematrix(CD,'CD.txt')
        %writematrix(CL,'CL.txt')
        
    end
    
    if mod(t,500) == 0
        %spy(flipud(isn'))
        %title(sprintf('t = %d',t))
        %axis equal
        %fName = sprintf('isn_f%d.png',t);
        %saveas(gca,fName)
    end
    
    % Torque at boundary nodes
    Tq = 0; % Initialize
    for n = 1:nb
        xbVec = [xb(n), yb(n), 0];
        tempV = cross((xbVec - [X 0]), [Fbx(n) Fby(n) 0]);
        Tq = Tq + tempV(3);
    end
    
    
    % Calculate Total Force and Torque
    if t ~= 1
        FRX(t) = 1/2*(sum(Fbx) + FXTemp) + FA; % Add Force;
        FRY(t) = 1/2*(sum(Fby) + FYTemp);
        TqR(t) = 1/2*(Tq + TqTemp);
    end
    FXTemp = sum(Fbx);
    FYTemp = sum(Fby);
    TqTemp = Tq; % Add Forced Torque
    
    % Update translational and rotational velocity of the cylinder
    if mod(t,2) == 0
        UCX(t+1) = UCX(t-1) + 2/M*FRX(t);
        UCX(t+2) = UCX(t+1);
        
        UCY(t+1) = UCY(t-1) + 2/M*FRY(t);
        UCY(t+2) = UCY(t+1);
        
        omega(t+1) = omega(t-1) + 2/I*TqR(t);
        omega(t+2) = omega(t+1);
    end
   
    % Update the position of the cylinder
    if mod(t,2) == 0
        ci = ci + 2*UCX(t);
        cj = cj + 2*UCY(t);
        theta = theta + 2*omega(t);
        [isn, icn, xb, yb, nb, X] = defineEllipse(a,b,nx,ny,ci,cj,theta);
    end
    
    
    % Update the boundary velocity
    for n = 1:nb
        xbVec = [xb(n), yb(n), 0];
        UBTemp = [UCX(t) UCY(t) 0] + cross([0 0 omega(t)], (xbVec - [X 0]));
        UBX(n) = UBTemp(1);
        UBY(n) = UBTemp(2);      
    end
   CX(t) = ci; CY(t) = cj;
end

% % Visualize the forces on the cylinder at each node
% FbR = sqrt(Fbx.^2+Fby.^2);
% Fbx = Fbx./FbR;
% Fby = Fby./FbR;
% 
% quiver(xb,yb,Fbx,Fby)
% axis equal
% hold on
% scatter(xb,yb)
