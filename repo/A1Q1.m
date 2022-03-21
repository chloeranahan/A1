m0 = 9.1093837015E-31;
mn = 0.26*m0;
T = 300;
kB = 1.38064852E-23;
tmn = 0.2E-12;

vth = ((2*kB*T)/mn)^0.5; %thermal velocity
%ask if we use the 2 or 3 version of this equation

mfp = vth*tmn; %mean free path

%list of requirements: 
%x-y plane 200nm by 100nm
%1000 particles min (final)
%each partcle has velocity vth
%each particle has random direction
%use Newton's Laws of motion and fixed delta t to update positions
%plot and trace out trajectories for ~10 particles
%2D plot of particle subset that updates each timestep
%use  pause command for update loop
%specular reflection at the y boundaries
%periodic boundary condition at the x boundaries
%use vectors Px Py Vx Vy

V = vth;
% theta = rand(Np,1)*2*pi;

%uncommenting lines  26, 37 & 38, then commenting lines 40 & 41 gives the
%version in which theta is used to determine position (rather than Vx/Vy),
%but it is unreliable and only runs sometimes on Np  = 10

xmax = 200E-9; %max positions
ymax = 100E-9;

Np = 1000; % # particles, want 1000-10000
Nplot = 20;
Px = xmax*rand(Np,1); %initial positions
Py = ymax*rand(Np,1);

% Vy = V*sin(theta);
% Vx = V*cos(theta);

Vy = V*(rand(Np,1)-0.5); %initial velocities
Vx = sqrt((V^2)-(Vy.^2));
%Vx = V*(rand(Np,1)-0.5); %-0.5 added so electrons have a chance at having a positive or negative velocity

dt = 0.05*(ymax/V); %time step

tstop = 200; %simulation time

px = Px;
py = Py;

c = hsv(Nplot);
%plot(Xt,Yt,'-','color',TrajColors(j,:),'linewidth',3)

for i = 1:tstop
    px = Px;
    Px = Px + Vx*dt;
    
    py = Py;
    Py = Py + Vy*dt;
    
    ix1 = Px < 0;
    Px(ix1) = Px(ix1) + xmax;
    px(ix1) = px(ix1) + xmax;
    
    ix2 = Px > xmax;
    Px(ix2) = Px(ix2) - xmax;
    px(ix2) = px(ix2) - xmax;
    
    iy = Py < 0 | Py > ymax;
    Vy(iy) = -Vy(iy);
    
    
    figPlot = figure(1);
    xlabel('X (m)')
    ylabel('Y (m)')
    hold on
    for j = 1:Nplot
    plot([px(j),Px(j)]',[py(j),Py(j)]','color',c(j,:));
    end
    axis([0 xmax 0 ymax]);
    
    
    vavg = mean(sqrt(Vx.^2 + Vy.^2)); %average velocity
    
    TSi = ((vavg.^2)*mn)/(2*kB); %temperature of the Si semiconductor
    
    title(sprintf('T_S_i = %.3d',TSi))
        
    figure(2)
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    hold on
    plot(i,TSi,'go')

    pause(0.0001)
end
