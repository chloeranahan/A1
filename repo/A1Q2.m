m0 = 9.1093837015E-31;
mn = 0.26*m0;
T = 300;
kB = 1.38064852E-23;
tmn1 = 0.2E-12;

vth = ((2*kB*T)/mn)^0.5; %thermal velocity

mfp = vth*tmn1; %mean free path

V = vth;

xmax = 200E-9; %max positions
ymax = 100E-9;

Np = 1000; % # particles, want 1000-10000
Nplot = 14;
Px = xmax*rand(Np,1); %initial positions
Py = ymax*rand(Np,1);

Px1 = Px;
Py1 = Py;

Vy = V*(randn(Np,1)-0.5); %initial velocities
Vx = V*(randn(Np,1)-0.5); %-0.5 added so electrons have a chance at having a positive or negative velocity

dt = 0.01*(ymax/V); %time step
a = 1; %no acceleration

tstop = 200; %simulation time

Ppx = Px;
Ppy = Py;

c = hsv(Nplot);

box1 = 0; %turning on box 1
box2 = 0; %turning on box 2

if box1 == 1
b1 = rectangle('Position',[0.8E-7,0,0.4E-7,0.4E-7]); %box1 position (bottom)
end

if box2 == 1
b2 = rectangle('Position',[0.8E-7,0.6E-7,0.4E-7,0.4E-7]); %box2 position (top)
end

inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
while sum(inbox) > 0
    Px(inbox) = rand(sum(inbox),1)*xmax;
    Py(inbox) = rand(sum(inbox),1)*ymax;
    inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
end
    
   
figPlot = figure (2);
h = histogram(sqrt(Vx.^2 + Vy.^2),50);
title('Maxwell-Boltzmann Distribution')
xlabel('Velocity (m/s)')
ylabel('Number of Particles')

    
nmPaths = 0;
PathDistSum = 0;

for i = 1:tstop
    Ppx = Px;    
    Px = Px + Vx*dt; 
    
    Ppy = Py;
    Py = Py + Vy*dt;
    
    ix1 = Px < 0;
    Px(ix1) = Px(ix1) + xmax;
    Ppx(ix1) = Ppx(ix1) + xmax;
    
    ix2 = Px > xmax;
    Px(ix2) = Px(ix2) - xmax;
    Ppx(ix2) = Ppx(ix2) - xmax;
    
    iy = Py < 0 | Py > ymax;
    Vy(iy) = -Vy(iy);
    
    Psc = 1 - exp(-(dt/tmn1));
    std = sqrt((kB*T)/mn);
    
    isc = Psc > rand(Np,1);
    Vx = Vx + a*dt;
    Vy = Vy + a*dt;
    Px = Px + Vx*dt + 0.5*a*(dt)^2;
    Py = Py + Vy*dt + 0.5*a*(dt)^2;
    Vx(isc) = randn(sum(isc),1)*std;
    Vy(isc) = randn(sum(isc),1)*std;
    
    nmPaths = nmPaths + sum(isc);
    dist = sqrt((Px1(isc)-Px(isc)).^2 + ((Py1(isc)-Py(isc)).^2));
    PathDistSum = PathDistSum + sum(dist);
    
    Px1(isc) = Px(isc);
    Py1(isc) = Py(isc);
    
    if box1 == 1
        inbox1 = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7; %in box1
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox1 & LorR) = Ppx(inbox1 & LorR);
        Vx(inbox1 & LorR) = -Vx(inbox1 & LorR);
        Py(inbox1 & ~LorR) = Ppy(inbox1 & ~LorR);
        Vy(inbox1 & ~LorR) = -Vy(inbox1 & ~LorR);
    end
    
    if box2 == 1
        inbox2 = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7; %in box2
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox2 & LorR) = Ppx(inbox2 & LorR);
        Vx(inbox2 & LorR) = -Vx(inbox2 & LorR);
        Py(inbox2 & ~LorR) = Ppy(inbox2 & ~LorR);
        Vy(inbox2 & ~LorR) = -Vy(inbox2 & ~LorR);
    end
    
    
    figPlot = figure(1);
    xlabel('X (m)')
    ylabel('Y (m)')
    hold on
    for j = 1:Nplot
    plot([Ppx(j),Px(j)]',[Ppy(j),Py(j)]','color',c(j,:));
    end
    axis([0 xmax 0 ymax]);
    
    vavg = mean(sqrt(Vx.^2 + Vy.^2)); %average velocity
    
    TSi = ((vavg.^2)*mn)/(2*kB); %temperature of the Si semiconductor
    
    title(sprintf('T_S_i = %.3d',TSi))
    
        
    figure(3)
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    hold on
    plot(i,TSi,'go')

    pause(0.001)
end

avgmfp = PathDistSum/nmPaths;
tmn2 = avgmfp/vth;

