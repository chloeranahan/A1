% if Py(i+1) == 0
%     %0 <= theta(i) && theta(i) <= pi
%     %c = c + 1;
%     theta(i+1) = pi - theta(i);
% elseif Py(i+1) == ymax
%     %pi <= theta(i) && theta(i) <= 2*pi
%     c = c + 1;
%     theta(i) = 2*pi - (theta(i-1) - pi);
% else
%     c = c +1;
%     Py(i) = Py(i-1) + (Vy(i-1)*dt);
% end


% if Px(i) == 0
%     %c = c + 1;
%     Px(i+1) = xmax;
%     Py(i+1) = Py(i);
% elseif Px(i) == xmax
%     %c = c + 1;
%     Px(i+1) = 0;
%     Py(i+1) = Py(i);
% else
%     %c = c +1;
%     Px(i+1) = Px(i) + (Vx(i)*dt);
% end


% figPlot = figure(2);
% subplot(2,1,1);
% plot(hx);
% hold on
% subplot(2,1,2);
% plot(hy);
% hold on


%     figPlot = figure(1);
%     title(sprintf('T_Si_ = %.3d',TSi))
%     hold on
%     for j = 1:Nplot
%     plot([px(j),Px(j)]',[py(j),Py(j)]','color',c(j,:));
%     end
    
    
% function [ output_args ] = AddElectrons(Px, Py, Vx, Vy)
% 
% for p = 0:Np - 1
%     nAtoms = nAtoms + 1;
%     x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
%     y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
%     AtomType(nAtoms) = mn;
% end
% 
% for p = 1:Np
%     Vx(nAtoms - Np + p) = V * cos(PartAng);
%     Vy(nAtoms - Np + p) = V * sin(PartAng);
% end
% 
% end
% 
% AddElectrons();
% 
%  figure(2)
%     hold on
%     plot(i,TSi,'go')
% 
%     pause(0.001)

%MFP
%for loop for all the particles
%doesn't thave to be n, just index for the number of interactions
%all initial positions are called P1(n)
%if multiple Px are equal, set P2(n) = Px

%distance between scattering events
%MFP(n) = P2 - P1

%Px = xmax*rand(Np,1) - xmax;

% inbox2LR = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7;
%     if px < 0.8E-7 | px > 1.2E-7
%         Vx(inbox2) = -Vx(inbox2);
%     else
%         Vy(inbox2) = -Vy(inbox2);
%     end

% inbox1T = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7 & Ppy > 0.4E-7; %in box1, coming from top
%     Vy(inbox1T) = -Vy(inbox1T);
%     
%     inbox2LR = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7 & (Ppx < 0.8E-7 | Ppx > 1.2E-7); %in box2, coming from L or R
%     Vx(inbox2LR) = -Vx(inbox2LR);
%     
%     inbox2B = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7 & Ppy < 0.6E-7; %in box2, coming from bottom
%     Vy(inbox2B) = -Vy(inbox2B);


% inbox1 = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7; %in box1
%     LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
%     Px(inbox1 & LorR) = Ppx(inbox1 & LorR);
%     Vx(inbox1 & LorR) = -Vx(inbox1 & LorR);
%     Py(inbox1 & ~LorR) = Ppy(inbox1 & ~LorR);
%     Vy(inbox1 & ~LorR) = -Vy(inbox1 & ~LorR);
%     
%     inbox2 = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7; %in box2
%     LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
%     Px(inbox2 & LorR) = Ppx(inbox2 & LorR);
%     Vx(inbox2 & LorR) = -Vx(inbox2 & LorR);
%     Py(inbox2 & ~LorR) = Ppy(inbox2 & ~LorR);
%     Vy(inbox2 & ~LorR) = -Vy(inbox2 & ~LorR);



A3- FRI 5pm
% m0 = 9.1093837015E-31;
% mn = 0.26*m0;
% T = 300;
% kB = 1.38064852E-23;
% tmn = 0.2E-12;
% 
% vth = ((2*kB*T)/mn)^0.5; %thermal velocity
% 
% mfp = vth*tmn; %mean free path
% 
% V = vth;
% 
% xmax = 200E-9; %max positions
% ymax = 100E-9;
% 
% Np = 1000; % # particles, want 1000-10000
% Nplot = 14;
% Px = zeros(Np,1); %initial positions shifted so they come from the left side of the screen
% Py = ymax*rand(Np,1);
% 
% Vy = V*(randn(Np,1)-0.5); %initial velocities
% Vx = V*(abs(randn(Np,1))); %took away the -0.5 and changed randn back to rand so that all the inital x velocities are positive
% 
% dt = 0.01*(ymax/V); %time step
% a = 1; %no acceleration
% 
% tstop = 500; %simulation time
% 
% Ppx = Px; %previous postions
% Ppy = Py;
% 
% c = hsv(Nplot);
% 
% box1 = 0; %turning on box 1
% box2 = 0; %turning on box 2
% 
% if box1 == 1
% b1 = rectangle('Position',[0.8E-7,0,0.4E-7,0.4E-7]); %box1 position (bottom)
% end
% 
% if box2 == 1
% b2 = rectangle('Position',[0.8E-7,0.6E-7,0.4E-7,0.4E-7]); %box2 position (top)
% end
% 
% inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
% while sum(inbox) > 0
%     Px(inbox) = rand(sum(inbox),1)*xmax;
%     Py(inbox) = rand(sum(inbox),1)*ymax;
%     inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
% end
% 
% 
% for i = 1:tstop
%     Ppx = Px;
%     Ppy = Py;
%     
%     Px = Px + Vx*dt;
%     Py = Py + Vy*dt;
%     
%     
%     if box1 == 1
%         inbox1 = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7; %in box1
%         LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
%         Px(inbox1 & LorR) = Ppx(inbox1 & LorR);
%         Vx(inbox1 & LorR) = -Vx(inbox1 & LorR);
%         Py(inbox1 & ~LorR) = Ppy(inbox1 & ~LorR);
%         Vy(inbox1 & ~LorR) = -Vy(inbox1 & ~LorR);
%     end
%     
%     if box2 == 1
%         inbox2 = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7; %in box2
%         LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
%         Px(inbox2 & LorR) = Ppx(inbox2 & LorR);
%         Vx(inbox2 & LorR) = -Vx(inbox2 & LorR);
%         Py(inbox2 & ~LorR) = Ppy(inbox2 & ~LorR);
%         Vy(inbox2 & ~LorR) = -Vy(inbox2 & ~LorR);
%     end
%     
%     
%     %xBC = 1; %x boundary conditions ON
%     xBC = 0; %x boundary conditions OFF
%    
%     if xBC == 1
%         ix1 = Px < 0;
%         Px(ix1) = Px(ix1) + xmax;
%         Ppx(ix1) = Ppx(ix1) + xmax;
%        
%         ix2 = Px > xmax;
%         Px(ix2) = Px(ix2) - xmax;
%         Ppx(ix2) = Ppx(ix2) - xmax;
%     elseif xBC == 0
%         Px(ix1) = Px(ix1);
%         Ppx(ix1) = Ppx(ix1);
%         
%         Px(ix2) = Px(ix2);
%         Ppx(ix2) = Ppx(ix2);
%     end
%     
%     iy1 = Py < 0 | Py > ymax;
%     Vy(iy1) = -Vy(iy1);
% 
%     
%     %scatter = 1; %scattering is ON
%     scatter = 0; %scattering is OFF
%     
%     if scatter == 1
%         Psc = 1 - exp(-(dt/tmn));
%     elseif scatter == 0
%         Psc = 0;
%     end
%     
%     std = sqrt((kB*T)/mn);
%     
%     isc = Psc > rand(Np,1);
%     Vx = Vx + a*dt;
%     Vy = Vy + a*dt;
%     Px = Px + Vx*dt + 0.5*a*(dt)^2;
%     Py = Py + Vy*dt + 0.5*a*(dt)^2;
%     Vx(isc) = randn(sum(isc),1)*std;
%     Vy(isc) = randn(sum(isc),1)*std;
%     
%     
%     figPlot = figure(1);
%     hold on
%     axis([0 xmax 0 ymax]);
%     
%     for j = 1:Nplot
%     plot([Ppx(j),Px(j)]',[Ppy(j),Py(j)]','color',c(j,:));
%     end
%     
%     
%     vavg = mean(sqrt(Vx.^2 + Vy.^2)); %average velocity
%     
%     TSi = ((vavg.^2)*mn)/(2*kB); %temperature of the Si semiconductor
%     
%     title(sprintf('T_S_i = %.3d',TSi))
%     
%     figure(2)
%     hold on
%     plot(i,TSi,'go')
%     
%     pause(0.001)
% end

% figPlot = figure (2);
% PxPy = [Px Py];
% ElecDens = hist3(PxPy,'Nbins',[20 20]);
% title('Electron Density Map')
% xlabel('Px')
% ylabel('Py')

% hx = histogram(Vx);
% title('Maxwell-Boltzmann Distribution')
% xlabel('X Velocity (m/s)')
% ylabel('Number of Particles')
% figPlot = figure (3);
% hy = histogram(Vy);
% title('Maxwell-Boltzmann Distribution')
% xlabel('Y Velocity (m/s)')
% ylabel('Number of Particles')


% ElecDens = hist3(PxPy,'Nbins',[20 20]);

%,[20 10]