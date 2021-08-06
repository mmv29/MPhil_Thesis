%% Model inputs
% epiboly_progress range 0-1
% start from 20% epiboly (2hps) to 100% epiboly (6hps) for control embryo
% beta = 0.20 i.e. 80% in 4 hours

% define inputs (initial conditions)
[embryo]=[  ...
    0.10,   ... %alpha          (division rate R_y)         [1/time]          [1/hours]
    30,     ... %gamma          (decay rate)                [1/length]        [1/mm]
    3000,   ... % N0            (initial cell number)       [none]            [cells]
    2,      ... % t0            (start time)                [time]            [hours]
    0.1,    ... % dt            (time step)                 [time]            [hours]
    6,      ... % t1            (end time)                  [time]            [hours]
    0.33,   ... % x0            (yolk radius R_y)           [length]          [mm]
    0.003,  ... % dx            (thickness step)            [length]          [mm]
    5.7e+4, ... % density       (cell density)              [1/volumen]       [cells/mm^3]
    0.2,    ... % beta          (epiboly rate)              [1/time]          [%epi/hour]
    0.2     ... % ini epi       (initial epiboly)           [none]            [%epi] 
    ];
% specify initial condition
input=embryo;
% specify what to plot
selector=3; % 2: cell number; 3: epiboly thickness; 4:avg div rate; 5:div rate/dt

%% Model outputs
alpha=input(1);     % alpha         (division rate at R_y)      [/hour]
gamma=input(2);     % gamma         (decay rate)                [1/mm]
Ni=input(3);        % N0            (initial cell number)       [cells]
t0=input(4);        % t0            (start time)                [hours]
dt=input(5);        % dt            (time step)                 [hours]
t1=input(6);        % t1            (end time)                  [hours] 
x0=input(7);        % x0            (yolk radius R_y)           [mm]
dx=input(8);        % dx            (thickness step)            [mm]
density=input(9);   % density       (cell density)              [cells/mm^3]
beta=input(10);     % beta          (epiboly rate)              [%epi/hour]
epi0=input(11);     % ini epi       (initial epiboly)           [%epi]

t=t0:dt:t1; % time vector
data=zeros(5,size(t,2)+1); %storage matrix
    data(1,:)=[t t1+dt]; %time storage
    data(2,1)=Ni; %cell number storage   
    for i=1:size(t,2) %time loop
        epiboly_percentage = epi0+beta*(t(i)-t0);   %assumption linear epiboly progress
        blastoderm_volume = (data(2,i)/density);    %blastoderm cells we have/cell density = volume blastoderm [mm^3]
        spherical_shell_volume = blastoderm_volume/epiboly_percentage;  %fit blastoderm into portion of spherical shell [mm^3]
        xn = ((spherical_shell_volume+4/3*pi*x0^3)/(4/3*pi))^(1/3); %xn (R_e) of spherical shell volume
        thickness = xn-x0;  %subtract yolk radius R_y
        avg_cell_number_change=0;   %clear and initialize per timestep
        for x=x0:dx:xn
            division_rate=alpha*exp((-x+x0)*gamma); %inverse of tau
            data(5,i)=division_rate;
            avg_cell_number_change=avg_cell_number_change+4*pi*(x^2)*dx*division_rate*epiboly_percentage*density;   %integrate over x, number of cells divided per dx shell
        end
            data(2,i+1)=data(2,i)+avg_cell_number_change*dt;    %add previous cells to new cells*time that has passed
            %data(2,i+1)=data(2,i); %no cell division (control)
            data(3,i) = thickness;  %epiboly thickness 
            data(4,i)=(data(2,i+1)/data(2,i))-1;    %avg division rate (% cells)
            
   end
    
%% Plot simulation

plots(data,selector);
hold on;
    
    function plots(mat,ind) 

    plot(mat(1,1:end-1),mat(ind,1:end-1), 'input plot colour/marker');
    ylabel('input y-axis label');
    xlabel('input x-axis label');
    legend('');
    end
 