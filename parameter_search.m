clear;
%% Experimental data 
experimental_time = [2 3 4 5 6];
experimental_thickness_control= [0.1316,0.0607,0.0446,0.0411,0.0431];
experimental_thickness_big = [0.29,0.1178,0.0718,0.0554,0.0529];
%% Choose min and max values for each parameter
%parameter         abbrev       min    max   step

%division rate (div_rate)
alpha_min=0.001;
alpha_max=2;
alpha_steps=50;
div_rate = linspace(alpha_min,alpha_max,alpha_steps);
%decay rate (decay_rate)
gamma_min=1;
gamma_max=40;
gamma_steps=50;
decay_rate = linspace(gamma_min,gamma_max,gamma_steps);
%initial cell number (N0)
N0_min=2000;
N0_max=4000;
N0_steps=5;
ini_N0 = linspace(N0_min,N0_max,N0_steps);

%% Store parameter vectors
test_matrix=transpose(combvec(div_rate,decay_rate,ini_N0));
%all # combinations possible from 3 parameters, based on #steps/parameter

results_matrix = [test_matrix zeros(length(test_matrix),2)];
%save parameter combination used/iteration + 2 extra tables for grading

%% Run model script for total possible # recombinations 
for g = 1:length(test_matrix)
    disp(g*100/length(test_matrix));
[embryo]=[ ...
    test_matrix(g,1),   ... %alpha          (division rate)             [1/time]          [1/hours]
    test_matrix(g,2),   ... %gamma          (decay rate)                [1/length]        [1/mm]
    test_matrix(g,3),   ... % N0            (initial cell number)       [none]            [cells]
    2,                  ... % t0            (start time)                [time]            [hours]
    0.1,                ... % dt            (time step)                 [time]            [hours]
    6,                  ... % t1            (end time)                  [time]            [hours]
    0.33,               ... % x0            (yolk radius)               [length]          [mm]
    0.003,              ... % dx            (thickness step)            [length]          [mm]
    5.7e+4,             ... % density       (cell density)              [1/volumen]       [cells/mm^3]
    0.2,                ... % beta          (epiboly rate)              [1/time]          [%epi/hour]
    0.2                 ... % ini epi       (initial epiboly)           [none]            [%epi] 
    ];
% specify initial condition
input=embryo;

%% Model output/iteration
alpha=input(1); % alpha      (division rate)             [/hour]
gamma=input(2); % gamma      (decay rate)                [1/mm]
N0=input(3);    % N0         (initial cell number)       [cells]
t0=input(4);    % t0         (start time)                [hours]
dt=input(5);    % dt         (time step)                 [hours]
t1=input(6);    % t1         (end time)                  [hours] 
x0=input(7);    % x0         (yolk radius)               [mm]
dx=input(8);    % dx         (thickness step)            [mm]
density=input(9); % density  (cell density)              [cells/mm^3]
beta=input(10); % beta       (epiboly rate)              [%epi/hour]
epi0=input(11); %ini epi     (initial epiboly)           [%epi]

t=t0:dt:t1; % time vector
data=zeros(5,size(t,2)+1); %storage matrix
    data(1,:)=[t t1+dt]; %time storage
    data(2,1)=N0; %cell number storage   
    for i=1:size(t,2) %time loop
        epiboly_percentage = epi0+beta*(t(i)-t0); %assumption linear epiboly progress
        blastoderm_volume = (data(2,i)/density); %cells we have/density=volume blastoderm [mm^3]
        spherical_shell_volume = blastoderm_volume/epiboly_percentage; %fit blastoderm into portion of a spherical shell [mm^3]
        xn = ((spherical_shell_volume+4/3*pi*x0^3)/(4/3*pi))^(1/3); %xn from spherical shell volume
        thickness = xn-x0; %subtract yolk radius
        avg_cell_number_change=0; %clear and initialize per timestep
        for x=x0:dx:xn
            division_rate=alpha*exp((-x+x0)*gamma); %inverse of tau
            avg_cell_number_change=avg_cell_number_change+4*pi*(x^2)*dx*division_rate*epiboly_percentage*density;   %integrate over x, number of cells divided per dx shell
        end
            data(2,i+1)=data(2,i)+avg_cell_number_change*dt;    %add previous cells to new cells*time that has passed
            data(3,i) = thickness; %epiboly thickness
            data(4,i)=(data(2,i+1)/data(2,i))-1; %avg division rate
            
   end
    
    %% Grading each combination (control embryo)
    % select similar time points as experimental data
    model_thickness_control=interp1(data(1,:),data(3,:),experimental_time);
    
    grade_control=(experimental_thickness_control-model_thickness_control).^2;

    results_matrix(g,4) = mean(grade_control);
    results_matrix(g,5) = std(grade_control);
end

%% 3D scatter plots: comparing grades between 3 parameters
x_by = 1;   %x-axis 1st parameter
y_by= 2;    %y-axis 2nd parameter
colour_by = 3;  %colour 3rd parameter
colour_list=['r','b','g','y','m','c','k']; 

colour_by_values = unique(results_matrix(:,colour_by)); %colour/unique parameter value
label_list = ["alpha" "gamma" "N0"]; %x/y/z axis labels

for i=1:length(colour_by_values) 
    trimmed_results_matrix=results_matrix(results_matrix(:,colour_by)==colour_by_values(i),:);  %colour by 3rd parameter from list of colours
    scatter3(trimmed_results_matrix(:,x_by),trimmed_results_matrix(:,y_by),trimmed_results_matrix(:,5),colour_list(i)); %3D scatter plot all 3 parameters
    legend(num2str(colour_by_values));  %legend 3rd parameter
    
    xlabel(label_list(x_by));
    ylabel(label_list(y_by));
    hold on;
end
