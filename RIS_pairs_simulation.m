close all
clear
clc

%% switches & options...
postprocessing_only = 0;
use_pml = 0;         % use pml boundaries instead of mur
cal=0; %% if it is zero, no calculation is done.Only the structure is presented

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

%% configuration of RISs %%%%
dim_meta=3; %%% dimension of metasurface-elements per column and row. ONLY ODD NUMBERS!!!

RIS1_x=0; %%RIS1 position%%
RIS1_y=5;
RIS1_z=0;

RIS2_x=0; %%RIS2 position%%
RIS2_y=-5;
RIS2_z=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup FDTD parameter & excitation function
max_timesteps = 500000;
min_decrement = 1e-5; % equivalent to -50 dB
f0 = 8e9; % center frequency
fc = 2e9; % 20 dB corner frequency
FDTD = InitFDTD( 'NrTS', max_timesteps, 'EndCriteria', min_decrement );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
if (use_pml>0)
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % use pml instead of mur
end
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
substrate_srr.epsR = 2.2;

% set the resolution for the finer structures, e.g. the S-SRR's gap
max_res = c0 / (f0 + fc) / sqrt(substrate_srr.epsR) / unit /40;
% set the resolution for the coarser structures, e.g. the surrounding air
coarseResolution = c0/(f0 + fc) / unit / 20;
CSX = InitCSX();

%%  the dimensions of fine tuning procedure
L1=10.9; %length of the outer ring
L2=10.36; %width of the outer ring
G1=1.05; % width of the outer gap
width_outer=1.05; %%width outer patch
width_iner=width_outer;
G3=1.05; %distance between iner and outer
G2=G1; % width of inner gap
srr_thickness=0.15;


%% parameters for load patches
width_patch=1.2;
patch_thickness=srr_thickness;
dx=10;% x-axis distance between elements
dz=10; % z-axis distance between elements

%% substrate setup
distance_x=L2+dx;
distance_z=L1+dz;
substrate_srr.width = dim_meta*distance_x;
substrate_srr.length =dim_meta*distance_z;
substrate_srr.thickness = 2.1;
substrate_srr.cells = 4;
tand = 0.024;
substrate_srr.kappa = tand*2*pi*f0*EPS0*substrate_srr.epsR;

feed.R = 50; % feed resistance

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'patch_ant.xml';
if (postprocessing_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
end


%% create materials
CSX = AddMaterial(CSX,'substrate_srr'); %create substrate
CSX = SetMaterialProperty(CSX,'substrate_srr','Epsilon',substrate_srr.epsR,'Kappa',substrate_srr.kappa);
CSX = AddMetal(CSX,'groundplane'); %create groundplane
CSX = AddMetal( CSX, 'patch');%create patches
CSX = AddMetal( CSX, 'SRR'); %create S-SRRs

%% First RIS
%% create substrate
start = [-substrate_srr.length/2+RIS1_x, RIS1_y,-substrate_srr.width/2+RIS1_z];
stop = [ substrate_srr.length/2+RIS1_x, substrate_srr.thickness+RIS1_y,substrate_srr.width/2+RIS1_z];
CSX = AddBox(CSX,'substrate_srr',1,start,stop);
%% create groundplane
start = [ -substrate_srr.length/2+RIS1_x, RIS1_y+substrate_srr.thickness,-substrate_srr.width/2+RIS1_z];
stop = [substrate_srr.length/2+RIS1_x, RIS1_y+substrate_srr.thickness, substrate_srr.width/2+RIS1_z];
CSX = AddBox(CSX,'groundplane',2,start,stop);
%create unit cells-S-SRRs
number_of_elements=0;

for i=-((dim_meta-1)/2):((dim_meta-1)/2)
    x=i*distance_x;
    for j=-((dim_meta-1)/2):((dim_meta-1)/2)
        number_of_elements=number_of_elements+1;
        z=j*distance_z;
        
        %% outer ring
        start = [L2/2+RIS1_x+x, RIS1_y+srr_thickness/2,-L1/2+RIS1_z+width_outer+z];
        stop = [L2/2-width_outer+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2+RIS1_z-width_outer+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [L2/2+G1/2-width_outer/2+RIS1_x+x, RIS1_y+srr_thickness/2, -L1/2+RIS1_z+z];
        stop = [G1/2+RIS1_x+x, RIS1_y-srr_thickness/2,-L1/2+width_outer+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2-G1/2+width_outer/2+RIS1_x+x, RIS1_y+srr_thickness/2, -L1/2+RIS1_z+z];
        stop = [-G1/2+RIS1_x+x, RIS1_y-srr_thickness/2,-L1/2+width_outer+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+RIS1_x+x, RIS1_y+srr_thickness/2, L1/2+RIS1_z+z];
        stop = [ L2/2+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2-width_outer+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+RIS1_x+x, RIS1_y+srr_thickness/2,-L1/2+RIS1_z+width_outer+z];
        stop = [-L2/2+width_outer+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2+RIS1_z-width_outer+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        
        %% iner ring
        start = [ L2/2-width_iner-G3-width_iner+RIS1_x+x, RIS1_y+srr_thickness/2, -L1/2+width_iner+G3+RIS1_z+z];
        stop = [-L2/2+width_iner+G3+width_iner+RIS1_x+x, RIS1_y-srr_thickness/2,-L1/2+width_iner+G3+width_iner+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-G2/2+RIS1_x+x, RIS1_y+srr_thickness/2,L1/2-width_iner-G3+RIS1_z+z];
        stop = [ -L2/2+width_iner+width_iner+G3+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2-width_iner-G3-width_iner+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [G2/2+RIS1_x+x, RIS1_y+srr_thickness/2,L1/2-width_iner-G3+RIS1_z+z];
        stop = [ L2/2-width_iner-G3-width_iner+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2-width_iner-G3-width_iner+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+width_iner+G3+RIS1_x+x, RIS1_y+srr_thickness/2,-L1/2+width_iner+G3+RIS1_z+z];
        stop = [-L2/2+width_iner+G3+width_iner+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2-width_iner-G3+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [L2/2-width_iner-G3+RIS1_x+x, RIS1_y+srr_thickness/2,-L1/2+width_iner+G3+RIS1_z+z];
        stop = [L2/2-width_iner-G3-width_iner+RIS1_x+x, RIS1_y-srr_thickness/2,L1/2-width_iner-G3+RIS1_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        
        % create patches
        if i~=(dim_meta-1)/2
            start = [L2/2+RIS1_x+i*(dx+L2), RIS1_y+patch_thickness/2, L1/4+RIS1_z+z];
            stop = [ L2/2+RIS1_x+dx+i*(dx+L2), RIS1_y-srr_thickness/2,L1/4-width_patch+RIS1_z+z];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
            start = [L2/2+RIS1_x+i*(dx+L2), RIS1_y+patch_thickness/2, -L1/4+RIS1_z+z];
            stop = [ L2/2+RIS1_x+dx+i*(dx+L2), RIS1_y-patch_thickness/2,-L1/4+width_patch+RIS1_z+z];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
        end
        
        if j~=(dim_meta-1)/2
            start = [L2/4+RIS1_x+x, RIS1_y+patch_thickness/2,L1/2+RIS1_z+j*(L1+dz)];
            stop = [L2/4+width_patch+RIS1_x+x, RIS1_y-patch_thickness/2,L1/2+RIS1_z+dz+j*(L1+dz)];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
            start = [-L2/4+RIS1_x+x, RIS1_y+patch_thickness/2,L1/2+RIS1_z+j*(L1+dz)];
            stop = [-L2/4-width_patch+RIS1_x+x, RIS1_y-patch_thickness/2,L1/2+RIS1_z+dz+j*(L1+dz)];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
        end
        
        for k=number_of_elements
            start = [-G1/2+RIS1_x+x, RIS1_y, -L1/2+RIS1_z+z];
            stop  = [G1/2+RIS1_x+x, substrate_srr.thickness+RIS1_y, -L1/2+width_outer+RIS1_z+z];
            [CSX, port{k}] = AddLumpedPort(CSX, 1 ,k ,feed.R, start, stop, [0 1 0], true);
            
        end
    end
    
end


%% Second Metasurface
%% create substrate
start = [-substrate_srr.length/2+RIS2_x, RIS2_y,-substrate_srr.width/2+RIS2_z];
stop = [ substrate_srr.length/2+RIS2_x, substrate_srr.thickness+RIS2_y,substrate_srr.width/2+RIS2_z];
CSX = AddBox(CSX,'substrate_srr',1,start,stop);
%% create groundplane
start = [ -substrate_srr.length/2+RIS2_x, RIS2_y,-substrate_srr.width/2+RIS2_z];
stop = [substrate_srr.length/2+RIS2_x, RIS2_y, substrate_srr.width/2+RIS2_z];
CSX = AddBox(CSX,'groundplane',2,start,stop);

number_of_elements_2=0;
x=0;
z=0;

for i=-((dim_meta-1)/2):((dim_meta-1)/2)
    x=i*distance_x;
    for j=-((dim_meta-1)/2):((dim_meta-1)/2)
        number_of_elements_2=number_of_elements_2+1;
        z=j*distance_z;
        %% outer ring
        CSX = AddMetal( CSX, 'SRR');
        start = [L2/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+RIS2_z+z+width_outer];
        stop = [L2/2-width_outer+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+z-width_outer];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [L2/2+G1/2-width_outer/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+RIS2_z+z];
        stop = [G1/2+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_outer+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2-G1/2+width_outer/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+RIS2_z+z];
        stop = [-G1/2+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_outer+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness, L1/2+RIS2_z+z];
        stop = [ L2/2+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_outer+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+RIS2_z+z+width_outer];
        stop = [-L2/2+width_outer+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+z-width_outer];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        
        
        %% iner ring
        start = [ L2/2-width_iner-G3-width_iner+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+width_iner+G3+RIS2_z+z];
        stop = [-L2/2+width_iner+G3+width_iner+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+width_iner+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-G2/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+RIS2_z+z];
        stop = [ -L2/2+width_iner+width_iner+G3+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3-width_iner+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [G2/2+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+RIS2_z+z];
        stop = [ L2/2-width_iner-G3-width_iner+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3-width_iner+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [-L2/2+width_iner+G3+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+RIS2_z+z];
        stop = [-L2/2+width_iner+G3+width_iner+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        start = [L2/2-width_iner-G3+RIS2_x+x, RIS2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+RIS2_z+z];
        stop = [L2/2-width_iner-G3-width_iner+RIS2_x+x, RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+RIS2_z+z];
        CSX = AddBox(CSX,'SRR',10,start,stop);
        
        if i~=(dim_meta-1)/2
            start = [L2/2+RIS2_x+i*(dx+L2), RIS2_y+patch_thickness/2+substrate_srr.thickness, L1/4+RIS2_z+z];
            stop = [ L2/2+RIS2_x+dx+i*(dx+L2), RIS2_y-srr_thickness/2+substrate_srr.thickness,L1/4-width_patch+RIS2_z+z];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
            start = [L2/2+RIS2_x+i*(dx+L2), RIS2_y+patch_thickness/2+substrate_srr.thickness, -L1/4+RIS2_z+z];
            stop = [ L2/2+RIS2_x+dx+i*(dx+L2), RIS2_y-srr_thickness/2+substrate_srr.thickness,-L1/4+width_patch+RIS2_z+z];
            CSX = AddBox(CSX,'patch',4,start,stop);
        end
        
        if j~=(dim_meta-1)/2
            start = [L2/4+RIS2_x+x, RIS2_y+patch_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+j*(L1+dz)];
            stop = [L2/4+width_patch+RIS2_x+x, RIS2_y-patch_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+dz+j*(L1+dz)];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
            start = [-L2/4+RIS2_x+x, RIS2_y+patch_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+j*(L1+dz)];
            stop = [-L2/4-width_patch+RIS2_x+x, RIS2_y-patch_thickness/2+substrate_srr.thickness,L1/2+RIS2_z+dz+j*(L1+dz)];
            CSX = AddBox(CSX,'patch',4,start,stop);
            
            for k=number_of_elements_2
                start = [-G1/2+RIS2_x+x, RIS2_y, -L1/2+RIS2_z+z];
                stop  = [G1/2+RIS2_x+x, RIS2_y+substrate_srr.thickness, -L1/2+width_outer+RIS2_z+z];
                [CSX, port{k+number_of_elements}] = AddLumpedPort(CSX, 5 ,k+number_of_elements ,feed.R, start, stop, [0 1 0], false);
                
            end
            
        end
        
    end
end

% setup a mesh
mesh.x = [];
mesh.y = [];

% two mesh lines for the metal coatings of teh substrate
mesh.z = linspace(-substrate_srr.thickness+RIS1_y, 0, substrate_srr.cells +1);

% find optimal mesh lines for the patch and ground, not yes the microstrip line
mesh = DetectEdges(CSX, mesh, 'SetProperty',{'groundplane', 'patch'}, '2D_Metal_Edge_Res', max_res/2);

%replace gap mesh lines which are too close by a single mesh line in each
%axis
tooclose = find (diff(mesh.y) < max_res/2);
if ~isempty(tooclose)
    mesh.y(tooclose) = (mesh.y(tooclose) + mesh.y(tooclose+1))/2;
    mesh.y(tooclose + 1) = [];
end

tooclose = find (diff(mesh.x) < max_res/2);
if ~isempty(tooclose)
    mesh.x(tooclose) = (mesh.x(tooclose) + mesh.x(tooclose+1))/2;
    mesh.x(tooclose + 1) = [];
end

tooclose = find (diff(mesh.z) < max_res/2);
if ~isempty(tooclose)
    mesh.z(tooclose) = (mesh.z(tooclose) + mesh.z(tooclose+1))/2;
    mesh.z(tooclose + 1) = [];
end

% store the microstrip  edges in a temporary variable
meshline = DetectEdges(CSX, [], 'SetProperty', 'SRR', '2D_Metal_Edge_Res', max_res/2);
% as well as the edges of the substrate (without 1/3 - 2/3 rule)
meshsubstrate = DetectEdges(CSX, [], 'SetProperty', 'substrate_srr');
% add only the x mesh lines of the microstrip
mesh.x = [mesh.x meshline.x];
% and only the top of the substrate, the other edges are covered by the ground plane
mesh.y = [mesh.y, meshsubstrate.y]; % top of substrate

% for now we have only the edges, now calculate mesh lines inbetween
mesh = SmoothMesh(mesh, max_res);

% add the outer boundary
mesh.x = [mesh.x,-dim_meta*20, +dim_meta*20];
mesh.y = [mesh.y, -dim_meta*10, dim_meta*10];
mesh.z = [mesh.z, -dim_meta*20, +dim_meta*20];

% add coarse mesh lines for the free space
%mesh = SmoothMesh(mesh, coarseResolution,'algorithm',[1 3]);
mesh = SmoothMesh(mesh, coarseResolution);
% define the grid
CSX = DefineRectGrid( CSX, unit, mesh);

start = [mesh.x(2)     mesh.y(2)     mesh.z(2)];
stop  = [mesh.x(end-1) mesh.y(end-1) mesh.z(end-1)];
CSX = AddDump(CSX,'Ef', 'DumpType', 10, 'Frequency',(f0));
CSX = AddBox(CSX,'Ef',10,start, stop); %assign box

%% add a nf2ff calc box,size is 3 cells away from bound cond
if (use_pml == 0)
    start = [mesh.x(4) mesh.y(4) mesh.z(4)];
    stop = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
else
    start = [mesh.x(12) mesh.y(12) mesh.z(12)];
    stop = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
end
[CSX, nf2ff] = CreateNF2FFBox(CSX,'nf2ff',start,stop);

%% prepare and run simulation folder
Sim_Path = 'tmp_ris_pairs';
Sim_CSX = 'ris_pairs.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);

if (cal==1)
    %% run openEMS
    RunOpenEMS(Sim_Path,Sim_CSX);
    
    %% POST PROCESSING AND DO THE PLOTS
    %Postprocessing &Plots
    freq = linspace(f0-fc,f0+fc,1001);
    
    
    for i=1:2*dim_meta^2
        port{i} = calcPort(port{i},Sim_Path,freq);
        Zin{i} = port{i}.uf.tot ./ port{i}.if.tot;
        P_in{i} = real(0.5 * port{i}.uf.tot .* conj( port{i}.if.tot ));
        Pincoming{i}=port{i}.P_inc;
        Preflected{i}=port{i}.P_ref;
        Paccepted{i}=port{i}.P_acc;   %%%incoming minus reflected, may be negative for passive port
    end
    
    
    for i=1:2*dim_meta^2
        for j=1:2*dim_meta^2
            s{j,i} = port{j}.uf.ref ./ port{i}.uf.inc;
        end
    end
    
    
    for i=1:2*dim_meta^2
        for j=1:2*dim_meta^2
            if j==i && i<=dim_meta^2
                f_res{i}=find(s{j,i}==min(s{j,i}));
            elseif i>dim_meta^2
                f_coupling{j,i}=find(s{j,i}==max(s{j,i}));
            end
        end
    end
    
    
    for i=1:dim_meta^2
        figure
        for j=1:2*dim_meta^2
            if i==j
                plot(freq/1e6, 20*log10(abs(s{j,i})), 'b-', 'Linewidth', 4 );
                hold on
            elseif j>dim_meta^2 && i~=j
                plot( freq/1e6, 20*log10(abs(s{j,i})), 'k-', 'Linewidth', 0.25 );
                hold on
            elseif j<dim_meta^2   && i~=j
                plot( freq/1e6, 20*log10(abs(s{j,i})), 'r-', 'Linewidth', 0.25 );
                hold on
            end
            grid on
            title( 'Simulated S-parameters' );
            xlabel( 'frequency f / MHz' );
            ylabel( 'S-parameters(dB)' );
        end
        drawnow
    end
    
end
