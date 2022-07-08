close all
clear
clc

%% switches & options...
postprocessing_only = 0;
use_pml = 0;         % use pml boundaries instead of mur
cal=1;
active=1; % declares what S-SRR is active and what passive

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

%% setup FDTD parameter & excitation function
max_timesteps = 150000;
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
substrate.epsR = 2.2;

% set the resolution for the finer structures, e.g. the S-SRR's gap
max_res = c0 / (f0 + fc) / sqrt(substrate.epsR) / unit /40;% cell size: lambda/40
% set the resolution for the coarser structures, e.g. the surrounding air
coarseResolution = c0/(f0 + fc) / unit / 20; % cell size: lambda/20

CSX = InitCSX();

%% S-SRR's dimensions
L1=10.9; %length of the outer ring
L2=10.36; %width of the outer ring
G1=1.05; % width of the outer gap
width_outer=1.05; %%width outer patch
width_iner=width_outer; 
G3=1.05; %distance between iner and outer
G2=G1; % width of inner gap
srr_thickness=0.15;

%substrate setup
dx=1.85; 
dz=1.85+L2-L1;
substrate_srr.width = L2+dx;
substrate_srr.length =L1+dz;
substrate_srr.thickness = 4.8;
substrate_srr.cells = 4;
tand = 0.024;
substrate_srr.kappa = tand*2*pi*f0*EPS0*substrate.epsR;
feed.R = 50;% feed resistance

%% S-SRRs positions %%
ssrr_x=0; 
ssrr_y=5; 
ssrr_z=0;

ssrr2_x=0; 
ssrr2_y=-5; 
ssrr2_z=0; 


%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'patch_ant.xml';
if (postprocessing_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
end


%% create substrate
CSX = AddMaterial(CSX,'substrate_srr');
CSX = SetMaterialProperty(CSX,'substrate_srr','Epsilon',substrate.epsR,'Kappa',substrate_srr.kappa);
start = [-substrate_srr.length/2+ssrr_x, ssrr_y,-substrate_srr.width/2+ssrr_z];
stop = [ substrate_srr.length/2+ssrr_x, substrate_srr.thickness+ssrr_y,substrate_srr.width/2+ssrr_z];
CSX = AddBox(CSX,'substrate_srr',1,start,stop);

%% create groundplane
CSX = AddMetal(CSX,'groundplane'); %create a PEC
start = [ -substrate_srr.length/2+ssrr_x,ssrr_y+substrate_srr.thickness,-substrate_srr.width/2+ssrr_z];
stop = [substrate_srr.length/2+ssrr_x, ssrr_y+substrate_srr.thickness, substrate_srr.width/2+ssrr_z];
CSX = AddBox(CSX,'groundplane',10,start,stop);


%% outer ring
     CSX = AddMetal( CSX, 'SRR');
     start = [L2/2+ssrr_x, ssrr_y+srr_thickness/2,-L1/2+ssrr_z+width_outer]; 
     stop = [L2/2-width_outer+ssrr_x, ssrr_y-srr_thickness/2,L1/2+ssrr_z-width_outer];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [L2/2+G1/2-width_outer/2+ssrr_x, ssrr_y+srr_thickness/2, -L1/2+ssrr_z]; 
     stop = [G1/2+ssrr_x, ssrr_y-srr_thickness/2,-L1/2+width_outer+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2-G1/2+width_outer/2+ssrr_x, ssrr_y+srr_thickness/2, -L1/2+ssrr_z]; 
     stop = [-G1/2+ssrr_x, ssrr_y-srr_thickness/2,-L1/2+width_outer+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+ssrr_x, ssrr_y+srr_thickness/2, L1/2+ssrr_z];
     stop = [ L2/2+ssrr_x, ssrr_y-srr_thickness/2,L1/2-width_outer+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+ssrr_x, ssrr_y+srr_thickness/2,-L1/2+ssrr_z+width_outer]; 
     stop = [-L2/2+width_outer+ssrr_x, ssrr_y-srr_thickness/2,L1/2+ssrr_z-width_outer]; 
     CSX = AddBox(CSX,'SRR',10,start,stop);

 %%iner ring
     start = [ L2/2-width_iner-G3-width_iner+ssrr_x, ssrr_y+srr_thickness/2, -L1/2+width_iner+G3+ssrr_z];
     stop = [-L2/2+width_iner+G3+width_iner+ssrr_x, ssrr_y-srr_thickness/2,-L1/2+width_iner+G3+width_iner+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-G2/2+ssrr_x, ssrr_y+srr_thickness/2,L1/2-width_iner-G3+ssrr_z];
     stop = [ -L2/2+width_iner+width_iner+G3+ssrr_x, ssrr_y-srr_thickness/2,L1/2-width_iner-G3-width_iner+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [G2/2+ssrr_x, ssrr_y+srr_thickness/2,L1/2-width_iner-G3+ssrr_z];
     stop = [ L2/2-width_iner-G3-width_iner+ssrr_x, ssrr_y-srr_thickness/2,L1/2-width_iner-G3-width_iner+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+width_iner+G3+ssrr_x, ssrr_y+srr_thickness/2,-L1/2+width_iner+G3+ssrr_z];
     stop = [-L2/2+width_iner+G3+width_iner+ssrr_x, ssrr_y-srr_thickness/2,L1/2-width_iner-G3+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [L2/2-width_iner-G3+ssrr_x, ssrr_y+srr_thickness/2,-L1/2+width_iner+G3+ssrr_z];
     stop = [L2/2-width_iner-G3-width_iner+ssrr_x, ssrr_y-srr_thickness/2,L1/2-width_iner-G3+ssrr_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
  
     
%%% S-SRR1 %%%%%%

%% create substrate
CSX = SetMaterialProperty(CSX,'substrate_srr','Epsilon',substrate.epsR,'Kappa',substrate_srr.kappa);
start = [-substrate_srr.length/2+ssrr2_x, ssrr2_y,-substrate_srr.width/2+ssrr2_z];
stop = [ substrate_srr.length/2+ssrr2_x, substrate_srr.thickness+ssrr2_y,substrate_srr.width/2+ssrr2_z];
CSX = AddBox(CSX,'substrate_srr',1,start,stop);

%% create groundplane
CSX = AddMetal(CSX,'groundplane'); %create a PEC
start = [ -substrate_srr.length/2+ssrr2_x,ssrr2_y,-substrate_srr.width/2+ssrr2_z];
stop = [substrate_srr.length/2+ssrr2_x, ssrr2_y, substrate_srr.width/2+ssrr2_z];
CSX = AddBox(CSX,'groundplane',10,start,stop);

    %% outer ring
     CSX = AddMetal( CSX, 'SRR');
     start = [L2/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+ssrr2_z+width_outer]; 
     stop = [L2/2-width_outer+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2+ssrr2_z-width_outer];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [L2/2+G1/2-width_outer/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+ssrr2_z];  
     stop = [G1/2+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_outer+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2-G1/2+width_outer/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+ssrr2_z];
     stop = [-G1/2+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_outer+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness, L1/2+ssrr2_z];
     stop = [ L2/2+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_outer+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+ssrr2_z+width_outer]; 
     stop = [-L2/2+width_outer+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2+ssrr2_z-width_outer]; 
     CSX = AddBox(CSX,'SRR',10,start,stop);

 %%iner ring
     start = [ L2/2-width_iner-G3-width_iner+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness, -L1/2+width_iner+G3+ssrr2_z];
     stop = [-L2/2+width_iner+G3+width_iner+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+width_iner+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-G2/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+ssrr2_z];
     stop = [ -L2/2+width_iner+width_iner+G3+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3-width_iner+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [G2/2+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+ssrr2_z];
     stop = [ L2/2-width_iner-G3-width_iner+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3-width_iner+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [-L2/2+width_iner+G3+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+ssrr2_z];
     stop = [-L2/2+width_iner+G3+width_iner+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
     start = [L2/2-width_iner-G3+ssrr2_x, ssrr2_y+srr_thickness/2+substrate_srr.thickness,-L1/2+width_iner+G3+ssrr2_z];
     stop = [L2/2-width_iner-G3-width_iner+ssrr2_x, ssrr2_y-srr_thickness/2+substrate_srr.thickness,L1/2-width_iner-G3+ssrr2_z];
     CSX = AddBox(CSX,'SRR',10,start,stop);
    
 if (active==1)
     start = [-G1/2+ssrr_x, ssrr_y, -L1/2+ssrr_z];
     stop  = [G1/2+ssrr_x, substrate_srr.thickness+ssrr_y, -L1/2+width_outer+ssrr_z];
     [CSX,port{1}] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 1 0], true); 
     start = [-G1/2+ssrr2_x, ssrr2_y, -L1/2+ssrr2_z];
     stop  = [G1/2+ssrr2_x, substrate_srr.thickness+ssrr2_y, -L1/2+width_outer+ssrr2_z];
     [CSX,port{2}] = AddLumpedPort(CSX, 5 ,2 ,feed.R, start, stop, [0 1 0], false); 
 elseif (active==2)
     start = [-G1/2+ssrr_x, ssrr_y, -L1/2+ssrr_z];
     stop  = [G1/2+ssrr_x, substrate_srr.thickness+ssrr_y, -L1/2+width_outer+ssrr_z];
     [CSX, port{1}] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 1 0], false); 
     start = [-G1/2+ssrr2_x, ssrr2_y, -L1/2+ssrr2_z];
     stop  = [G1/2+ssrr2_x, substrate_srr.thickness+ssrr2_y, -L1/2+width_outer+ssrr2_z];
     [CSX,port{2}] = AddLumpedPort(CSX, 5 ,2 ,feed.R, start, stop, [0 1 0], true); 
 end  

% setup a mesh
mesh.x = [];
mesh.y = [];

% two mesh lines for the metal coatings of teh substrate
mesh.z = linspace(-substrate_srr.thickness+ssrr_y, 0, substrate_srr.cells +1);

% find optimal mesh lines for the patch and ground, not yes the microstrip line
mesh = DetectEdges(CSX, mesh, 'SetProperty','groundplane', '2D_Metal_Edge_Res', max_res/2);
%replace gap mesh lines which are too close by a single mesh line
tooclose = find (diff(mesh.y) < max_res/4);
if ~isempty(tooclose)
  mesh.y(tooclose) = (mesh.y(tooclose) + mesh.y(tooclose+1))/2;
  mesh.y(tooclose + 1) = [];
end

tooclose = find (diff(mesh.x) < max_res/4);
if ~isempty(tooclose)
  mesh.x(tooclose) = (mesh.x(tooclose) + mesh.x(tooclose+1))/2;
  mesh.x(tooclose + 1) = [];
end

tooclose = find (diff(mesh.z) < max_res/4);
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
mesh.x = [mesh.x -60, 60];
mesh.y = [mesh.y, -60, 65];
mesh.z = [mesh.z, -45, 45];

% add coarse mesh lines for the free space
mesh = SmoothMesh(mesh, coarseResolution);

% define the grid
CSX = DefineRectGrid( CSX, unit, mesh);


%% add a nf2ff calc box,size is 3 cells away from bound cond
if (use_pml == 0)
    start = [mesh.x(4) mesh.y(4) mesh.z(4)];
    stop = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
else
    start = [mesh.x(12) mesh.y(12) mesh.z(12)];
    stop = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
end
[CSX, nf2ff] = CreateNF2FFBox(CSX,'nf2ff',start,stop);


%% Paraview
CSX = AddDump(CSX,'Ef','DumpType',10,'Frequency',(5e9));
CSX = AddBox(CSX,'Ef',10,[-substrate_srr.width -substrate_srr.length -10*substrate_srr.thickness],[substrate_srr.width substrate_srr.length 10*substrate_srr.thickness]); %assign box

%% prepare and run simulation folder
Sim_Path = 'tmp_ssrr';
Sim_CSX = 'ssrr.xml';

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

port{1} = calcPort(port{1},Sim_Path,freq);
port{2} = calcPort(port{2},Sim_Path,freq);

Pincoming_1=port{1}.P_inc;
Preflected_1=port{1}.P_ref;
Paccepted_1=port{1}.P_acc;   %%%incoming minus reflected, may be negative for passive ports
Pincoming_2=port{2}.P_inc;
Preflected_2=port{2}.P_ref;
Paccepted_2=port{2}.P_acc;   %%%incoming minus reflected, may be negative for passive ports

Zin1 = port{1}.uf.tot ./ port{1}.if.tot;
Zin2 = port{2}.uf.tot ./ port{2}.if.tot;

if (active == 1)
    s11 = port{1}.uf.ref ./ port{1}.uf.inc;
    s21 = port{2}.uf.ref ./ port{1}.uf.inc;
elseif (active == 2)
    s22 = port{2}.uf.ref ./ port{2}.uf.inc;
    s12 = port{1}.uf.ref ./ port{2}.uf.inc;
end    
if (active == 1)
    % plot feed point impedance
    figure
    plot( freq/1e9, real(Zin1), 'k-', 'Linewidth', 2 );
    hold on
    grid on
    plot( freq/1e9, imag(Zin1), 'r--', 'Linewidth', 2 );
    title( 'Feed Point Impedance' );
    xlabel( 'frequency f /  GHz' );
    ylabel( 'impedance Z_{in} / Ohm' );
    legend( 'real', 'imag' );

    % plot refl coefficient S11,S-parameters S21
    figure
    plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
    hold on
    grid on
    plot( freq/1e9, 20*log10(abs(s21)), 'b-', 'Linewidth', 2 );
    title( 'Simulated S-parameters' );
    xlabel( 'frequency f / GHz' );
    ylabel( 'S-parameters(dB)' );
    legend('S11','S21');
    drawnow
    
end


if (active == 2)
    % plot feed point impedance
    figure
    plot( freq/1e9, real(Zin2), 'k-', 'Linewidth', 2 );
    hold on
    grid on
    plot( freq/1e9, imag(Zin2), 'r--', 'Linewidth', 2 );
    title( 'Feed Point Impedance' );
    xlabel( 'frequency f / MHz' );
    ylabel( 'impedance Z_{in} / Ohm' );
    legend( 'real', 'imag' );

    % plot reflection coefficient S11,S-parameters S21,S31,S41
    figure
    plot( freq/1e9, 20*log10(abs(s22)), 'k-', 'Linewidth', 2 );
    hold on
    grid on
    plot( freq/1e9, 20*log10(abs(s12)), 'b-', 'Linewidth', 2 );
    title( 'Simulated S-parameters' );
    xlabel( 'frequency f / GHz' );
    ylabel( 'S-parameters(dB)' );
    legend('S22','S12');
    
    drawnow
    
    
    
end


   % plot reflection coefficient S11,S-parameters S21,S31,S41

    figure
    plot( freq/1e9,Pincoming_1, 'k-', 'Linewidth', 2 );
    hold on
    plot( freq/1e9,Preflected_1, 'r-', 'Linewidth', 2 );
    hold on
    plot( freq/1e9,Paccepted_1, 'b-', 'Linewidth', 2 );
    grid on
    title( 'Power' );
    xlabel( 'frequency f / GHz' );
    ylabel( 'Power' );
    legend('Power Incoming 1','Power reflected 1', 'Power accepted 1');
    drawnow
    
    figure
    plot( freq/1e9,Pincoming_2, 'k-', 'Linewidth', 2 );
    hold on
    plot( freq/1e9,Preflected_2, 'r-', 'Linewidth', 2 );
    hold on
    plot( freq/1e9,Paccepted_2, 'b-', 'Linewidth', 2 );
    grid on
    title( 'Power' );
    xlabel( 'frequency f / GHz' );
    ylabel( 'Power' );
    legend('Power Incoming 2','Power reflected 2', 'Power accepted 2');
    drawnow


%% NFFF contour plots
%find resonance frequncy from s11
if (active == 1)
    f_res_ind = find(s11==min(s11));
    f_res = freq(f_res_ind);
elseif (active == 2)
    f_res_ind = find(s22==min(s22));
    f_res = freq(f_res_ind);
end

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, (-180:2:180)*pi/180, [0 90]*pi/180, 'Mode', 1);


% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );

if (active == 1)
    disp( ['efficiency(ant1) @30GHz = ' num2str(100*nf2ff.Prad./port{1}.P_inc(f_res_ind)) ' %']);
elseif (active == 2)
    disp( ['efficiency(ant2) @30GHz = ' num2str(100*nf2ff.Prad./port{2}.P_inc(f_res_ind)) ' %']);
end

% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'normalize',1)

% log-scale directivity plot
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2])
% conventional plot approach
% plot( nf2ff.theta*180/pi, 20*log10(nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:)))+10*log10(nf2ff.Dmax));

drawnow

% Show 3D pattern
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);


E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);
end
