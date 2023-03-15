function obj = romsComp(obj,file,plotchoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to automatically generate comparison plots between ROMS simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% - obj = romsComp(obj,file1,file2,plotchoice)
%
% Inputs:
% - obj   = roms objects (2) that are already initialized
%			... comparisons will be obj(1) - obj(2)
% - file  = cell array of files to load  
% -         ... file{1} = obj(1) files etc.
% - plotchoice (see below) 
%
%	------------------------------
%	---- PHYSICAL DIAGNOSTICS ----
%	------------------------------
%	1  == 3d physical diagnostics
%	2  == latitude physical slices
%	3  == surface physical fields
%	4  == longitude physical slices
%	5  == u velocity at latitude
%	-------------------------
%	---- BGC DIAGNOSTICS ----
%	-------------------------
%	6  == 3d bgc diagnostics
%	7  == latitude bgc tracer slices
%	8  == surface bgc fields
%	9  == longitude bgc tracer slices
%	10 == surface chlA
%	11 == OMZ thickness
%	12 == Integrated rates 
%	13 == Integrated tracers
%	14 == Depth slice map
%	15 == latitude bgc rate slices
%	16 == longitude bgc rate slices
%   17 == latitude N-loss
%   18 == longitude N-loss
%       19 == Fe depth slices  
%       20 == POC flux comparisons at 75m   
%	------------------------------
%	------ OTHER DIAGNOSTICS -----
%	------------------------------
%
% Example:
% - obj = romsComp(obj,{[1],[1]},[2:11]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle
addpath /data/project1/demccoy/ROMS/validation/n2o

% Clear objects
for i = 1:length(obj)
	obj(i) = clearROMS(obj(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get plotchoice
tmpchoice = zeros(1,100);
tmpchoice(plotchoice) = 1;
plotchoice = tmpchoice;

% PHYSICAL PLOT options
% (1) PHYS: 3D zslices
pltcnt = 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1];
	plots(pltcnt).vars  = {'temp','salt','sigma'};
	plots(pltcnt).cmaps = {'thermal','haline','dense'};
	plots(pltcnt).zdeps = [0 50 150 300];
	plots(pltcnt).levs  = {linspace(3,30,40),linspace(3,30,40),linspace(0,25,40),linspace(2,18,40);...
						   linspace(30,37,40),linspace(30,37,40),linspace(33.4,36.5,40),linspace(33.5,35.5,40);...
						   linspace(19,27,40),linspace(19,27,40),linspace(24,28,40),linspace(26.5,28.6,40)};
	plots(pltcnt).dlevs = {linspace(-3,3,40),linspace(-0.6,0.6,40),linspace(-0.6,0.6,40)};

% (2) PHYS: Latitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1];
	plots(pltcnt).vars   = {'temp','salt','sigma'};
	plots(pltcnt).cmaps  = {'thermal','haline','dense'};
	plots(pltcnt).choice = 'lat'; 
	plots(pltcnt).lats   = 0;
	plots(pltcnt).xlims  = [140 260];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).levs   = {linspace(3,30,40),linspace(34.5,35.5,40),linspace(22,32,40)};
	plots(pltcnt).dlevs  = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (3) PHYS: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1];
	plots(pltcnt).vars  = {'MLD','SSH'};
	plots(pltcnt).cmaps = {'deep','deep'};
	plots(pltcnt).bal   = [2 2];
	plots(pltcnt).levs  = {linspace(0,150,40),linspace(0,1.5,40)};
	plots(pltcnt).dlevs = {linspace(-60,60,40),linspace(-0.3,0.3,40)};

% (4) PHYS: Longitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1];
	plots(pltcnt).vars   = {'temp','salt','sigma'};
	plots(pltcnt).cmaps  = {'thermal','haline','dense'};
	plots(pltcnt).choice = 'lon';
	plots(pltcnt).lons   = [210 255 272];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).xlims  = [-30 50;-40 20;-40 10];
	plots(pltcnt).levs   = {linspace(3,30,40),linspace(3,30,40),linspace(3,30,40);...
						    linspace(33,36.5,40),linspace(33.5,36.5,40),linspace(33.5,36.5,40);...
						    linspace(22,32,40),linspace(22,32,40),linspace(22,32,40)};
	plots(pltcnt).dlevs  = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (5) PHYS: Zonal equator velocity
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).zlims = [50 500];
	plots(pltcnt).xlims = [-20 20];
	plots(pltcnt).vars  = {'u'};
	plots(pltcnt).cmaps = {'balance'};
	plots(pltcnt).levs  = {linspace(-0.3,0.3,40)};
	plots(pltcnt).dlevs = {linspace(-0.3,0.3,40)};

% BIOGEOCHEMICAL PLOT options
% (6) BGC: 3D tracer zslices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1 1 1 1 1];
	plots(pltcnt).vars  = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).zdeps = [0 50 150 300 450];
	plots(pltcnt).bal   = [2 2 2 2 2 2];
    plots(pltcnt).levs  = {linspace(0,350,40),linspace(0,350,40),linspace(0,350,40),linspace(0,350,40),linspace(0,350,40);...
                           linspace(0,40,40),linspace(0,40,40),linspace(0,40,40),linspace(0,40,40),linspace(0,40,40);...
                           linspace(0,3,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
                           linspace(0,1,40),linspace(0,1,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
                           linspace(0,0.02,40),linspace(0,0.05,40),linspace(0,0.08,40),linspace(0,0.08,40),linspace(0,0.08,40);...
						   linspace(0,3,40),linspace(0,3,40),linspace(0,1,40),linspace(0,1,40),linspace(0,1,40);...
						   linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40)};
    plots(pltcnt).dlevs = {linspace(-10,10,41),linspace(-5,5,41),linspace(-0.1,0.1,41),...
						  linspace(-1,1,41),linspace(-0.02,0.02,41),linspace(-0.5,0.5,41),linspace(-5,5,41)};

% (7) BGC: Latitude tracer slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1 1 1 1 1];
	plots(pltcnt).vars   = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps  = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).choice = 'lat';
	plots(pltcnt).lats   = 0;
	plots(pltcnt).xlims  = [140 260];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).bal    = [2 2 2 2 2 2];
	plots(pltcnt).levs   = {linspace(0,300,40);...
						    linspace(0,45,40);...
						    linspace(0,3,40);...
						    linspace(0,1,40);...
						    linspace(0,0.08,40);...
						    linspace(0,3,40);...
						    linspace(-25,25,41)};
	plots(pltcnt).dlevs  = {linspace(-10,10,41),linspace(-5,5,41),linspace(-0.1,0.1,41),...
						    linspace(-1,1,41),linspace(-0.02,0.02,41),linspace(-0.5,0.5,41),linspace(-5,5,41)};

% (8) BGC: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1];
	plots(pltcnt).bal   = [2];
	plots(pltcnt).vars  = {'NPP'};
	plots(pltcnt).cmaps = {'algae'};
	plots(pltcnt).levs  = {linspace(0,1200,40)};
	plots(pltcnt).dlevs = {linspace(-100,100,41)};

% (9) BGC: Longitude tracer slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1 1 1 1 1];
	plots(pltcnt).vars   = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps  = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).choice = 'lon';
	plots(pltcnt).lons   = [210 255 272];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).xlims  = [-30 50;-40 20;-40 10];
	plots(pltcnt).bal    = [2 2 2 2 2 2];
	plots(pltcnt).levs   = {linspace(0,300,40),linspace(0,300,40),linspace(0,300,40);...
						    linspace(0,45,40),linspace(0,45,40),linspace(0,45,40);...
						    linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						    linspace(0,5,40),linspace(0,5,40),linspace(0,5,40);...
						    linspace(0,0.08,40),linspace(0,0.08,40),linspace(0,0.08,40);...
						    linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						    linspace(-25,25,41),linspace(-25,25,41),linspace(-25,25,41)};
	plots(pltcnt).dlevs  = {linspace(-10,10,41),linspace(-5,5,41),linspace(-0.1,0.1,41),...
						    linspace(-1,1,41),linspace(-0.02,0.02,41),linspace(-0.5,0.5,41),linspace(-5,5,41)};
 
% (10) BGC: surface chlA
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).vars      = {'SFC_CHL'};
	plots(pltcnt).cmaps     = {'algae'};
	plots(pltcnt).absLevs   = log10([0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0 6.0 10.0]);
	plots(pltcnt).absLbls   = [0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0 6.0 10.0];
	plots(pltcnt).absCaxis  = real(log10([0.01 10]));
	plots(pltcnt).diffLevs  = [0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0];
	plots(pltcnt).diffLbls  = [0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0];
	plots(pltcnt).diffLevs  = [-fliplr(plots(pltcnt).diffLevs) plots(pltcnt).diffLevs];
	plots(pltcnt).diffLbls  = [-fliplr(plots(pltcnt).diffLbls) plots(pltcnt).diffLbls];
	plots(pltcnt).diffLevs  = romsMaster.dfloglevs(plots(pltcnt).diffLevs,0.001);
	plots(pltcnt).diffCaxis = [plots(pltcnt).diffLevs(1) plots(pltcnt).diffLevs(end)];

% (11) BGC: OMZ thickness
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).cmap      = cbrewer('seq','YlOrRd',40);
	plots(pltcnt).cmap(1,:) = [1 1 1];
	plots(pltcnt).omzthresh = [20 50];
	plots(pltcnt).levs      = {linspace(0,1000,40),linspace(0,1000,40)};
	plots(pltcnt).dlevs     = {linspace(-100,100,41),linspace(-100,100,41)};;


% (12) BGC: Integrated rates
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).vars    = {'AMMOX','NITROX','DENITRIF1','DENITRIF2','DENITRIF3','ANAMMOX','N2OAMMOX'};
	plots(pltcnt).opt     = [1 1 1 1 1 1 1];
	plots(pltcnt).convert = [1 1 1 1 2 2 2];
	plots(pltcnt).colr    = {rgb('DarkMagenta'),rgb('Goldenrod')}; 
	plots(pltcnt).cmaps   = {'tempo','tempo','tempo','tempo','tempo','tempo','tempo'};
	plots(pltcnt).levs    = {linspace(0,5e-11,64),linspace(0,5e-11,64),linspace(0,1e-10,64),linspace(0,1e-10,64),...
				             linspace(1,1e-10,64),linspace(0,1e-10,64),linspace(0,5e-14,64)};
	plots(pltcnt).dlevs   = {linspace(-5e-12,5e-12,65),linspace(-5e-12,5e-12,65),linspace(-5e-12,5e-12,65),linspace(-5e-12,5e-12,65),...
							 linspace(-5e-12,5e-12,65),linspace(-5e-12,5e-12,65),linspace(-5e-14,5e-14,65)};

% (13) BGC: Integrated tracers
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).vars    = {'O2','NO3','PO4','NO2','N2O','NH4'};
	plots(pltcnt).colr    = {rgb('DarkMagenta'),rgb('Goldenrod')}; 
	plots(pltcnt).cmaps   = {'-ice','tempo','tempo','tempo','tempo','tempo'};
	plots(pltcnt).levs    = {linspace(0,1e6,64),linspace(0,2e5,64),linspace(0,1.5e4,64),linspace(0,2000,64),linspace(0,150,64),linspace(0,500,64)};
	plots(pltcnt).dlevs   = {linspace(-1500,1500,65),linspace(-200,200,65),linspace(-10,10,65),linspace(-500,500,65),linspace(-5,5,65),linspace(-150,150,65)};

% (14) Other: depth slice locations
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).lats = plots(2).lats;
	plots(pltcnt).lons = plots(4).lons;

% (15) BGC: Latitude rate slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt     = [1 1 1 1 1 1 1];
	plots(pltcnt).vars    = {'AMMOX','NITROX','DENITRIF1','DENITRIF2','DENITRIF3','ANAMMOX','N2OAMMOX'};
	plots(pltcnt).cmaps   = {'tempo','tempo','tempo','tempo','tempo','tempo','tempo'}; 
    plots(pltcnt).convert = [1 1 1 1 2 2 2];
	plots(pltcnt).choice  = 'lat';
	plots(pltcnt).lats    = 0;
	plots(pltcnt).xlims   = [140 260];
	plots(pltcnt).zlims   = [0 1000];
	plots(pltcnt).bal     = [2 2 2 2 2 2];
	plots(pltcnt).levs    = {linspace(0,0.02,128);linspace(0,0.02,128);linspace(0,0.05,128);...  
						     linspace(0,0.05,128);linspace(0,0.05,128);linspace(0,0.05,128); linspace(0,5e-5,128)};    
	plots(pltcnt).dlevs   = {linspace(-0.02,0.02,128);linspace(-0.02,0.02,128);linspace(-0.05,0.05,128);...
					         linspace(-0.05,0.05,128);linspace(-0.05,0.05,128);linspace(-0.02,0.02,128);linspace(-5e-5,5e-5,128)};

% (16) BGC: Longitude rate slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt     = [1 1 1 1 1 1 1];
	plots(pltcnt).vars    = {'AMMOX','NITROX','DENITRIF1','DENITRIF2','DENITRIF3','ANAMMOX','N2OAMMOX'};
	plots(pltcnt).cmaps   = {'tempo','tempo','tempo','tempo','tempo','tempo','tempo'}; 
    plots(pltcnt).convert = [1 1 1 1 2 2 2];
	plots(pltcnt).choice  = 'lon';
    plots(pltcnt).lons    = [210 255 272];
    plots(pltcnt).zlims   = [0 1000];
    plots(pltcnt).xlims   = [-30 50;-40 20;-40 10];
    plots(pltcnt).bal     = [2 2 2 2 2 2];
	plots(pltcnt).levs    = {linspace(0,0.02,128);linspace(0,0.02,128);linspace(0,0.05,128);...  
						     linspace(0,0.05,128);linspace(0,0.05,128);linspace(0,0.05,128); linspace(0,5e-5,128)};    
	plots(pltcnt).dlevs   = {linspace(-0.02,0.02,128);linspace(-0.02,0.02,128);linspace(-0.05,0.05,128);...
					         linspace(-0.05,0.05,128);linspace(-0.05,0.05,128);linspace(-0.02,0.02,128);linspace(-5e-5,5e-5,128)};

% (17) BGC: Latitude rate slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt     = [1 1 1];
	plots(pltcnt).vars    = {'DENITRIF2','ANAMMOX','N2OAMMOX'};
	plots(pltcnt).cmaps   = {'-ice','-ice','-ice'}; 
    plots(pltcnt).convert = [1 2 2];
	plots(pltcnt).choice  = 'lat';
	plots(pltcnt).lats    = 0;
	plots(pltcnt).xlims   = [140 260];
	plots(pltcnt).zlims   = [0 1000];
	plots(pltcnt).bal     = [2 2 2];
	plots(pltcnt).levs    = {linspace(0,100,101);linspace(0,100,101);linspace(0,100,101);...  
						     linspace(0,100,101);linspace(0,100,101);linspace(0,100,101); linspace(0,100,101)};    
	plots(pltcnt).dlevs   = {linspace(-100,100,101);linspace(-100,100,101);linspace(-100,100,101);...  
					         linspace(-100,100,101);linspace(-100,100,101);linspace(-100,100,101); linspace(-100,100,101)};

% (18) BGC: Longitude rate slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt     = [1 1 1];
	plots(pltcnt).vars    = {'DENITRIF2','ANAMMOX','N2OAMMOX'};
	plots(pltcnt).cmaps   = {'-ice','-ice','-ice'}; 
    plots(pltcnt).convert = [1 2 2];
	plots(pltcnt).choice  = 'lon';
    plots(pltcnt).lons    = [210 255 272];
    plots(pltcnt).zlims   = [0 1000];
    plots(pltcnt).xlims   = [-30 50;-40 20;-40 10];
    plots(pltcnt).bal     = [2 2 2];
	plots(pltcnt).levs    = {linspace(0,100,101);linspace(0,100,101);linspace(0,100,101);...  
						     linspace(0,100,101);linspace(0,100,101);linspace(0,100,101); linspace(0,100,101)};    
	plots(pltcnt).dlevs   = {linspace(-100,100,101);linspace(-100,100,101);linspace(-100,100,101);...  
					         linspace(-100,100,101);linspace(-100,100,101);linspace(-100,100,101); linspace(-100,100,101)};

% (19) BGC: Fe depth slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
        plots(pltcnt).opt   = [1];
        plots(pltcnt).vars  = {'Fe'};
        plots(pltcnt).cmaps = {'thermal'};
        plots(pltcnt).zdeps = [0 50 150 300 450 600 800 1000 1500 2000];
        plots(pltcnt).bal   = [2];
        plots(pltcnt).levs  = {linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),...
                                                   linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),...
                                                   linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40)};
    plots(pltcnt).dlevs = {linspace(-1e-3,1e-3,41)};

% (20) POC_FLUX_IN comparisons
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
        plots(pltcnt).cmaps = {'deep'};
        plots(pltcnt).levs  = {linspace(0,5e-4,40)};
        plots(pltcnt).dlevs = {linspace(-5e-4,5e-4,41)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change to diagDir and load eoverrides options
compDir = ['/data/project1/demccoy/ROMS/',obj(1).info.simName,'/analysis/comp/'];
mkdir(compDir)
cd(compDir);
%try; run(['comp_overrides.m']);
%catch; disp('No overrides'); 
%end
comppath = [compDir,'plots/',obj(1).info.runName,'_VS_',obj(2).info.runName,'/'];
mkdir(comppath)

% Clear workspace and begin
clearvars -except obj plots comppath file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PHYSICAL DIAGNOSTICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D physical diagnostics (zslices)
% P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	% Inter-ROMS comparisons
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for o = 1:length(obj)
			obj(o) = clearROMS(obj(o));
			obj(o) = zslice(obj(o),vars(v),zdeps,file{o});
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		for z = 1:length(zdeps)
			close all
			roms1dat    = nanmean(squeeze(tmp{1}(:,:,z,:)),3);
			roms2dat    = nanmean(squeeze(tmp{2}(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms1']);
			close(figs(1));
			% ROMS2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms2']);
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms_diff']);
			close(figs(3));
		end
	end
	% Clear data
	for o = 1:length(obj)
		obj(o) = clearROMS(obj(o));
	end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical latitude section plots
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lats,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}; 
		roms2dat = tmp{2}; 
		for l = 1:length(lats)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface physical diagnostics
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	for o = 1:length(obj)
		obj(o) = clearROMS(obj(o));
		obj(o) = computeVar(obj(o),vars(opt==1),file{o});
	end
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		close all
		roms1dat = nanmean(obj(1).data.avg.(vars{v}).data,3);
		roms2dat = nanmean(obj(2).data.avg.(vars{v}).data,3);
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'bal',bal(v),'levels',levs{v},'difflevels',dlevs{v});
	   % ROMS figure
		set(0,'CurrentFigure',figs(1));
		title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms1']);
		close(figs(1));
		set(0,'CurrentFigure',figs(2));
		title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms2']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms_diff']);
		close(figs(3));
	end
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical longitude section plots
% P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lons,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}; 
		roms2dat = tmp{2}; 
		for l = 1:length(lons)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% latitude u-velocity slices 
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	for o = 1:length(obj)
		obj(o) = equatorUcmp(obj(o),'179E_160W');
	end
	roms1dat = obj(1).data.avg.u.slice;
	roms2dat = obj(2).data.avg.u.slice;
	[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'zlims',zlims,'xlims',xlims,...
		'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
	% ROMS figure
	set(0,'CurrentFigure',figs(1));
	title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{1}).name,': ',num2str(obj(1).slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(1),obj(1).data.avg.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[comppath,'equ_roms_roms1']);
	close(figs(1))
	% ROMS figure
	set(0,'CurrentFigure',figs(2));
	title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{1}).name,': ',num2str(obj(2).slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(2),obj(2).data.avg.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[comppath,'equ_roms_roms2']);
	close(figs(2))
	% Difference figure
	set(0,'CurrentFigure',figs(3));
	title(['ROMS Difference: ',num2str(obj(1).slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(3),obj(1).data.avg.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[compath,'equ_roms_diff']);
	close(figs(3))
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% BGC DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3d bgc diagnostics
% P6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	% Inter-ROMS comparisons
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for o = 1:length(obj)
			obj(o) = clearROMS(obj(o));
			obj(o) = zslice(obj(o),vars(v),zdeps,file{o});
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		for z = 1:length(zdeps)
			close all
			roms1dat    = nanmean(squeeze(tmp{1}(:,:,z,:)),3);
			roms2dat    = nanmean(squeeze(tmp{2}(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms1']);
			close(figs(1));
			% ROMS2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms2']);
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms_diff']);
			close(figs(3));
		end
	end
	% Clear data
	for o = 1:length(obj)
		obj(o) = clearROMS(obj(o));
	end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc latitude section plots
% P7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lats,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}; 
		roms2dat = tmp{2}; 
		for l = 1:length(lats)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface field plots
% P8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	for o = 1:length(obj)
		obj(o) = clearROMS(obj(o));
		obj(o) = loadData(obj(o),vars(opt==1),file{o});
	end
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		close all
		roms1dat = nanmean(obj(1).data.avg.(vars{v}).data,3);
		roms2dat = nanmean(obj(2).data.avg.(vars{v}).data,3);
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'bal',bal(v),'levels',levs{v},'difflevels',dlevs{v});
	   % ROMS figure
		set(0,'CurrentFigure',figs(1));
		title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms1']);
		close(figs(1));
		set(0,'CurrentFigure',figs(2));
		title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms2']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms_diff']);
		close(figs(3));
	end
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc longitude section plots
% P9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lons,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}; 
		roms2dat = tmp{2}; 
		for l = 1:length(lons)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface chla
% P10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	for o = 1:length(obj)
		obj(o) = loadData(obj(o),vars,file{o});
	end
	roms1dat = nanmean(obj(1).data.avg.(vars{1}).data,3);
	roms2dat = nanmean(obj(2).data.avg.(vars{1}).data,3);
	diffdat  = roms1dat - roms2dat;
	roms1dat = real(log10(roms1dat));
	roms2dat = real(log10(roms2dat));
	diffdat  = romsMaster.dfloglevs(diffdat,0.001);
	% Plot
	[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
	% ROMS1 figure
	set(0,'CurrentFigure',figs(1));
	title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{1}).name,': sfc'],'Interpreter','Latex');
	ylabel(cbs(1),obj(1).data.avg.(vars{1}).units,'Interpreter','Latex');
	cbs(1).XTickLabel = absLbls; 
	export_fig('-png',[comppath,vars{1},'_roms1']);
	% ROMS2 figure
	set(0,'CurrentFigure',figs(2));
	title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{1}).name,': sfc'],'Interpreter','Latex');
	ylabel(cbs(2),obj(2).data.avg.(vars{1}).units,'Interpreter','Latex');
	cbs(2).XTickLabel = absLbls; 
	export_fig('-png',[comppath,vars{1},'_roms2']);
	close(figs(2));
	% Diff figure
	set(0,'CurrentFigure',figs(3));
	title(['ROMS Difference'],'Interpreter','Latex');
	ylabel(cbs(3),obj(1).data.avg.(vars{1}).units,'Interpreter','Latex');
	cbs(3).XTick = diffLevs;
	cbs(3).XTickLabel = diffLbls;
	cbs(3).Limits = diffCaxis;
	export_fig('-png',[comppath,vars{1},'_roms_diff']);
	close(figs(3));
    
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
    clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMZ thickness
% P11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Get OMZ thickness
	for o = 1:length(obj)
		obj(o) = computeVar(obj(o),{'OMZ'},file{o},'thresh',omzthresh);
	end
    % Make comparison plots
	for th = 1:length(omzthresh)
		roms1dat = nanmean(squeeze(obj(1).data.avg.OMZ.int(:,:,th,:)),3);
		roms2dat = nanmean(squeeze(obj(2).data.avg.OMZ.int(:,:,th,:)),3);
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'levels',levs{th},'difflevels',dlevs{th});
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title([obj(1).info.runTit,' ',obj(1).data.avg.OMZ.name,'($O_2$ $<$ ',num2str(omzthresh(th)),' $mmol$ $m^{-3}$)'],...
			'Interpreter','Latex');
		ylabel(cbs(1),obj(1).data.avg.OMZ.units,'Interpreter','Latex')
		set(gcf,'ColorMap',cmap);
		export_fig('-png',[comppath,'OMZ_roms_th',num2str(omzthresh(th)),'_roms1']);
		close(figs(1));
		% ROMS figure
		set(0,'CurrentFigure',figs(2));
		title([obj(2).info.runTit,' ',obj(2).data.avg.OMZ.name,'($O_2$ $<$ ',num2str(omzthresh(th)),' $mmol$ $m^{-3}$)'],...
			'Interpreter','Latex');
		ylabel(cbs(2),obj(2).data.avg.OMZ.units,'Interpreter','Latex')
		set(gcf,'ColorMap',cmap);
		export_fig('-png',[comppath,'OMZ_roms_th',num2str(omzthresh(th)),'_roms2']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).data.avg.OMZ.units,'Interpreter','Latex')
		export_fig('-png',[comppath,'OMZ_roms_diff_th',num2str(omzthresh(th))]);
		close all
	end
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
    clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrated rates 
% P12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpactStruct(plots(pltcnt));
	% Convert mmol N/s to TgN/yr
	mmolNps_to_TgNpy = [(10^-3)*14*3600*24*365.25]/(1e12);
	% Go through each variable
	for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			if isempty(obj(o).grid.avg)
				obj(o) = loadDepth(obj(o),file{o});
			end
			obj(o) = loadData(obj(o),vars(v),file{o});
			obj(o) = intVar(obj(o),vars(v));
			out(o).int = obj(o).data.avg.(vars{v}).int;
			out(o).tot = obj(o).data.avg.(vars{v}).tot;
		end
		
		% Get total figure
        fig = piofigs('mfig',0.66);
        sb(1) = subplot(1,10,[1:7]);
        set(fig,'CurrentAxes',sb(1));
        for i = 1:length(colr)
            plot(NaN,NaN,'-','color',colr{i}); hold on
        end
        l = legend({obj(1).info.runTit,obj(2).info.runTit},'location','southoutside'); 
		l.Box = 'off'; l.AutoUpdate = 'off'; l.FontSize = 8;
        l.Interpreter = 'none'; l.Orientation = 'horizontal';
        all = [];
        for o = 1:length(obj);
            % Get data
            data = out(o).tot.*mmolNps_to_TgNpy.*convert(v);
            set(fig,'CurrentAxes',sb(1));
            plot(1:12,data(1,:),'-','color',colr{o},'linewidth',2);
            hold on
            all(o) = nanmean(sum(data,1));
        end
        xlim([1 12]);
        yl = get(gca,'YLim'); yl = [0 yl(end)]; ylim(yl);
        title(['Total ',lower(vars{v})],'FontSize',8);
        ylabel('$TgN$ $yr^{-1}$','Interpreter','Latex');
        %xlabel(['Model year ',num2str(obj(1).info.runYear)]);

        % Side panel
        sb(2) = subplot(1,10,[8:10]);
        b    = bar(1:length(all),all,'facecolor','flat');
        xlim([0.5 length(all)+0.5]);
        ylim([yl]);
        sb(2).Box = 'off'; sb(2).YTickLabel = []; sb(2).FontSize = 8;
        for o = 1:length(obj)
            b.CData(o,:) = colr{o};
        end
        b.EdgeColor = 'k'; b.LineWidth = 1;
		set(sb(2),'XTick',[1 2]);
        set(sb(2),'XTickLabel',[]);
        title(['Annual Mean'],'FontSize',8);
        xtickangle(45);
        sb(2).Position([2 4]) = [sb(2).Position(2)+0.1 sb(2).Position(4)-0.1];
        sb(1).Position([2 4]) = sb(2).Position([2 4]);
        % Save figure
        export_fig('-png',[comppath,'tot_',lower(vars{v})]);
        close all

		% Get integrated figure
		roms1dat = nanmean(out(1).int,3).*mmolNps_to_TgNpy.*convert(v);
		roms2dat = nanmean(out(2).int,3).*mmolNps_to_TgNpy.*convert(v);
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs{v},'difflevels',dlevs{v},'units','$TgN m^{-2} yr^{-1}$'); 
		% ROMS1 figure
		set(0,'CurrentFigure',figs(1));
		title({obj(1).info.runTit,['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms1']);
		close(figs(1));
		% ROMS2 figure
		set(0,'CurrentFigure',figs(2));
		title({obj(2).info.runTit,['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms2']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title({'ROMS Difference',['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms_diff']);
		close(figs(3));		
	end
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
    clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrated tracers 
% P13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpactStruct(plots(pltcnt));
	% Go through each variable
	for v = 1:length(vars)
		for o = 1:length(obj)
			if isempty(obj(o).grid.avg)
				obj(o) = loadDepth(obj(o),file{o});
			end
			obj(o) = loadData(obj(o),vars(v),file{o});
			obj(o) = intVar(obj(o),vars(v));
			out(o).int = obj(o).data.avg.(vars{v}).int;
			out(o).tot = obj(o).data.avg.(vars{v}).tot;
		end

		% Get total figure
        fig = piofigs('mfig',0.66);
        sb(1) = subplot(1,10,[1:7]);
        set(fig,'CurrentAxes',sb(1));
        for i = 1:length(colr)
            plot(NaN,NaN,'-','color',colr{i}); hold on
        end
        l = legend({obj(1).info.runTit,obj(2).info.runTit},'location','southoutside'); 
		l.Box = 'off'; l.AutoUpdate = 'off'; l.FontSize = 8;
        l.Interpreter = 'none'; l.Orientation = 'horizontal';
        all = [];
        for o = 1:length(obj);
            % Get data
            data = out(o).tot;
            set(fig,'CurrentAxes',sb(1));
            plot(1:12,data(1,:),'-','color',colr{o},'linewidth',2);
            hold on
            all(o) = nanmean(sum(data,1));
        end
        xlim([1 12]);
        yl = get(gca,'YLim'); yl = [0 yl(end)]; ylim(yl);
        title(['Total ',lower(vars{v})],'FontSize',8);
        ylabel(obj(1).data.avg.(vars{v}).totunits,'Interpreter','Latex');
        %xlabel(['Model year ',num2str(obj(1).info.runYear)]);

        % Side panel
        sb(2) = subplot(1,10,[8:10]);
        b    = bar(1:length(all),all,'facecolor','flat');
        xlim([0.5 length(all)+0.5]);
        ylim([yl]);
        sb(2).Box = 'off'; sb(2).YTickLabel = []; sb(2).FontSize = 8;
        for o = 1:length(obj)
            b.CData(o,:) = colr{o};
        end
        b.EdgeColor = 'k'; b.LineWidth = 1;
		set(sb(2),'XTick',[1 2]);
        set(sb(2),'XTickLabel',[]);
        title(['Annual Mean'],'FontSize',8);
        xtickangle(45);
        sb(2).Position([2 4]) = [sb(2).Position(2)+0.1 sb(2).Position(4)-0.1];
        sb(1).Position([2 4]) = sb(2).Position([2 4]);
        % Save figure
        export_fig('-png',[comppath,'tot_',lower(vars{v})]);
        close all

		% Get integrated figure
		roms1dat = nanmean(out(1).int,3);
		roms2dat = nanmean(out(2).int,3);
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs{v},...
			'difflevels',dlevs{v},'units',obj(1).data.avg.(vars{v}).intunits); 
		% ROMS1 figure
		set(0,'CurrentFigure',figs(1));
		title({obj(1).info.runTit,['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms1']);
		close(figs(1));
		% ROMS2 figure
		set(0,'CurrentFigure',figs(2));
		title({obj(2).info.runTit,['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms2']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title({'ROMS Difference',['Integrated ',lower(vars{v})]},'Interpreter','none');
		export_fig('-png',[comppath,'int_',lower(vars{v}),'_roms_diff']);
		close(figs(3));		
	end
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
    clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps
% P14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
	% Get lon lines
	for i = 1:length(lons)
		lonidx  = lons(i);
		lonx{i} = obj(1).grid.lon_rho(lonidx,:);
		lony{i} = obj(1).grid.lat_rho(lonidx,:);
	end
	% Get lat lines
	for i = 1:length(lats)
		latidx  = lats(i);
		latx{i} = obj(1).grid.lon_rho(:,latidx);
		laty{i} = obj(1).grid.lat_rho(:,latidx);
	end
	% Plot
	[fig] = quickMap(obj(1),'ticks',1,'fontsize',8);
	hold on
	m_plot(obj(1).grid.polygon(:,1),obj(1).grid.polygon(:,2),'k','linewidth',2);
	for i = 1:length(lons)
		m_plot(lonx{i},lony{i},'--k');
	end
	for i = 1:length(lats)
		m_plot(latx{i},laty{i},'--k');
	end
	title(['Location of depth slices'],'Interpreter','Latex');
	export_fig('-png',[comppath,'trans_locations']);
    % Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
    clearvars -except obj plots comppath file pltcnt	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latitude slice rates 
% P15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lats,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}.*86400.*convert(v); 
		roms2dat = tmp{2}.*86400.*convert(v); 
		for l = 1:length(lats)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitude slice rates 
% P16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for o = 1:length(obj)
			obj(o) = sliceROMS(obj(o),vars(v),choice,lons,file{o},'type','z_avg','zlim',zlims(end));
			tmp{o} = obj(o).data.avg.(vars{v}).slice;
		end
		roms1dat = tmp{1}.*86400.*convert(v); 
		roms2dat = tmp{2}.*86400.*convert(v); 
		for l = 1:length(lons)
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),'$mmol$ $N$ $m^{-3}$ $d^{-1}$','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latitude slice N-loss 
% P17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
	for o = 1:length(obj)
		% Load vars, calculate total N-loss
		obj(o) = sliceROMS(obj(o),vars,choice,lats,file{o},'type','z_avg','zlim',zlims(end));
		obj(o).data.avg.NLoss.slice = obj(o).data.avg.N2OAMMOX.slice.*2 + ...
									  obj(o).data.avg.DENITRIF2.slice.*1 + ...
									  obj(o).data.avg.ANAMMOX.slice.*2; 
		% Blank low rates to clean up plots
		obj(o).data.avg.NLoss.slice((obj(o).data.avg.NLoss.slice.*86400<0.001)) = NaN;
		obj(o).data.avg.NLoss.name  = 'averaged total N-loss'; 
		obj(o).data.avg.NLoss.units = '$mmol$ $N$ $m^{-3}$ $s^-1$';
	end
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		roms1dat  = obj(1).data.avg.(vars{v}).slice.*86400.*convert(v); 
		roms2dat  = obj(2).data.avg.(vars{v}).slice.*86400.*convert(v); 
		roms1frac = (roms1dat./(obj(1).data.avg.NLoss.slice.*86400)).*100;
		roms2frac = (roms2dat./(obj(2).data.avg.NLoss.slice.*86400)).*100;
		for l = 1:length(lats)
			[figs,cbs] = sliceCmp(obj(1),roms1frac,roms2frac,0,'cmap',cmaps{v},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,' fraction: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lat',num2str(lats(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,' fraction: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lat',num2str(lats(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lat',num2str(lats(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitude slice N-loss 
% P18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpactStruct(plots(pltcnt));
    % Inter-ROMS comparisons
	for o = 1:length(obj)
		% Load vars, calculate total N-loss
		obj(o) = sliceROMS(obj(o),vars,choice,lons,file{o},'type','z_avg','zlim',zlims(end));
		obj(o).data.avg.NLoss.slice = obj(o).data.avg.N2OAMMOX.slice.*2 + ...
									  obj(o).data.avg.DENITRIF2.slice.*1 + ...
									  obj(o).data.avg.ANAMMOX.slice.*2; 
		% Blank low rates to clean up plots
		obj(o).data.avg.NLoss.slice((obj(o).data.avg.NLoss.slice.*86400<0.001)) = NaN;
		obj(o).data.avg.NLoss.name  = 'averaged total N-loss'; 
		obj(o).data.avg.NLoss.units = '$mmol$ $N$ $m^{-3}$ $s^-1$';
	end
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		roms1dat  = obj(1).data.avg.(vars{v}).slice.*86400.*convert(v); 
		roms2dat  = obj(2).data.avg.(vars{v}).slice.*86400.*convert(v); 
		roms1frac = (roms1dat./(obj(1).data.avg.NLoss.slice.*86400)).*100;
		roms2frac = (roms2dat./(obj(2).data.avg.NLoss.slice.*86400)).*100;		
		for l = 1:length(lons)
			[figs,cbs] = sliceCmp(obj(1),roms1frac,roms2frac,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,' fraction: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lon',num2str(lons(l)),'_roms1']);
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,' fraction: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lon',num2str(lons(l)),'_roms2']);
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),'Percent','Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_frac_lon',num2str(lons(l)),'_roms_diff']);
			close(figs(3))
		end
	end
	% Clear data
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
	clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe depth slices 
% P19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
        unpactStruct(plots(pltcnt));
        % Go through each compare with diagnostics
        % Variable loop
        for v = 1:length(vars)
                if ~opt(v)
                        disp(['...skipping ',vars{v},'...']);
                        continue
                end
                for o = 1:length(obj)
                        obj(o) = clearROMS(obj(o));
                        obj(o) = zslice(obj(o),vars(v),zdeps,file{o});
                        tmp{o} = obj(o).data.avg.(vars{v}).slice;
                end
                % Depth loop
                for z = 1:length(zdeps)
                        close all
                        roms1dat    = nanmean(squeeze(tmp{1}(:,:,z,:)),3);
                        roms2dat    = nanmean(squeeze(tmp{2}(:,:,z,:)),3);
                        [figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
                        % ROMS1 figure
                        set(0,'CurrentFigure',figs(1));
                        title([obj(1).info.runTit,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                        ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
                        export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms1'],'-m2.5');
                        close(figs(1));
                        % ROMS2 figure
                        set(0,'CurrentFigure',figs(2));
                        title([obj(2).info.runTit,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                        ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
                        export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms2'],'-m2.5');
                        close(figs(2));
                        % Diff figure
                        set(0,'CurrentFigure',figs(3));
                        title(['ROMS Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                        ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
                        export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms_diff'],'-m2.5');
                        close(figs(3));
                end
        end
    for o = 1:length(obj)
        obj(o) = clearROMS(obj(o));
    end
        clearvars -except obj plots comppath file pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps of POC_FLUX_IN
% P20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
        unpactStruct(plots(pltcnt));
      
        % Load POC FLUX IN
        for o = 1:length(obj)
       	  obj(o) = clearROMS(obj(o));
       	  obj(o) = zslice(obj(o),{'POC_FLUX_IN'},75,file{o});
        end
                close all
                roms1dat    = nanmean(obj(1).data.avg.POC_FLUX_IN.slice,3);
                roms2dat    = nanmean(obj(2).data.avg.POC_FLUX_IN.slice,3);
                [figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{1},'levels',levs{1},'difflevels',dlevs{1});
                % ROMS1 figure
                set(0,'CurrentFigure',figs(1));
                title([obj(1).info.runTit,' ',obj(1).data.avg.POC_FLUX_IN.name,': 75m'],'Interpreter','Latex');
                ylabel(cbs(1),obj(1).data.avg.POC_FLUX_IN(1).units,'Interpreter','Latex');
                export_fig('-png',[comppath,'POC_FLUX_IN_roms1'],'-m2.5');
                close(figs(1));
                % ROMS2 figure
                set(0,'CurrentFigure',figs(2));
                title([obj(2).info.runTit,' ',obj(2).data.avg.POC_FLUX_IN.name,': 75m'],'Interpreter','Latex');
                ylabel(cbs(2),obj(2).data.avg.POC_FLUX_IN(1).units,'Interpreter','Latex');
                export_fig('-png',[comppath,'POC_FLUX_IN_roms2'],'-m2.5');
                close(figs(2));
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['ROMS Difference'],'Interpreter','Latex');
                ylabel(cbs(3),obj(1).data.avg.POC_FLUX_IN(1).units,'Interpreter','Latex');
                export_fig('-png',[comppath,'POC_FLUX_IN_diff'],'-m2.5');



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
