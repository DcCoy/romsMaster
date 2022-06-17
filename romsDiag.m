function obj = romsDiag(obj,plotchoice,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to automatically generate diagnostic plots for individual ROMS simulations 
% or compare two simulations with one-another.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% - obj = romsDiag(obj,plotchoice)
%
% Inputs:
% - obj = roms object(s) that are already initialized
%		... if length(obj) == 1, then validation will occur
%		... if length(obj) == 2, then intercomparison will occur
%			... comparisons will be obj - obj(2)
% - plotchoice (see below) 
%	------------------------------');
%	---- PHYSICAL DIAGNOSTICS ----');
%	------------------------------');
%	1  == 3d physical diagnostics');
%	2  == equator physical slices');
%	3  == surface physical fields');
%	4  == longitude physical slices');
%	5  == u velocity at equator');
%	-------------------------');
%	---- BGC DIAGNOSTICS ----');
%	-------------------------');
%	6  == 3d bgc diagnostics');
%	7  == equator bgc slices');
%	8  == surface bgc fields');
%	9  == longitude bgc slices');
%	10 == surface chlA');
%	11 == OMZ thickness');
%	12 == N-cycle profile comparisons');
%	13 == ROMS vs OCIM (N2O)
%	------------------------------');
%	------ OTHER DIAGNOSTICS -----');
%	------------------------------');
%	14 == slice degree maps');
%
% Optional Inputs (varargin);
% - comp = Only perform inter-run comparisons (1)
% - runNames = cell array of runNames to compare
% - runYears = cell array of runYears to compare (NOT WORKING)
%
% Example:
% - obj = romsDiag(obj,[1:11]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle
addpath /data/project1/demccoy/ROMS/validation/n2o

% Process optional inputs
A.runNames = {obj.info.runName};
A.runYears = obj.info.runYear;
A.comp  = [0];
A.opt   = [];
A.vars  = [];
A.cmaps = [];
A.zdeps = [];
A.levs  = [];
A.dlevs = [];
A.lats  = [];
A.lons  = [];
A.xlims = [];
A.ylims = [];
A.zlims = [];
A.bal   = [];
A = parse_pv_pairs(A,varargin);

% Clear objects
obj = clearROMS(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHYSICAL PLOT options
% (1) PHYS: 3D zslices
pltcnt = 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1];
	plots(pltcnt).vars  = {'temp','salt','sigma'};
	plots(pltcnt).cmaps = {'thermal','haline','dense'};
	plots(pltcnt).zdeps = [0 50 150 300];
	plots(pltcnt).levs  = {linspace(3,30,40),linspace(3,30,40),linspace(0,25,40),linspace(2,18,40);...
						   linspace(30,37,40),linspace(30,37,40),linspace(33.4,36.5,40),linspace(33.5,35.5,40);...
						   linspace(19,27,40),linspace(19,27,40),linspace(24,28,40),linspace(26.5,28.6,40)};
	plots(pltcnt).dlevs = {linspace(-3,3,40),linspace(-0.6,0.6,40),linspace(-0.6,0.6,40)};

% (2) PHYS: Equator slices
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1];
	plots(pltcnt).vars  = {'temp','salt','sigma'};
	plots(pltcnt).cmaps = {'thermal','haline','dense'};
	plots(pltcnt).lats  = 0;
	plots(pltcnt).xlims = [140 260];
	plots(pltcnt).zlims = [0 1000];
	plots(pltcnt).levs  = {linspace(3,30,40),linspace(34.5,35.5,40),linspace(22,32,40)};
	plots(pltcnt).dlevs = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (3) PHYS: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1];
	plots(pltcnt).vars = {'MLD','SSH'};
	plots(pltcnt).cmaps = {'deep','deep'};
	plots(pltcnt).bal   = [2 2];
	plots(pltcnt).levs  = {linspace(0,150,40),linspace(0,1.5,40)};
	plots(pltcnt).dlevs = {linspace(-60,60,40),linspace(-0.3,0.3,40)};

% (4) PHYS: Longitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1];
	plots(pltcnt).vars  = {'temp','salt','sigma'};
	plots(pltcnt).cmaps = {'thermal','haline','dense'};
	plots(pltcnt).lons  = [210 255 272];
	plots(pltcnt).zlims = [0 1000];
	plots(pltcnt).xlims = [-30 50;-40 20;-40 10];
	plots(pltcnt).levs  = {linspace(3,30,40),linspace(3,30,40),linspace(3,30,40);...
						   linspace(33,36.5,40),linspace(33.5,36.5,40),linspace(33.5,36.5,40);...
						   linspace(22,32,40),linspace(22,32,40),linspace(22,32,40)};
	plots(pltcnt).dlevs = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (5) PHYS: Zonal equator velocity
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).zlims = [50 500];
	plots(pltcnt).xlims = [-20 20];
	plots(pltcnt).vars  = {'u'};
	plots(pltcnt).cmaps = {'balance'};
	plots(pltcnt).levs  = {linspace(-0.3,0.3,40)};
	plots(pltcnt).dlevs = {linspace(-0.3,0.3,40)};

% BIOGEOCHEMICAL PLOT options
% (6) BGC: 3D zslices
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1 1 1 1 1];
	plots(pltcnt).vars  = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).zdeps = [0 50 150 300 450];
	plots(pltcnt).bal   = [2 2 2 2 2 2];
    plots(pltcnt).levs  = {linspace(0,350,21),linspace(0,350,21),linspace(0,350,21),linspace(0,350,21),linspace(0,350,21);...
                           linspace(0,40,21),linspace(0,40,21),linspace(0,40,21),linspace(0,40,21),linspace(0,40,21);...
                           linspace(0,3,21),linspace(0,3,21),linspace(0,3,21),linspace(0,3,21),linspace(0,3,21);...
                           linspace(0,1,21),linspace(0,1,21),linspace(0,3,21),linspace(0,3,21),linspace(0,3,21);...
                           linspace(0,0.02,21),linspace(0,0.05,21),linspace(0,0.08,21),linspace(0,0.08,21),linspace(0,0.08,21);...
						   linspace(0,3,21),linspace(0,3,21),linspace(0,1,21),linspace(0,1,21),linspace(0,1,21);...
						   linspace(-25,25,21),linspace(-25,25,21),linspace(-25,25,21),linspace(-25,25,21),linspace(-25,25,21)};
    plots(pltcnt).dlevs = {linspace(-80,80,21),linspace(-10,10,21),linspace(-0.6,0.6,21),...
                           linspace(-0.6,0.6,21),linspace(-0.05,0.05,21),linspace(-0.6,0.6,21),linspace(-10,10,21)};

% (7) BGC: Equator slices
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1 1 1 1];
	plots(pltcnt).vars  = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).lats  = 0;
	plots(pltcnt).xlims = [140 260];
	plots(pltcnt).zlims = [0 1000];
	plots(pltcnt).bal   = [2 2 2 2 2 2];
	plots(pltcnt).levs  = {linspace(0,300,40);...
						   linspace(0,45,40);...
						   linspace(0,3,40);...
						   linspace(0,1,40);...
						   linspace(0,0.08,40);...
						   linspace(0,3,40);...
						   linspace(-25,25,41)};
	plots(pltcnt).dlevs = {linspace(-80,80,41),linspace(-15,15,41),linspace(-0.5,0.5,41),...
						   linspace(-1,1,41),linspace(-0.05,0.05,41),linspace(-1,1,41),linspace(-10,10,41)};

	% Overrides
	plots(pltcnt).opt   = [0 0 0 0 0 0 1];

% (8) BGC: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).vars  = {'NPP'};
	plots(pltcnt).cmaps = {'algae'};
	plots(pltcnt).levs  = {linspace(0,1000,40)};
	plots(pltcnt).dlevs = {linspace(-700,700,40)};

% (9) BGC: Longitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1 1 1 1];
	plots(pltcnt).vars  = {'O2','NO3','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).lons  = [210 255 272];
	plots(pltcnt).zlims = [0 1000];
	plots(pltcnt).xlims = [-30 50;-40 20;-40 10];
	plots(pltcnt).bal   = [2 2 2 2 2 2];
	plots(pltcnt).levs  = {linspace(0,300,40),linspace(0,300,40),linspace(0,300,40);...
						   linspace(0,45,40),linspace(0,45,40),linspace(0,45,40);...
						   linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						   linspace(0,0.4,40),linspace(0,0.4,40),linspace(0,0.4,40);...
						   linspace(0,0.08,40),linspace(0,0.08,40),linspace(0,0.08,40);...
						   linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						   linspace(-25,25,41),linspace(-25,25,41),linspace(-25,25,41)};
	plots(pltcnt).dlevs = {linspace(-80,80,41),linspace(-15,15,41),linspace(-0.5,0.5,41),...
						   linspace(-1,1,41),linspace(-0.05,0.05,41),linspace(-1,1,41),linspace(-10,10,41)};

% (10) BGC: surface chlA
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
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
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).cmap      = cbrewer('seq','YlOrRd',21);
	plots(pltcnt).cmap(1,:) = [1 1 1];
	plots(pltcnt).omzthresh = [20 50];
	plots(pltcnt).levs      = {linspace(0,1000,21),linspace(0,2000,21)};
	plots(pltcnt).dlevs     = {linspace(-1000,1000,41),linspace(-2000,2000,41)};;


% (12) BGC: Ncycle profiles
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).region_name  = {     'Eq',            'ETSP',              'ETNP',         'CCS',           'NPSG'   }; 
	plots(pltcnt).obs_regions  = {[217 222 -2 3],[278.0 283.0 -17 -12],  [250 255 14 19],[237 242 30 35],[206 211 28 33]};
	plots(pltcnt).roms_regions = {[217 222 -2 3],[273.5 278.5 -9.5 -4.5],[255 260 13 18],[237 242 30 35],[206 211 28 33]};
	plots(pltcnt).tracer_lims  = {    [0 50],            [0 50],              [0 50],		[0 50],          [0 50];...    % NO3
									  [0 1.5],		     [0 7],				  [0 6],		[0 0.5],         [0 0.5];...   % NO2
									  [0 50],		     [0 160],			  [0 70],		[0 60],          [0 60];...    % N2O*100
									  [0 1.5],		     [0 3.5],			  [0 1.5],		[0 2.5],         [0 0.5];...   % NH4
									  [0 250],		     [0 250],			  [0 220],		[0 300],         [0 300]};     % O2
	plots(pltcnt).ratio_lims   = {	  [0 1.5],			 [0 7],				  [0 6],        [0 0.5],         [0 0.5];...   % NO2
									  [0 1.5],			 [0 3.5],			  [0 1.5],      [0 2.5],         [0 0.5];...   % NH4
									  [1E-3 1E3] ,		 [1E-3 1E3],		  [1E-3 1E3],   [1E-3 1E3],      [1E-3 1E3];...% NO2vNH4
									  [1E-3 1E3] ,		 [1E-3 1E3],		  [1E-3 1E3],   [1E-3 1E3],      [1E-3 1E3];...% NH4vNO2
									  [0 210],			 [0 250],			  [0 220],      [0 300],         [0 300]};	   % O2

% (13) BGC: ROMS vs OCIM
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).regions   = {[225 300 -30 30]};

% (14) OTHER: Slice maps
pltcnt = pltcnt + 1;
plots(pltcnt).choice = plotchoice(pltcnt);
	plots(pltcnt).lons = [210 255 272];
	plots(pltcnt).lats = [0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get plotchoice
tmpchoice = zeros(1,100);
tmpchoice(plotchoice) = 1;
plotchoice = tmpchoice;

% Load overrides options
run(['/data/project1/demccoy/ROMS/',obj.info.simName,'/analysis/diag/diag_overrides.m']);

% Reset pltcnt
pltcnt = 0;

% Clear workspace and begin
clearvars -except obj A varargin plots pltcnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PHYSICAL DIAGNOSTICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D physical diagnostics (zslices)
% P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Diagnostic plots
	if A.comp == 0
		% runName loop
		for i = 1:length(runNames)	
			obj = changeInputs(obj,'runName',runNames{i});
			% Variable loop
			for j = 1:length(vars)
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = clearROMS(obj);
				obj = loadData(obj,vars(j),'type','z_avg','depth',zdeps);
				obj = loadDiag(obj,vars(j),zdeps);
				% Depth loop
				for k = 1:length(zdeps)
					close all
					romsdat    = nanmean(squeeze(obj.romsData.(vars{j})(1).data(:,:,k,:)),3);
					diagdat    = nanmean(squeeze(obj.diagData.(vars{j})(1).data(:,:,k,:)),3);
					[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{j},'levels',levs{j,k},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_roms']);
					close(figs(1));
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_diag']);
					close(figs(2));	
					% Diff figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_diff']);	
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end	
	% Inter-ROMS comparisons
	elseif A.comp == 1;	
		for i = 1:length(vars)
			obj = clearROMS(obj);
			if ~opt(i)
				disp(['...skipping ',vars{i},'...']);
				continue
			end
			for j = 1:length(runNames)
				obj = changeInputs(obj,'runName',runNames{j});
				obj = loadData(obj,vars{i},'type','z_avg','depth',zdeps);
				tmp{j} = obj.romsData.(vars{i}).data;
				if j == 1
					comppath = obj.paths.plots.comp;
				end
			end
			for j = 1:length(zdeps)
				close all
				roms1dat    = nanmean(squeeze(tmp{1}(:,:,j,:)),3);
				roms2dat    = nanmean(squeeze(tmp{2}(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'cmap',cmaps{i},'levels',levs{i,k},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_z',num2str(zdeps(j)),'_roms_diff']);
				close(figs(3));
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical equatorial section plots
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Diagnostics
	if A.comp == 0
        % runName loop
        for i = 1:length(runNames)
            obj = changeInputs(obj,'runName',runNames{i});
			for j = 1:length(vars);
				obj = clearROMS(obj);
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = sliceROMS(obj,vars(j),'lat',lats);
				obj = sliceDiag(obj,vars(j),'lat',lats);
				for k = 1:length(lats)
					romsdat = obj.romsData.(vars{j}).slice;
					diagdat = obj.diagData.(vars{j}).slice;
					[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{j},'xlims',xlims(k,:),'zlims',zlims,...
						'figdim',0.5,'slice',k,'levels',levs{j,k},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_roms']);
					close(figs(1))
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_diag']);
					close(figs(2))
					% Difference figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_diff']);
					close(figs(3))
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end
	% Inter-run comparisons
	elseif A.comp == 1
        for i = 1:length(vars)
			obj = clearROMS(obj);
            if ~opt(i)
                disp(['...skipping ',vars{i},'...']);
                continue
            end
            for j = 1:length(runNames)
                obj = changeInputs(obj,'runName',runNames{j});
				obj = sliceROMS(obj,vars(i),'lat',lats);
                tmp{j} = obj.romsData.(vars{i}).slice;
                if j == 1
                    comppath = obj.paths.plots.comp;
                end
            end
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
			for j = 1:length(lats)
				[figs,cbs] = sliceCmp(obj,roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
					'figdim',0.5,'slice',j,'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lats(j)),'$^oN$'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_lat',num2str(lats(j)),'_roms_diff']);
				close(figs(3))
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface physical diagnostics
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	for i = 1:length(obj)
		obj = loadData(obj,vars(opt==1),'type','raw');
	end
	if A.comp == 0
		obj = loadDiag(obj,vars(opt==1),0);
	end
	for i = 1:length(vars)
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		obj = clearROMS(obj);
		for j = 1:length(obj.diagData.(vars{i}))
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(obj.romsData.(vars{i})(1).data,3);
				diagdat    = nanmean(obj.diagData.(vars{i})(j).data,3);
				[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj.romsData.(vars{i})(1).name],'Interpreter','Latex');
				ylabel(cbs(1),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj.diagData.(vars{i})(j).name],'Interpreter','Latex');
				ylabel(cbs(2),obj.diagData.(vars{i})(j).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_diag_',num2str(j)]);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_diff_',num2str(j)]);
				close(figs(3));
			end
		end
		if length(obj)==2
			close all
			roms1dat    = nanmean(obj.romsData.(vars{i})(1).data,3);
			roms2dat    = nanmean(obj(2).romsData.(vars{i})(1).data,3);
			[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1));
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
			export_fig('-pdf',[obj.paths.plots.comp,vars{i},'_roms_diff']);
			close(figs(3));
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical longitude section plots
% P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Diagnostics
	if A.comp == 0
        % runName loop
        for i = 1:length(runNames)
            obj = changeInputs(obj,'runName',runNames{i});
			for j = 1:length(vars);
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = clearROMS(obj);
				obj = sliceROMS(obj,vars(j),'lon',lons);
				obj = sliceDiag(obj,vars(j),'lon',lons);
				for k = 1:length(lons)
					romsdat = obj.romsData.(vars{j}).slice;
					diagdat = obj.diagData.(vars{j}).slice;
					[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{j},'xlims',xlims(k,:),'zlims',zlims,...
						'figdim',0.5,'slice',k,'levels',levs{j,k},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_roms']);
					close(figs(1))
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_diag']);
					close(figs(2))
					% Difference figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_diff']);
					close(figs(3))
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end
	% Inter-run comparisons
	elseif A.comp == 1
        for i = 1:length(vars)
            if ~opt(i)
                disp(['...skipping ',vars{i},'...']);
                continue
            end
			obj = clearROMS(obj);
            for j = 1:length(runNames)
                obj = changeInputs(obj,'runName',runNames{j});
				obj = sliceROMS(obj,vars(i),'lon',lons);
                tmp{j} = obj.romsData.(vars{i}).slice;
                if j == 1
                    comppath = obj.paths.plots.comp;
                end
            end
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
			for j = 1:length(lons)
				[figs,cbs] = sliceCmp(obj,roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims(j,:),'zlims',zlims,...
					'figdim',0.5,'slice',j,'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_lon',num2str(lons(j)),'_roms_diff']);
				close(figs(3))
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equator u-velocity slices 
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	for i = 1:length(obj)
		obj = equatorUcmp(obj,'179E_160W');
	end
	for i = 1:length(obj)
		% Skip validations?
		if A.comp > 0
			continue
		end
		diagdat = obj.diagData.u.slice;
		romsdat = nanmean(obj.romsData.u.slice,3);
		[figs,cbs] = sliceCmp(obj,romsdat,diagdat,1,'zlims',zlims,'xlims',xlims,...
			'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj.romsData.(vars{1}).name,': ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(1),obj.romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj.paths.plots.diag,'equ_roms']);
		close(figs(1))
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title([obj.diagData.(vars{1}).name,': ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(2),obj.diagData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj.paths.plots.diag,'equ_diag']);
		close(figs(2))
		% Difference figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference: ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(3),obj.romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj.paths.plots.diag,'equ_diff']);
		close(figs(3))
	end
	if length(obj)==2
		roms1dat = obj.romsData.u.slice;
		roms2dat = obj(2).romsData.u.slice;
		[figs,cbs] = sliceCmp(obj,roms1dat,roms2dat,0,'zlims',zlims,'xlims',xlims,...
			'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
		close(figs(1));
		close(figs(2));
		% Difference figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference: ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(3),obj.romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj.paths.plots.comp,'equ_roms_diff']);
		close(figs(3))
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% BGC DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3d bgc diagnostics
% P6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Go through each compare with diagnostics
	if A.comp == 0
		% runName loop
		for i = 1:length(runNames)	
			obj = changeInputs(obj,'runName',runNames{i});
			% Variable loop
			for j = 1:length(vars)
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = clearROMS(obj);
				obj = loadData(obj,vars(j),'type','z_avg','depth',zdeps);
				obj = loadDiag(obj,vars(j),zdeps);
				% Depth loop
				for k = 1:length(zdeps)
					close all
					romsdat    = nanmean(squeeze(obj.romsData.(vars{j})(1).data(:,:,k,:)),3);
					diagdat    = nanmean(squeeze(obj.diagData.(vars{j})(1).data(:,:,k,:)),3);
					[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{j},'levels',levs{j,k},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_roms']);
					close(figs(1));
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_diag']);
					close(figs(2));	
					% Diff figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(zdeps(k)),'m'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_z',num2str(zdeps(k)),'_diff']);	
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end	
	elseif A.comp == 1;	
		for i = 1:length(vars)
			if ~opt(i)
				disp(['...skipping ',vars{i},'...']);
				continue
			end
			obj = clearROMS(obj);
			for j = 1:length(runNames)
				obj = changeInputs(obj,'runName',runNames{j});
				obj = loadData(obj,vars(i),'type','z_avg','depth',zdeps);
				tmp{j} = obj.romsData.(vars{i}).data;
				if j == 1
					comppath = obj.paths.plots.comp;
				end
			end
			for j = 1:length(zdeps)
				close all
				roms1dat    = nanmean(squeeze(tmp{1}(:,:,j,:)),3);
				roms2dat    = nanmean(squeeze(tmp{2}(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'cmap',cmaps{i},'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_z',num2str(zdeps(j)),'_roms_diff']);
				close(figs(3));
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc equatorial section plots
% P7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Diagnostics
	if A.comp == 0
        % runName loop
        for i = 1:length(runNames)
            obj = changeInputs(obj,'runName',runNames{i});
			for j = 1:length(vars);
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = clearROMS(obj);
				obj = sliceROMS(obj,vars(j),'lat',lats);
				obj = sliceDiag(obj,vars(j),'lat',lats);
				for k = 1:length(lats)
					romsdat = obj.romsData.(vars{j}).slice;
					diagdat = obj.diagData.(vars{j}).slice;
					[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{j},'xlims',xlims(k,:),'zlims',zlims,...
						'figdim',0.5,'slice',k,'levels',levs{j},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_roms']);
					close(figs(1))
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_diag']);
					close(figs(2))
					% Difference figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(lats(k)),'$^oN$'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lat',num2str(lats(k)),'_diff']);
					close(figs(3))
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end
	% Inter-run comparisons
	elseif A.comp == 1
        for i = 1:length(vars)
            if ~opt(i)
                disp(['...skipping ',vars{i},'...']);
                continue
            end
			obj = clearROMS(obj);
            for j = 1:length(runNames)
                obj = changeInputs(obj,'runName',runNames{j});
				obj = sliceROMS(obj,vars(i),'lat',lats);
                tmp{j} = obj.romsData.(vars{i}).slice;
                if j == 1
                    comppath = obj.paths.plots.comp;
                end
            end
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
			for j = 1:length(lats)
				[figs,cbs] = sliceCmp(obj,roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
					'figdim',0.5,'slice',j,'levels',levs{i},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lats(j)),'$^oN$'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_lat',num2str(lats(j)),'_roms_diff']);
				close(figs(3))
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface field plots
% P8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	for i = 1:length(obj)
		obj = loadData(obj,vars,'type','raw');
	end
	if A.comp == 0
		obj = loadDiag(obj,vars,0);
	end
	for i = 1:length(vars)
		obj = clearROMS(obj);
		for j = 1:length(obj.diagData.(vars{i}))
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(obj.romsData.(vars{i})(1).data,3);
				diagdat    = nanmean(obj.diagData.(vars{i})(j).data,3);
				[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{i},'levels',levs{i},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj.romsData.(vars{i})(1).name],'Interpreter','Latex');
				ylabel(cbs(1),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj.diagData.(vars{i})(j).name],'Interpreter','Latex');
				ylabel(cbs(2),obj.diagData.(vars{i})(j).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_diag_',num2str(j)]);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj.paths.plots.diag,vars{i},'_diff_',num2str(j)]);
				close(figs(3));
			end
		end
		if length(obj)==2
			close all
			roms1dat    = nanmean(obj.romsData.(vars{i})(1).data,3);
			roms2dat    = nanmean(obj(2).romsData.(vars{i})(1).data,3);
			[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'cmap',cmaps{i},'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1));	
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj.romsData.(vars{i})(1).units,'Interpreter','Latex');
			export_fig('-pdf',[obj.paths.plots.comp,vars{i},'_roms_diff']);
			close(figs(3));
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc longitude section plots
% P9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Diagnostics
	if A.comp == 0
        % runName loop
        for i = 1:length(runNames)
            obj = changeInputs(obj,'runName',runNames{i});
			for j = 1:length(vars);
				if ~opt(j)
					disp(['...skipping ',vars{j},'...']);
					continue
				end
				obj = clearROMS(obj);
				obj = sliceROMS(obj,vars(j),'lon',lons);
				obj = sliceDiag(obj,vars(j),'lon',lons);
				for k = 1:length(lons)
					romsdat = obj.romsData.(vars{j}).slice;
					diagdat = obj.diagData.(vars{j}).slice;
					[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{j},'xlims',xlims(k,:),'zlims',zlims,...
						'figdim',0.5,'slice',k,'levels',levs{j,k},'difflevels',dlevs{j});
					% ROMS figure
					set(0,'CurrentFigure',figs(1));
					title(['ROMS ',obj.romsData.(vars{j}).name,': ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(1),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_roms']);
					close(figs(1))
					% Diag figure
					set(0,'CurrentFigure',figs(2));
					title([obj.diagData.(vars{j}).name,': ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(2),obj.diagData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_diag']);
					close(figs(2))
					% Difference figure
					set(0,'CurrentFigure',figs(3));
					title(['Difference: ',num2str(lons(k)-360),'$^oW$'],'Interpreter','Latex');
					ylabel(cbs(3),obj.romsData.(vars{j}).units,'Interpreter','Latex');
					export_fig('-pdf',[obj.paths.plots.diag,vars{j},'_lon',num2str(lons(k)),'_diff']);
					close(figs(3))
				end
				obj.romsData = [];
				obj.diagData = [];
			end
		end
	% Inter-run comparisons
	elseif A.comp == 1
        for i = 1:length(vars)
			obj = clearROMS(obj);
            if ~opt(i)
                disp(['...skipping ',vars{i},'...']);
                continue
            end
            for j = 1:length(runNames)
                obj = changeInputs(obj,'runName',runNames{j});
				obj = sliceROMS(obj,vars(i),'lon',lons);
                tmp{j} = obj.romsData.(vars{i}).slice;
                if j == 1
                    comppath = obj.paths.plots.comp;
                end
            end
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
			for j = 1:length(lons)
				[figs,cbs] = sliceCmp(obj,roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims(j,:),'zlims',zlims,...
					'figdim',0.5,'slice',j,'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[comppath,vars{i},'_lon',num2str(lons(j)),'_roms_diff']);
				close(figs(3))
			end
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface chla
% P10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	for i = 1:length(obj)
		obj = loadData(obj,vars,'type','raw');
	end
	if A.comp == 0
		obj = loadDiag(obj,vars,0);
	end
	% Reduce data
	for i = 1:length(obj)
		% Skip validations?
		if A.comp > 0
			continue
		end
		% Reduce data
		romsdat = nanmean(obj.romsData.(vars{1})(1).data,3);
		diagdat = nanmean(obj.diagData.(vars{1})(1).data,3);
		diffdat = romsdat - diagdat;
		romsdat = real(log10(romsdat));
		diagdat = real(log10(diagdat));
		diffdat = romsMaster.dfloglevs(diffdat,0.01);
		% Plot
		[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj.romsData.(vars{1}).name,': sfc'],'Interpreter','Latex');
		ylabel(cbs(1),obj.romsData.(vars{1}).units,'Interpreter','Latex');
		cbs(1).XTickLabel = absLbls; 
		export_fig('-pdf',[obj.paths.plots.diag,vars{1},'_roms']);
		close(figs(1));
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title([obj.diagData.(vars{1}).name,': sfc'],'Interpreter','Latex');
		ylabel(cbs(2),obj.diagData.(vars{1}).units,'Interpreter','Latex');
		cbs(2).XTickLabel = absLbls; 
		export_fig('-pdf',[obj.paths.plots.diag,vars{1},'_diag']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj.romsData.(vars{1}).units,'Interpreter','Latex');
		cbs(3).XTick = diffLevs;
		cbs(3).XTickLabel = diffLbls; 
		cbs(3).Limits = diffCaxis;
		export_fig('-pdf',[obj.paths.plots.diag,vars{1},'_diff']);
		close(figs(3));
	end
	if length(obj)==2
		roms1dat = nanmean(obj.romsData.(vars{1})(1).data,3);
		roms2dat = nanmean(obj(2).romsData.(vars{1})(1).data,3);
		diffdat = roms1dat - roms2dat;
		roms1dat = real(log10(roms1dat));
		roms2dat = real(log10(roms2dat));
		diffdat = romsMaster.dfloglevs(diffdat,0.001);
		% Plot
		[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
		close(figs(1));
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj.romsData.(vars{1}).units,'Interpreter','Latex');
        cbs(3).XTick = diffLevs;
        cbs(3).XTickLabel = diffLbls;
        cbs(3).Limits = diffCaxis;
		export_fig('-pdf',[obj.paths.plots.comp,vars{1},'_roms_diff']);
		close(figs(3));
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMZ thickness
% P11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
    % Get OMZ thickness
	for i = 1:length(obj)
		obj = OMZthick(obj,omzthresh,1);
	end
    % Make comparison plots
	for i = 1:length(obj.diagData.OMZ);
		for j = 1:length(omzthresh)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				romsdat = nanmean(squeeze(obj.romsData.OMZ.int(:,:,:,j)),3);
				diagdat = squeeze(obj.diagData.OMZ(i).int(:,:,j));
				[figs,cbs] = mapCmp(obj,romsdat,diagdat,'levels',levs{j},'difflevels',dlevs{j});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS: ',obj.romsData.OMZ.name,'($O_2$ $<$ ',num2str(omzthresh(j)),' $mmol$ $m^{-3}$)'],'Interpreter','Latex');
				ylabel(cbs(1),obj.romsData.OMZ.units,'Interpreter','Latex')
				set(gcf,'ColorMap',cmap);
				export_fig('-pdf',[obj.paths.plots.diag,'OMZ_roms_th',num2str(omzthresh(j))]);
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj.diagData.OMZ(i).name],'Interpreter','Latex');
				ylabel(cbs(2),obj.diagData.OMZ(i).units,'Interpreter','Latex')
				set(gcf,'ColorMap',cmap);
				export_fig('-pdf',[obj.paths.plots.diag,'OMZ_diag_',num2str(i),'_th',num2str(omzthresh(j))]);
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj.romsData.OMZ.units,'Interpreter','Latex')
				export_fig('-pdf',[obj.paths.plots.diag,'OMZ_diff_',num2str(i),'_th',num2str(omzthresh(j))]);
				close all
			end
		end
	end
	if length(obj)==2
		for j = 1:length(omzthresh)
			roms1dat = squeeze(obj.romsData.OMZ.int(:,:,j));
			roms2dat = squeeze(obj(2).romsData.OMZ.int(:,:,j));
			[figs,cbs] = mapCmp(obj,roms1dat,roms2dat,'levels',levs{j},'difflevels',dlevs{j});
			close(figs(1));
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj.romsData.OMZ.units,'Interpreter','Latex')
			export_fig('-pdf',[obj.paths.plots.comp,'OMZ_roms_diff_th',num2str(omzthresh(j))]);
			close all
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ncycle subsurface comparisons
% P12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Call data_locations
	if exist([obj.paths.plots.diag,'prof_data_locations.pdf'])==0
		fig = data_locations(obj,obs_regions,roms_regions);
		export_fig('-pdf',[obj.paths.plots.diag,'prof_data_locations']);
		close all
	end
	% Call extract_obs	
    for i = 1:length(obs_regions)
		if exist(['data/obs_tracers_region_',region_name{i},'.mat'])==0
			disp('Extracting obs data');
			extract_obs_ncycle(obs_regions{i},region_name{i});
		end
    end
	% Call extract roms
	for i = 1:length(roms_regions)
		if exist(['data/roms_tracers_region_',region_name{i},'.mat'])==0
			disp('Extracting ROMS data');
			extract_roms_ncycle(obj,vars,roms_regions{i},region_name{i});
		end
	end
	% Call obs_vs_roms
	obs_vs_roms
	% Call obs_vs_roms_ratios
	obs_vs_roms_ratios
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROMS vs OCIM N2O 
% P13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Call ocim_locations
	if exist([obj.paths.plots.diag,'ocim_locations.pdf'])==0
		fig = ocim_locations(obj,regions);
		export_fig('-pdf',[obj.paths.plots.diag,'ocim_locations']);
		close all
	end
	% Extract OCIM	
	for i = 1:length(regions)
		OCIM_n2o_yield(obj,regions{i},i);
		ROMS_n2o_yield(obj,regions{i},i);
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OTHER DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps
% P14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).choice;
	unpactStruct(plots(pltcnt));
	tmpfields = fields(A);
	for i = 1:length(tmpfields);
		if ~isempty(A.(tmpfields{i}))
			eval([tmpfields{i},'=A.(tmpfields{i});']);
		end
	end
	% Get lon lines
	for i = 1:length(lons)
		lony{i} = [-90:0.1:90];
		lonx{i} = [lons(i)*ones(size(lony{i}))];
	end
	% Get lat lines
	for i = 1:length(lats)
		latx{i} = [0:0.1:360];
		laty{i} = [lats(i)*ones(size(latx{i}))];
	end
	% Plot
	[fig] = quickMap(obj,'ticks',1,'fontsize',8);
	hold on
	m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'k','linewidth',2);
	for i = 1:length(lons)
		[in,~] = inpolygon(lonx{i},lony{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
		m_plot(lonx{i}(in==1),lony{i}(in==1),'--k');
	end
	for i = 1:length(lats)
		[in,~] = inpolygon(latx{i},laty{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
		m_plot(latx{i}(in==1),laty{i}(in==1),'--k');
	end
	title(['Location of depth slices'],'Interpreter','Latex');
	export_fig('-pdf',[obj.paths.plots.diag,'trans_locations']);
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj A varargin plots pltcnt
end

