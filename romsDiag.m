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
%			... comparisons will be obj(1) - obj(2)
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
%	------------------------------');
%	------ OTHER DIAGNOSTICS -----');
%	------------------------------');
%	13 == slice degree maps');
%
% Optional Inputs (varargin);
% - comp = Only perform inter-run comparisons (1)
%
%
% Example:
% - obj = romsDiag(obj,[1:11]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle

% Process optional inputs
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
for i = 1:length(obj)
	% initiate plots
	obj(i) = initDiag(obj(i));
	% clear variables
	obj(i) = clearROMS(obj(i));
end

% Get plotchoice
tmpchoice = zeros(1,100);
tmpchoice(plotchoice) = 1;
plotchoice = tmpchoice;

% Load diag options
if strcmp(obj(1).info.r_tag,'pacmed_0p25');
	run /data/project1/demccoy/ROMS/pacmed_0p25/analysis/diag/pacmed_diag_opts.m
elseif strcmp(obj(1).info.r_tag,'peru_chile_0p1');
	run /data/project1/demccoy/ROMS/peru_chile_0p1/analysis/diag/peru_diag_opts.m
end

% Start pltcnt
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
	for i = 1:length(obj)
		obj(i) = loadData(obj(i),vars(opt==1),'type','z_avg','depth',zdeps);
	end
	if A.comp == 0
		obj(1) = loadDiag(obj(1),vars(opt==1),zdeps);
	end
	for i = 1:length(vars)
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(zdeps)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(squeeze(obj(k).romsData.(vars{i})(1).data(:,:,j,:)),3);
				diagdat    = nanmean(squeeze(obj(1).diagData.(vars{i})(1).data(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj(k),romsdat,diagdat,'cmap',cmaps{i},'levels',levs{i,j},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i}).name,': ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i}).name,': ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_diag']);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_diff']);
			end
			if length(obj)==2
				close all
				roms1dat    = nanmean(squeeze(obj(1).romsData.(vars{i})(1).data(:,:,j,:)),3);
				roms2dat    = nanmean(squeeze(obj(2).romsData.(vars{i})(1).data(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{i},'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_z',num2str(zdeps(j)),'_roms_diff']);
				close(figs(3));
			end
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	for i = 1:length(obj)
		obj(i) = sliceROMS(obj(i),vars(opt==1),'lat',lats);
	end
	if A.comp == 0
		obj(1) = sliceDiag(obj(1),vars(opt==1),'lat',lats);
	end
	for i = 1:length(vars);
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(obj)
			% Skip validations?
			if A.comp > 0
				continue
			end
			romsdat = obj(j).romsData.(vars{i}).slice;
			diagdat = obj(1).diagData.(vars{i}).slice;
			[figs,cbs] = sliceCmp(obj(j),romsdat,diagdat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'levels',levs{i},'difflevels',dlevs{i});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj(j).romsData.(vars{i}).name,': ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(j).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_roms']);
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj(1).diagData.(vars{i}).name,': ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_diag']);
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(j).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_diff']);
			close(figs(3))
		end
		if length(obj)==2
			roms1dat = obj(1).romsData.(vars{i}).slice;
			roms2dat = obj(2).romsData.(vars{i}).slice;
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1));
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_lat',num2str(lats),'_roms_diff']);
			close(figs(3))
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
		obj(i) = loadData(obj(i),vars(opt==1),'type','raw');
	end
	if A.comp == 0
		obj(1) = loadDiag(obj(1),vars(opt==1),0);
	end
	for i = 1:length(vars)
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(obj(1).diagData.(vars{i}))
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(obj(k).romsData.(vars{i})(1).data,3);
				diagdat    = nanmean(obj(1).diagData.(vars{i})(j).data,3);
				[figs,cbs] = mapCmp(obj(k),romsdat,diagdat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i})(1).name],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i})(j).name],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i})(j).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_diag_',num2str(j)]);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_diff_',num2str(j)]);
				close(figs(3));
			end
		end
		if length(obj)==2
			close all
			roms1dat    = nanmean(obj(1).romsData.(vars{i})(1).data,3);
			roms2dat    = nanmean(obj(2).romsData.(vars{i})(1).data,3);
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1));
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).romsData.(vars{i})(1).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_roms_diff']);
			close(figs(3));
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	for i = 1:length(obj)
		obj(i) = sliceROMS(obj(i),vars(opt==1),'lon',lons);
	end
	if A.comp == 0
		obj(1) = sliceDiag(obj(1),vars(opt==1),'lon',lons);
	end
	for i = 1:length(vars);
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(lons)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				% Grab data
				romsdat = obj(k).romsData.(vars{i}).slice;
				diagdat = obj(1).diagData.(vars{i}).slice;
				[figs,cbs] = sliceCmp(obj(k),romsdat,diagdat,0,'cmap',cmaps{i},'zlims',zlims,'xlims',xlims(j,:),...
					'figdim',0.5,'slice',j,'levels',levs{i,j},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i}).name,': ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_roms']);
				close(figs(1))
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i}).name,': ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_diag']);
				close(figs(2))
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_diff']);
				close(figs(3))
			end
			if length(obj)==2
				% Grab data
				roms1dat = squeeze(obj(1).romsData.(vars{i}).slice(:,:,:,j));
				roms2dat = squeeze(obj(2).romsData.(vars{i}).slice(:,:,:,j));
				[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{i},'zlims',zlims,'xlims',xlims(j,:),...
					'figdim',0.5,'slice',j,'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));
				close(figs(2));
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_lon',num2str(lons(j)),'_roms_diff']);
				close(figs(3))
			end
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
		obj(i) = equatorUcmp(obj(i),'179E_160W');
	end
	for i = 1:length(obj)
		% Skip validations?
		if A.comp > 0
			continue
		end
		diagdat = obj(i).diagData.u.slice;
		romsdat = nanmean(obj(i).romsData.u.slice,3);
		[figs,cbs] = sliceCmp(obj(i),romsdat,diagdat,1,'zlims',zlims,'xlims',xlims,...
			'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj(i).romsData.(vars{1}).name,': ',num2str(obj(i).slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(1),obj(i).romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj(i).paths.plots.diag,'equ_roms']);
		close(figs(1))
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title([obj(i).diagData.(vars{1}).name,': ',num2str(obj(i).slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(2),obj(i).diagData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj(i).paths.plots.diag,'equ_diag']);
		close(figs(2))
		% Difference figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference: ',num2str(obj(i).slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(3),obj(i).romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj(i).paths.plots.diag,'equ_diff']);
		close(figs(3))
	end
	if length(obj)==2
		roms1dat = obj(1).romsData.u.slice;
		roms2dat = obj(2).romsData.u.slice;
		[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'zlims',zlims,'xlims',xlims,...
			'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
		close(figs(1));
		close(figs(2));
		% Difference figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference: ',num2str(obj(1).slice.sect-360),'$^oW$'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).romsData.(vars{1}).units,'Interpreter','Latex');
		export_fig('-pdf',[obj(1).paths.plots.comp,'equ_roms_diff']);
		close(figs(3))
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	for i = 1:length(obj)
		obj(i) = loadData(obj(i),vars(opt==1),'type','z_avg','depth',zdeps);
	end
	if A.comp == 0
		obj(1) = loadDiag(obj(1),vars(opt==1),zdeps);
	end
	for i = 1:length(vars)
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(zdeps)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(squeeze(obj(k).romsData.(vars{i})(1).data(:,:,j,:)),3);
				diagdat    = nanmean(squeeze(obj(1).diagData.(vars{i})(1).data(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj(k),romsdat,diagdat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i,j},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i}).name,': ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i}).name,': ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_diag']);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_z',num2str(zdeps(j)),'_diff']);
				close(figs(3));
			end
			if length(obj)==2
				close all
				roms1dat    = nanmean(squeeze(obj(1).romsData.(vars{i})(1).data(:,:,j,:)),3);
				roms2dat    = nanmean(squeeze(obj(2).romsData.(vars{i})(1).data(:,:,j,:)),3);
				[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{i},'bal',bal(i),'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1));	
				close(figs(2));
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(zdeps(j)),'m'],'Interpreter','Latex');
				ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_z',num2str(zdeps(j)),'_roms_diff']);
				close(figs(3));
			end
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	for i = 1:length(obj)
		obj(i) = sliceROMS(obj(i),vars(opt==1),'lat',lats);
	end
	if A.comp == 0
		obj(1) = sliceDiag(obj(1),vars(opt==1),'lat',lats);
	end
    for i = 1:length(vars);
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
		for j = 1:length(obj)
			% Skip validations?
			if A.comp > 0
				continue
			end
			romsdat = obj(j).romsData.(vars{i}).slice;
			diagdat = obj(1).diagData.(vars{i}).slice;
			[figs,cbs] = sliceCmp(obj(j),romsdat,diagdat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj(j).romsData.(vars{i}).name,': ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(j).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_roms']);
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj(1).diagData.(vars{i}).name,': ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_diag']);
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(j).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(j).paths.plots.diag,vars{i},'_lat',num2str(lats),'_diff']);
			close(figs(3))
		end
		if length(obj)==2
			roms1dat = obj(1).romsData.(vars{i}).slice;
			roms2dat = obj(2).romsData.(vars{i}).slice;
			[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{i},'xlims',xlims,'zlims',zlims,...
				'figdim',0.5,'bal',bal(i),'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1))
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_lat',num2str(lats),'_roms_diff']);
			close(figs(3))
		end
    end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
		obj(i) = loadData(obj(i),vars,'type','raw');
	end
	if A.comp == 0
		obj(1) = loadDiag(obj(1),vars,0);
	end
	for i = 1:length(vars)
		for j = 1:length(obj(1).diagData.(vars{i}))
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				close all
				romsdat    = nanmean(obj(k).romsData.(vars{i})(1).data,3);
				diagdat    = nanmean(obj(1).diagData.(vars{i})(j).data,3);
				[figs,cbs] = mapCmp(obj(k),romsdat,diagdat,'cmap',cmaps{i},'levels',levs{i},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i})(1).name],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_roms']);
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i})(j).name],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i})(j).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_diag_',num2str(j)]);
				close(figs(2));	
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i})(1).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_diff_',num2str(j)]);
				close(figs(3));
			end
		end
		if length(obj)==2
			close all
			roms1dat    = nanmean(obj(1).romsData.(vars{i})(1).data,3);
			roms2dat    = nanmean(obj(2).romsData.(vars{i})(1).data,3);
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{i},'levels',levs{i},'difflevels',dlevs{i});
			close(figs(1));	
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).romsData.(vars{i})(1).units,'Interpreter','Latex');
			export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_roms_diff']);
			close(figs(3));
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	for i = 1:length(obj)
		obj(i) = sliceROMS(obj(i),vars(opt==1),'lon',lons);
	end
	if A.comp == 0
		obj(1) = sliceDiag(obj(1),vars(opt==1),'lon',lons);
	end
    for i = 1:length(vars);
		if ~opt(i)
			disp(['...skipping ',vars{i},'...']);
			continue
		end
        for j = 1:length(lons)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				% Grab data
				romsdat = obj(k).romsData.(vars{i}).slice;
				diagdat = obj(1).diagData.(vars{i}).slice;
				[figs,cbs] = sliceCmp(obj(k),romsdat,diagdat,0,'cmap',cmaps{i},'zlims',zlims,'xlims',xlims(j,:),...
					'figdim',0.5,'slice',j,'bal',bal(i),'levels',levs{i,j},'difflevels',dlevs{i});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj(k).romsData.(vars{i}).name,': ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_roms']);
				close(figs(1))
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(1).diagData.(vars{i}).name,': ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(2),obj(1).diagData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_diag']);
				close(figs(2))
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(k).paths.plots.diag,vars{i},'_lon',num2str(lons(j)),'_diff']);
				close(figs(3))
			end
			if length(obj)==2
				% Grab data
				roms1dat = squeeze(obj(1).romsData.(vars{i}).slice(:,:,:,j));
				roms2dat = squeeze(obj(2).romsData.(vars{i}).slice(:,:,:,j));
				[figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{i},'zlims',zlims,'xlims',xlims(j,:),...
					'figdim',0.5,'slice',j,'bal',bal(i),'levels',levs{i,j},'difflevels',dlevs{i});
				close(figs(1))
				close(figs(2))
				% Difference figure
				set(0,'CurrentFigure',figs(3));
				title(['ROMS Difference: ',num2str(lons(j)-360),'$^oW$'],'Interpreter','Latex');
				ylabel(cbs(3),obj(1).romsData.(vars{i}).units,'Interpreter','Latex');
				export_fig('-pdf',[obj(1).paths.plots.comp,vars{i},'_lon',num2str(lons(j)),'_roms_diff']);
				close(figs(3))
			end
        end
    end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
		obj(i) = loadData(obj(i),vars,'type','raw');
	end
	if A.comp == 0
		obj(1) = loadDiag(obj(1),vars,0);
	end
	% Reduce data
	for i = 1:length(obj)
		% Skip validations?
		if A.comp > 0
			continue
		end
		% Reduce data
		romsdat = nanmean(obj(i).romsData.(vars{1})(1).data,3);
		diagdat = nanmean(obj(1).diagData.(vars{1})(1).data,3);
		diffdat = romsdat - diagdat;
		romsdat = real(log10(romsdat));
		diagdat = real(log10(diagdat));
		diffdat = romsMaster.dfloglevs(diffdat,0.01);
		% Plot
		[figs,cbs] = mapCmp(obj(i),romsdat,diagdat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj(i).romsData.(vars{1}).name,': sfc'],'Interpreter','Latex');
		ylabel(cbs(1),obj(i).romsData.(vars{1}).units,'Interpreter','Latex');
		cbs(1).XTickLabel = absLbls; 
		export_fig('-pdf',[obj(i).paths.plots.diag,vars{1},'_roms']);
		close(figs(1));
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title([obj(1).diagData.(vars{1}).name,': sfc'],'Interpreter','Latex');
		ylabel(cbs(2),obj(1).diagData.(vars{1}).units,'Interpreter','Latex');
		cbs(2).XTickLabel = absLbls; 
		export_fig('-pdf',[obj(i).paths.plots.diag,vars{1},'_diag']);
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(i).romsData.(vars{1}).units,'Interpreter','Latex');
		cbs(3).XTick = diffLevs;
		cbs(3).XTickLabel = diffLbls; 
		cbs(3).Limits = diffCaxis;
		export_fig('-pdf',[obj(i).paths.plots.diag,vars{1},'_diff']);
		close(figs(3));
	end
	if length(obj)==2
		roms1dat = nanmean(obj(1).romsData.(vars{1})(1).data,3);
		roms2dat = nanmean(obj(2).romsData.(vars{1})(1).data,3);
		diffdat = roms1dat - roms2dat;
		roms1dat = real(log10(roms1dat));
		roms2dat = real(log10(roms2dat));
		diffdat = romsMaster.dfloglevs(diffdat,0.001);
		% Plot
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
		close(figs(1));
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).romsData.(vars{1}).units,'Interpreter','Latex');
        cbs(3).XTick = diffLevs;
        cbs(3).XTickLabel = diffLbls;
        cbs(3).Limits = diffCaxis;
		export_fig('-pdf',[obj(1).paths.plots.comp,vars{1},'_roms_diff']);
		close(figs(3));
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
		obj(i) = OMZthick(obj(i),omzthresh,1);
	end
    % Make comparison plots
	for i = 1:length(obj(1).diagData.OMZ);
		for j = 1:length(omzthresh)
			for k = 1:length(obj)
				% Skip validations?
				if A.comp > 0
					continue
				end
				romsdat = nanmean(squeeze(obj(k).romsData.OMZ.int(:,:,:,j)),3);
				diagdat = squeeze(obj(k).diagData.OMZ(i).int(:,:,j));
				[figs,cbs] = mapCmp(obj(k),romsdat,diagdat,'levels',levs{j},'difflevels',dlevs{j});
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS: ',obj(k).romsData.OMZ.name,'($O_2$ $<$ ',num2str(omzthresh(j)),' $mmol$ $m^{-3}$)'],'Interpreter','Latex');
				ylabel(cbs(1),obj(k).romsData.OMZ.units,'Interpreter','Latex')
				set(gcf,'ColorMap',cmap);
				export_fig('-pdf',[obj(k).paths.plots.diag,'OMZ_roms_th',num2str(omzthresh(j))]);
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj(k).diagData.OMZ(i).name],'Interpreter','Latex');
				ylabel(cbs(2),obj(k).diagData.OMZ(i).units,'Interpreter','Latex')
				set(gcf,'ColorMap',cmap);
				export_fig('-pdf',[obj(k).paths.plots.diag,'OMZ_diag_',num2str(i),'_th',num2str(omzthresh(j))]);
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj(k).romsData.OMZ.units,'Interpreter','Latex')
				export_fig('-pdf',[obj(k).paths.plots.diag,'OMZ_diff_',num2str(i),'_th',num2str(omzthresh(j))]);
				close all
			end
		end
	end
	if length(obj)==2
		for j = 1:length(omzthresh)
			roms1dat = squeeze(obj(1).romsData.OMZ.int(:,:,j));
			roms2dat = squeeze(obj(2).romsData.OMZ.int(:,:,j));
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'levels',levs{j},'difflevels',dlevs{j});
			close(figs(1));
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).romsData.OMZ.units,'Interpreter','Latex')
			export_fig('-pdf',[obj(1).paths.plots.comp,'OMZ_roms_diff_th',num2str(omzthresh(j))]);
			close all
		end
	end
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
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
	% Reduce inputs?
	vars = vars(opt==1);
	units = units(opt==1);
	xlims = xlims(opt==1);
	% Call data_locations
	fig = data_locations(obj,obs_regions,roms_regions);
	export_fig('-pdf',[obj(1).paths.plots.diag,'prof_data_locations']);
	% Call obs_vs_roms
	for i = 1:length(obs_regions)
		disp(['Region ',num2str(i)])
		[fig] = obs_vs_roms(obj,vars,obs_regions{i},roms_regions{i},xlims,units);
		suptitle(['Region ',num2str(i)]);
		export_fig('-pdf',[obj(1).paths.plots.diag,'prof_region_',num2str(i)]);
		close all
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OTHER DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps
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
	[fig] = quickMap(obj(1),'ticks',1,'fontsize',8);
	hold on
	m_plot(obj(1).region.polygon(:,1),obj(1).region.polygon(:,2),'k','linewidth',2);
	for i = 1:length(lons)
		[in,~] = inpolygon(lonx{i},lony{i},obj(1).region.polygon(:,1),obj(1).region.polygon(:,2));
		m_plot(lonx{i}(in==1),lony{i}(in==1),'--k');
	end
	for i = 1:length(lats)
		[in,~] = inpolygon(latx{i},laty{i},obj(1).region.polygon(:,1),obj(1).region.polygon(:,2));
		m_plot(latx{i}(in==1),laty{i}(in==1),'--k');
	end
	title(['Location of depth slices'],'Interpreter','Latex');
	export_fig('-pdf',[obj(1).paths.plots.diag,'trans_locations']);
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
	clearvars -except obj A varargin plots pltcnt
end

