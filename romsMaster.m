%------------------------------------------------------------------------------------------------
classdef romsMaster
%-------------------
% Class used to load ROMS data and make various plots,
% comparison figures against validation data, and many
% other actions.
%
% Do begin, define the object and initialize
% - obj = romsMaster;
% - obj = init(obj,'sim','peru')
%
% To view available routines
% - methods(obj)
%
% For help
% - help obj.(method)
%------------------
	%----------------------------------------------------------------------------------------
	properties
		grid % contains grid data and extra fields
		info % contains variable information
		region % contains grid data in defined region
		budget % contains parameters for budget calculations
      		romsData % contains ROMS data for comparisons
      		diagData % contains validation data for comparions
      		paths % paths to data, directories, etc
      		plots % contains info for plots (where validation data is from etc) 
   	end
	%----------------------------------------------------------------------------------------
   	methods
		%--------------------------------------------------------------------------------
		function obj = init(obj,sim,varargin)
        		% -------------------
        		% Initialization method: gathers paths and coordinate variables 
			%
                        % Usage: 
                        % - obj = init(obj,sim)
			%
			% Inputs:
			% - sim = ROMS simulation (peru or pacmed only)
			%
			% Optional Inputs:
			% - region = 'full' will override default region in current_sims.mat
			%            'manual' will not auto-run defineRegion (do it after 
			%            calling init(obj,...))
                        %
			% Example:
			% - obj = init(obj,'peru')         <-- if obj defined
			% - obj = init(romsMaster,'peru')  <-- if obj undefined
			% -------------------
	               	 
			% - Begin
			disp(' ');
                        disp('---------------------------------');
                        disp('---------------------------------');
                        disp('            romsMaster           ');
                        disp('---------------------------------');
                        disp('---------------------------------');
                        disp(' ');
		
                        % - Suppress warning messages and addpath to Danny's scripts
                        warning off
                        addpath /data/project1/demccoy/matlab_scripts/
		
			% - Run current_sims
			current_sims;

			% - Check inputs
			if ~isempty(sim) % sim
				if ~strcmp(sim,'peru') & ~strcmp(sim,'pacmed') & ~strcmp(sim,'pacmed_2006') & ~strcmp(sim,'pacmed_2010');
					disp(' ');
					disp('sim must be peru, pacmed, pacmed06, or pacmed10');
					disp(' '); return
				end	
			end

			% - Optional inputs...
			A.region = [];
                        A        = parse_pv_pairs(A,varargin);
	
			% - Choose simulation if no input, then load info based on choice/input
			if isempty(sim)
                        	sim = input('---------------------------\nperu_chile(1) or pacmed06(2) or pacmed10(3)\n---------------------------\n>> ');
				if ~ismember(sim,[1 2 3])
					disp(' ');
					disp('Only 1,2,3 allowed');
					disp(' '); return
				end
			end
			if sim == 1
				sim = 'peru';
			elseif sim == 2
				sim = 'pacmed_2006';
			elseif sim == 3
				sim = 'pacmed_2010';
			end
			if strcmp(sim,'peru');
                                load([RMdir,'peru_sim.mat']);
                        elseif strcmp(sim,'pacmed_2006');
                                load([RMdir,'pacmed_sim_06.mat']);
                        elseif strcmp(sim,'pacmed_2010') | strcmp(sim,'pacmed');
                                load([RMdir,'pacmed_sim_10.mat']);
                        end
      
			% - grab paths according to inputs
			obj.paths.simPath = simPath; % root path

      			% - grab file paths
         		obj.paths.avg    = [obj.paths.simPath,'avg/',avgfile];
			obj.paths.his    = [obj.paths.simPath,'his/','his_',ext,'.nc'];
			obj.paths.flux   = [obj.paths.simPath,'phys_flux/','phys_flux_avg_',ext,'.nc'];
         		obj.paths.zavg 	 = [obj.paths.simPath,'z_avg/','z_',avgfile];
			obj.paths.grid   = gfile;

      			% - initiate directories if they dont exist
        		[pathstr, name, ext] = fileparts(obj.paths.avg);
			mkdir([obj.paths.simPath,'budgetFigures']);
            		mkdir([obj.paths.simPath,'diagFigures']);
            		mkdir([obj.paths.simPath,'diagFigures/',name]);
            		mkdir([obj.paths.simPath,'diagFigures/',name,'/surface_figs']);
            		mkdir([obj.paths.simPath,'diagFigures/',name,'/z_section_figs']);
			mkdir([obj.paths.simPath,'diagFigures/',name,'/z_surface_figs']);
            		mkdir([obj.paths.simPath,'diagFigures/',name,'/dens_section_figs']);
            		mkdir([obj.paths.simPath,'diagFigures/',name,'/int_figs']);
            		mkdir([obj.paths.simPath,'diagFigures/',name,'/ncycle_figs']);
			mkdir([obj.paths.simPath,'diagFigures/',name,'/mld_figs']);
			mkdir([obj.paths.simPath,'diagFigures/',name,'/omz_figs']);
                        mkdir([obj.paths.simPath,'romsFigures']);
                        mkdir([obj.paths.simPath,'romsFigures/',name]);
                        mkdir([obj.paths.simPath,'romsFigures/',name,'/surface_figs']);
                        mkdir([obj.paths.simPath,'romsFigures/',name,'/z_section_figs']);
                        mkdir([obj.paths.simPath,'romsFigures/',name,'/z_surface_figs']);
			mkdir([obj.paths.simPath,'romsFigures/',name,'/diagnostic']);

      			% - grab plot paths
			obj.paths.plots.budget            = [obj.paths.simPath,'budgetFigures/'];
         		obj.paths.plots.diag.surfacefigs  = [obj.paths.simPath,'diagFigures/',name,'/surface_figs/'];
         		obj.paths.plots.diag.zsecfigs     = [obj.paths.simPath,'diagFigures/',name,'/z_section_figs/'];
			obj.paths.plots.diag.zsurfacefigs = [obj.paths.simPath,'diagFigures/',name,'/z_surface_figs/'];
         		obj.paths.plots.diag.dsecfigs     = [obj.paths.simPath,'diagFigures/',name,'/dens_section_figs/'];
         		obj.paths.plots.diag.intfigs      = [obj.paths.simPath,'diagFigures/',name,'/int_figs/'];
         		obj.paths.plots.diag.ncyclefigs   = [obj.paths.simPath,'diagFigures/',name,'/ncycle_figs/'];
			obj.paths.plots.diag.mldfigs      = [obj.paths.simPath,'diagFigures/',name,'/mld_figs/'];	
			obj.paths.plots.diag.omzfigs      = [obj.paths.simPath,'diagFigures/',name,'/omz_figs/'];	
                        obj.paths.plots.roms.surfacefigs  = [obj.paths.simPath,'romsFigures/',name,'/surface_figs/'];
                        obj.paths.plots.roms.zsecfigs     = [obj.paths.simPath,'romsFigures/',name,'/z_section_figs/'];
                        obj.paths.plots.roms.zsurfacefigs = [obj.paths.simPath,'romsFigures/',name,'/z_surface_figs/'];
                        obj.paths.plots.roms.diagnostic   = [obj.paths.simPath,'romsFigures/',name,'/diagnostic/'];
                        obj.paths.plots.roms.tmpfigs      = ['/data/project1/demccoy/tmpfigs/'];

        	 	% - get grid coordinates
         		obj.grid.lon_rho  = double(ncread(obj.paths.grid,'lon_rho'));
         		obj.grid.lat_rho  = double(ncread(obj.paths.grid,'lat_rho'));
                        obj.grid.pm       = double(ncread(obj.paths.grid,'pm'));
                        obj.grid.pn       = double(ncread(obj.paths.grid,'pn'));
                        obj.grid.angle    = double(ncread(obj.paths.grid,'angle'));
                        obj.grid.mask_rho = double(ncread(obj.paths.grid,'mask_rho'));
                        obj.grid.h        = double(ncread(obj.paths.grid,'h'));
                        obj.grid.area     = (1./(obj.grid.pm .* obj.grid.pn));
         		obj.grid.dx       = 1./ncread(obj.paths.grid,'pm');
         		obj.grid.dy       = 1./ncread(obj.paths.grid,'pn');
         		obj.grid.area_rho = double(obj.grid.dx.*obj.grid.dy);
			try
				obj.grid.Hz       = double(ncread(obj.paths.avg,'Hz'));
				obj.grid.z_r      = double(ncread(obj.paths.avg,'z_r'));
				obj.grid.z_w      = double(ncread(obj.paths.avg,'z_w'));
			catch; end

		       	% - get more coord diags 
			tmpinfo = ncinfo(obj.paths.avg);
			tmpinfo = tmpinfo.Dimensions;
			xi      = strmatch('xi_rho',{tmpinfo.Name});
			yi      = strmatch('eta_rho',{tmpinfo.Name});
			zi      = strmatch('s_rho',{tmpinfo.Name});
			ti      = strmatch('time',{tmpinfo.Name});
       		  	obj.grid.nx         = tmpinfo(xi).Length;
       		  	obj.grid.ny         = tmpinfo(yi).Length;
       		  	obj.grid.nz         = tmpinfo(zi).Length;
       		  	obj.grid.nt         = tmpinfo(ti).Length;
         		obj.grid.ndim_xy    = [obj.grid.nx,obj.grid.ny];
         		obj.grid.ndim_xyz   = [obj.grid.nx,obj.grid.ny,obj.grid.nz];
         		obj.grid.ndim_xyt   = [obj.grid.nx,obj.grid.ny,obj.grid.nt];
         		obj.grid.ndim_xyzt  = [obj.grid.nx,obj.grid.ny,obj.grid.nz,obj.grid.nt];
         		obj.grid.minlon_rho = nanmin(obj.grid.lon_rho(:));
         		obj.grid.maxlon_rho = nanmax(obj.grid.lon_rho(:));
         		obj.grid.minlat_rho = nanmin(obj.grid.lat_rho(:));
         		obj.grid.maxlat_rho = nanmax(obj.grid.lat_rho(:));
         		tmp                  = load('/data/project1/demccoy/ROMS/validation/wcoord.mat');
         		obj.grid.woa0p25    = tmp.wcoord0p25;
         		obj.grid.woa1p0     = tmp.wcoord1p0;
                        obj.grid.z_avg_dep  = tmp.wcoord0p25.depth;

                        % Fix masks
                        obj.grid.mask_rho(obj.grid.mask_rho==0) = NaN;

		        % - Make coord data all single
            		obj.grid = romsMaster.struct2double(obj.grid);

			% - make z_avg file if it doesn't exist
			if exist([simPath,'z_avg/z_',avgfile]) ~= 2
                        	disp('---------------------------------');
        			disp('Create z_avg Climatology File');
                       		disp('---------------------------------');
        			obj = make_zavg(obj,'gridFile',gfile,'avg',avgfile);
			end	
	
			% - get info for ROMS variables
                        % - list variables
                        obj.info = ncinfo(obj.paths.avg);
                        cnt2d = 1;
                        cnt3d = 1;
                        for i = 1:length(obj.info.Variables)
                                if length(obj.info.Variables(i).Size)==3 & obj.info.Variables(i).Size(3) == obj.grid.nt;
                                        obj.info.var2d{cnt2d}  = obj.info.Variables(i).Name;
                                        obj.info.name2d{cnt2d} = obj.info.Variables(i).Attributes(1).Value;
                                        obj.info.unit2d{cnt2d} = obj.info.Variables(i).Attributes(2).Value;
                                        obj.info.idx2d(cnt2d)  = i;
                                        cnt2d                  = cnt2d + 1;
                                elseif length(obj.info.Variables(i).Size)==4
                                        obj.info.var3d{cnt3d} = obj.info.Variables(i).Name;
                                        obj.info.name3d{cnt3d} = obj.info.Variables(i).Attributes(1).Value;
                                        obj.info.unit3d{cnt3d} = obj.info.Variables(i).Attributes(2).Value;
                                        obj.info.idx3d(cnt3d)  = i;
                                        cnt3d                  = cnt3d + 1;
                                end
                        end

			% Add some other info
			obj.info.Ext = ext;
			if strcmp(sim,'peru');
				obj.info.sim = 'peru_chile_0p1';
                        elseif strcmp(sim,'pacmed');
				obj.info.sim = 'pacmed_0p25';
                        end

                        % - Grab time info from average file
                        fieldnames = {obj.info.Attributes.Name};
                        fieldvalue = {obj.info.Attributes.Value};
                        dt         = fieldvalue(strcmp(fieldnames,'dt'));
                        dt         = double(dt{1});
                        navg       = fieldvalue(strcmp(fieldnames,'navg'));
                        navg       = double(navg{1});
			ntimes     = fieldvalue(strcmp(fieldnames,'ntimes'));
			ntimes     = double(ntimes{1});

                        % - Compute dt according to output frequency
                        obj.info.Freq   = dt*navg/86400;
			obj.info.Ntimes = ((dt*ntimes)/86400)/(obj.info.Freq);
			if obj.info.Freq > 27 & obj.info.Freq < 32
				obj.info.time_string = ['Monthly'];
			elseif obj.info.Freq == 1
				obj.info.time_string = ['Daily'];
			elseif obj.info.Freq < 1
				obj.info.time_string = ['Hourly'];
			end

			% - Save subregion indices
			obj.region.lat_lim = rlati;
			obj.region.lon_lim = rloni;
			obj.region.dep_lim = rdepi;

			% - Grab subregion
			if strcmp(A.region,'full')
				lnl = [1 obj.grid.nx];
				ltl = [1 obj.grid.ny];
				dpl = [-inf inf];
				obj = defineRegion(obj,'lon_lim',[lnl],'lat_lim',[ltl],'dep_lim',[dpl]);
			elseif strcmp(A.region,'manual');
			else
				obj = defineRegion(obj);
			end

		end % end methods init

		%--------------------------------------------------------------------------------
		function [fig,ax] = gridView(obj,varargin)
                        % ----------------------
                        % Plots the lat/lon indices of a grid file
			%
			% Usage:
			% - [fig,ax] = gridView(obj,varargin)
			%
			% Inputs:
			% - dx = plot lon lines separated by dx (default = 20)
			% - dy = plot lat lines separated by dy (default = 20)
			%
			% Example:
			% - [fig,ax] = gridView(obj,'dx',20,'dy',20)
			% ----------------------

			% Grab inputs (varargin)
                        A.dx    = [20];
			A.dy    = [20];
			A.ticks = [0];
                        A       = parse_pv_pairs(A,varargin);

			% Plot lon/lat lines
			fig(1) = piofigs('lfig',1.5);
			ax(1)  = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks);	
			[a,b]  = size(obj.grid.lon_rho);
			for i = 1:A.dx:a
				m_plot(obj.grid.lon_rho(i,:),obj.grid.lat_rho(i,:),'r');
				m_text(obj.grid.lon_rho(i,end),obj.grid.lat_rho(i,end),num2str(i),'fontsize',8);
				for j = 1:A.dy:b
					hold on
					m_plot(obj.grid.lon_rho(:,j),obj.grid.lat_rho(:,j),'b');
					m_text(obj.grid.lon_rho(end,j),obj.grid.lat_rho(end,j),num2str(j),'fontsize',8);
				end
			end

			fname = [obj.info.sim,'_grid'];
			print('-djpeg',[obj.paths.simPath,fname]);

			% Optional plot of region
			if ~isempty(obj.region.lon_rho)
				% Generate whole map
				fig(2)  = piofigs('lfig',1.5);
				set(0,'CurrentFigure',fig(2));
				[ax(2)] = map_plot(fig(2),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks);
				hold on
				% Plot grid box
				m_plot(obj.grid.lon_rho(1,:),  obj.grid.lat_rho(1,:),'k','linewidth',2);
				m_plot(obj.grid.lon_rho(:,1),  obj.grid.lat_rho(:,1),'k','linewidth',2);
				m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'k','linewidth',2);
				m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'k','linewidth',2);
				% Plot region box
				m_plot(obj.region.lon_rho(1,:),obj.region.lat_rho(1,:),'--k','linewidth',2);
				m_plot(obj.region.lon_rho(:,1),obj.region.lat_rho(:,1),'--k','linewidth',2);
				m_plot(obj.region.lon_rho(end,:),obj.region.lat_rho(end,:),'--k','linewidth',2);
				m_plot(obj.region.lon_rho(:,end),obj.region.lat_rho(:,end),'--k','linewidth',2);
				fname = [obj.info.sim,'_region'];
				print('-djpeg',[obj.paths.simPath,fname]);
			end
		end % end method gridView

		%--------------------------------------------------------------------------------
		function obj = getBudg(obj,varname,varargin)
                        % --------------------
                        % Main method to perform budget analysis on a ROMS tracer (varname).
                        % Can also toggle to Redfield ratio (R_nc, 16N/116C (default)).
                        %
                        % Usage:
                        % - obj = getBudg(obj,varargin)
                        %
                        % Required Inputs:
                        % - varname = 'N','NO3','NH4','O','C' for different budgets
			%
			% Optional Inputs:
                        % - R_nc     = N:C Redfield ratio (default is 0.137 == 16/116)
                        %
                        % Example:
                        % - obj = getBudg(obj,'N2O')
			%
			% This will instrcut budget methods (i.e. intRates, computeDzDt) on which
			% tracer to perform budget on. Here, tracking N2O.
			%
			% Available tracer budgets:
			% NO3
			% NO2
			% N2O
			% N2O_decomp
                        % --------------------

                        disp('---------------------------------');
                        disp('Get parameters');

                        % - Toggles
                        A.R_nc     = [0.137];
                        A          = parse_pv_pairs(A,varargin);

			% Check inputs
			if isempty(varname)
				disp('varname must be defined, see romsMaster.getBudg')
				return
			end

                        % Initialize?
                        if isempty(obj.grid)
                        	disp('Initialize routine first (init)')
				return
			end
			if ~isfield(obj.region,'mask_rho3d') 
				obj = defineRegion(obj);
			end

                        % Get frequency,varname
                        obj.budget.varname     = varname;

                        % Get info for zlevs4 use
                        obj.budget.param.theta_s = [6.0];
                        obj.budget.param.theta_b = [3.0];
                        obj.budget.param.hc      = [250];
                        obj.budget.param.sc_type = ['new2012'];
                        obj.budget.param.NZ      = obj.region.nz;

			% Run budget scripts in order
			obj = getBudg(obj,varname)
			% Get dzdt terms
			obj = computeDzDt(obj);
			% Get dcdt terms
			obj = computeDcDt(obj);
			% Get concentrations
			obj = getConc(obj);
			% Get rates
			obj = getRates(obj);
			% Get fluxes
			obj = getFluxes(obj);
			% Get XYZ fluxes (advection)
			obj = computeXYZflux(obj);
			% Get sources-minus-sinks
			obj = computeSMS(obj);
			% Get remainder (net)
			obj = computeNet(obj);
			% Integrate concentrations vertically
			obj = intConc(obj);
			% Integrate rates vertically
			obj = intRates(obj);
			% Integrate vertically (and horizontally)
			obj = intBudg(obj); 
		end % end method getBudg

		%--------------------------------------------------------------------------------
                function obj = defineRegion(obj,varargin)
                        % ----------------------
                        % Defines a subregion for budgets/plotting.
                        % This will reduce data to a smaller grid and depth (if toggled).
                        %
                        % Usage:
                        % - obj = defineRegion(obj,varargin)
                        %
                        % Inputs:
                        % - lon_lim = x-grid indices (i.e. [200 400]) 
                        % - lat_lim = y-grid indices (i.e. [200 400])
			% - dep_lim = z-grid depth limits (i.e. [-600 0]) for 0 - 600m
                        %
                        % Example:
                        % - obj = defineRegion(obj,'lon_lim',[200 400],'lat_lim',[200 400],'dep_lim',[-600 0])
                        %   or simply...
                        % - obj = defineRegion(obj)
                        %   ...to use defaults defined in 'current_sims.m'
                        % ----------------------

                        disp('---------------------------------');
                        disp('Defining subregion for budget');

                        close all

                        % Toggles
                        A.lon_lim = [];
                        A.lat_lim = [];
			A.dep_lim = [];
                        A         = parse_pv_pairs(A,varargin);

                        % Initialize?
                        if isempty(obj.grid)
                                obj = initBudg(obj);
                        end
                        
			% Process inputs
			% Use defaults (current_sims.m)  if calling defineRegion without inputs
			if ~isempty(A.lon_lim)
                                obj.region.lon_lim = [A.lon_lim];
			elseif isempty(obj.region.lon_lim)
				obj.region.lon_lim = [1 obj.grid.nx];
			end
                        if ~isempty(A.lat_lim)
                                obj.region.lat_lim = [A.lat_lim];
			elseif isempty(obj.region.lat_lim)
				obj.region.lat_lim = [1 obj.grid.ny];
			end
			if ~isempty(A.dep_lim)
				obj.region.dep_lim = [A.dep_lim];
			elseif isempty(obj.region.dep_lim)
				obj.region.dep_lim = [-inf inf];
			end

                        % Save region
                        obj.region.lon_rho    = double(ncread(obj.paths.grid,'lon_rho', [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.lat_rho    = double(ncread(obj.paths.grid,'lat_rho', [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.pm         = double(ncread(obj.paths.grid,'pm',      [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.pn         = double(ncread(obj.paths.grid,'pn',      [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.angle      = double(ncread(obj.paths.grid,'angle',   [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.mask_rho   = double(ncread(obj.paths.grid,'mask_rho',[obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.h          = double(ncread(obj.paths.grid,'h',       [obj.region.lon_lim(1) obj.region.lat_lim(1)],...
                                                                                        [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1]));
                        obj.region.grid_area  = (1./(obj.region.pm .* obj.region.pn));
			try
				obj.region.z_r        = obj.grid.z_r(obj.region.lon_lim(1):obj.region.lon_lim(2),...
								     obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
				obj.region.z_w        = obj.grid.z_w(obj.region.lon_lim(1):obj.region.lon_lim(2),...
								     obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
				obj.region.Hz         = obj.grid.Hz(obj.region.lon_lim(1):obj.region.lon_lim(2),...
								     obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
			catch; end
			obj.region.nx         = diff(obj.region.lon_lim)+1;
			obj.region.ny         = diff(obj.region.lat_lim)+1;
			obj.region.nz         = obj.grid.nz;
			obj.region.nt         = obj.grid.nt;
                        obj.region.ndim_xy    = [obj.region.nx,obj.region.ny];
                        obj.region.ndim_xyz   = [obj.region.nx,obj.region.ny,obj.region.nz];
                        obj.region.ndim_xyt   = [obj.region.nx,obj.region.ny,obj.region.nt];
                        obj.region.ndim_xyzt  = [obj.region.nx,obj.region.ny,obj.region.nz,obj.region.nt];
                        obj.region.minlon_rho = nanmin(obj.region.lon_rho(:));
                        obj.region.maxlon_rho = nanmax(obj.region.lon_rho(:));
                        obj.region.minlat_rho = nanmin(obj.region.lat_rho(:));
                        obj.region.maxlat_rho = nanmax(obj.region.lat_rho(:));

                        % Fix masks
                        obj.region.mask_rho(obj.region.mask_rho==0) = NaN;

                        % Make 3D grid_area, grid_mask with depth limits applied
                        for z = 1:obj.grid.nz
				for t = 1:obj.grid.nt
					% Grab temporary mask/area
					tmp_mask = obj.region.mask_rho;
					tmp_area = obj.region.grid_area;
					% Apply dep limits?
					if isempty(obj.region.dep_lim)
						obj.region.mask_rho3d(:,:,z,t) = obj.region.mask_rho;
						obj.region.area3d(:,:,z)       = obj.region.grid_area;
					else
						tmp_mask(obj.region.z_r(:,:,z,t)<obj.region.dep_lim(1)) = NaN;
						tmp_mask(obj.region.z_r(:,:,z,t)>obj.region.dep_lim(2)) = NaN;
						tmp_area(obj.region.z_r(:,:,z,t)<obj.region.dep_lim(1)) = NaN;
						tmp_area(obj.region.z_r(:,:,z,t)>obj.region.dep_lim(2)) = NaN;
						obj.region.mask_rho3d(:,:,z,t) = tmp_mask;
						obj.region.area3d(:,:,z)       = tmp_area;
					end
				end
                        end
                	
			% Make 3D mask for z_avg data
			zind        = find(-obj.grid.z_avg_dep >= obj.region.dep_lim(1) & -obj.grid.z_avg_dep <= obj.region.dep_lim(2));
			mask_rhoz3d = nan(obj.region.nx,obj.region.ny,length(obj.grid.z_avg_dep));
			for z = 1:length(obj.grid.z_avg_dep)
				if ismember(z,zind)
					mask_rhoz3d(:,:,z) = obj.region.mask_rho;
				end
			end
			obj.region.mask_rhoz3d = repmat(mask_rhoz3d,[1 1 1 12]);
			
		end % end method defineRegion

		%--------------------------------------------------------------------------------
		function obj = getConc(obj)
			% -------------------
			% Grab budget concentration(s) at each timestep
			% Called in getBudg
			% End result is mmol/m3
			%
			% Usage:
			% - obj = getConc(obj)
			% -------------------

			% Load concentration based on input
			if strcmp(obj.budget.varname,'NO3');
				disp('...loading NO3...');	
				vars = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2');
				disp('...loading NO2...');
				vars = {'NO2'};
			elseif strcmp(obj.budget.varname,'NH4');
				disp('...loading NH4...');
				vars = {'NH4'};
			elseif strcmp(obj.budget.varname,'N2O');
				disp('...loading N2O...');
				vars = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp');
				disp('...loading N2O (decomp)...');
				vars = {'N2O','N2O_AO1','N2O_ATM','N2O_SIDEN','N2O_SODEN'};
			end

			% Load
			obj = loadData(obj,'vars',vars,'type','raw');
		
			% Apply 3D mask
			for i = 1:length(vars)
				obj.romsData.(vars{i}).data = obj.romsData.(vars{i}).data .* obj.region.mask_rho3d;
			end

			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				obj.romsData.N2O_DECOMP.data = obj.romsData.N2O_SODEN.data + obj.romsData.N2O_AO1.data + ...
						      	       obj.romsData.N2O_ATM.data   + obj.romsData.N2O_SIDEN.data;
			end
		end % end method getConc(obj)

		%--------------------------------------------------------------------------------
		function obj = getRates(obj);
			% -------------------
			% Grab 3D production/loss rates for budget 
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = getRates(obj)
			% -------------------

			disp('---------------------------------');
			disp('Get production/loss rates');
		
			% Load rates based on inputs from getParam	
			if strcmp(obj.budget.varname,'N');
				disp('...full nitrogen cycle rates...');
				vars = {'AMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2','DENITRIF3',...
					'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE',...	
					'SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE',...	
					'SP_NH4_UPTAKE','DIAT_NH4_UPTAKE','DIAZ_NH4_UPTAKE',...	
					'N2OAMMOX','N2OSODEN_CONS','N2OAO1_CONS','N2OATM_CONS',...
					'N2OSIDEN_CONS','N2O_PROD_NEV','N2O_CONS_NEV','POC_REMIN',...
					'TOT_PROD','DOC_REMIN','O2_CONSUMPTION','O2_PRODUCTION',...
					'GRAZE_SP','GRAZE_DIAT','GRAZE_DIAZ','SP_LOSS','DIAT_LOSS',...
					'ZOO_LOSS','DIAZ_LOSS','DON_PROD','DON_REMIN'};
			elseif strcmp(obj.budget.varname,'NO3');
				disp('...nitrate rates only...');
				vars = {'NITROX','DENITRIF1','SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE'};	
			elseif strcmp(obj.budget.varname,'NO2');
				disp('...nitrite rates only...');
				vars = {'AMMOX','N2OAMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2',...
					'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
			elseif strcmp(obj.budget.varname,'N2O');
				disp('...nitrous oxide rates only...');
				vars = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
			elseif strcmp(obj.budget.varname,'N2O_decomp');
				disp('...nitrous oxide rates (decomp) only...');
				vars = {'DENITRIF2','N2OAMMOX','DENITRIF3','N2OSODEN_CONS','N2OSIDEN_CONS',...
					'N2OAO1_CONS','N2OATM_CONS'};
			elseif strcmp(obj.budget.varname,'NH4');
				disp('...ammonium rates only...');
				return
			elseif strcmp(obj.budget.varname,'O2');
				disp('...oxygen rates only...');
				vars = {'O2_CONSUMPTION','O2_PRODUCTION'};
			end
			
			% Load rates, restrict to region
			obj = loadData(obj,'vars',vars,'type','raw');
	
			% Apply 3D Mask
			for i = 1:length(vars)
				obj.romsData.(vars{i}).data = obj.romsData.(vars{i}).data .* obj.region.mask_rho3d;
			end	

			% Add additional rates here
			if ismember('DIAT_NO3_UPTAKE',vars) & ismember('DIAZ_NO3_UPTAKE',vars) & ismember('SP_NO3_UPTAKE',vars)
				obj.romsData.photo_NO3.data = obj.romsData.DIAT_NO3_UPTAKE.data + obj.romsData.DIAZ_NO3_UPTAKE.data + ...
							      obj.romsData.SP_NO3_UPTAKE.data;
			end
			if ismember('DIAT_NO2_UPTAKE',vars) & ismember('DIAZ_NO2_UPTAKE',vars) & ismember('SP_NO2_UPTAKE',vars)
				obj.romsData.photo_NO2.data = obj.romsData.DIAT_NO2_UPTAKE.data + obj.romsData.DIAZ_NO2_UPTAKE.data + ...
							      obj.romsData.SP_NO2_UPTAKE.data;
			end

			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				obj.romsData.DENITRIF3_DECOMP.data = obj.romsData.N2OSODEN_CONS.data + obj.romsData.N2OAO1_CONS.data + ...
							             obj.romsData.N2OATM_CONS.data   + obj.romsData.N2OSIDEN_CONS.data;
			end

		end % end method getRates	

		%--------------------------------------------------------------------------------
		function obj = getFluxes(obj);
			% -------------------
			% Grab 2D fluxes for budget, convert to 3D 
			% Called in getBudg
			% Fluxes in mmol/m2/s
			%
			% Usage:
			% - obj = getFluxes(obj)
			% -------------------
			
			disp('---------------------------------');
			disp('Get 2D fluxes')
			
			% Load fluxes based on inputs from getParam
			if strcmp(obj.budget.varname,'N')
				disp('...full nitrogen cycle fluxes...');
				return
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate fluxes only...');
				vars = {'SED_DENITRIF'};
				lvls = {'sed'};
				disp('Sediment fluxes not coded yet'); return
				obj.romsData.NO3.fg = zeros(size(obj.region.mask_rho3d));
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...no nitrite fluxes...');
				obj.romsData.NO2.fg = zeros(size(obj.region.mask_rho3d));
				return
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide fluxes only...');
				vars = {'FG_N2O'};
				lvls = {'sfc'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide fluxes (decomp) only...');
				vars = {'FG_N2O','FG_N2O_ATM','FG_N2O_SIDEN','FG_N2O_SODEN','FG_N2O_AO1'};
				lvls = {'sfc','sfc','sfc','sfc','sfc'};
			elseif strcmp(obj.budget.varname,'O2')
				disp('...oxygen fluxes only...');
				vars = {'FG_O2'};
				lvls = {'sfc'};
				return
			end	
			
			% Load fluxes, restrict to region
			obj = loadData(obj,'vars',vars,'type','raw');

			% Convert 2D flux to 3D based on lvls
			for i = 1:length(vars)
				
				% Apply 2D mask to 2D data
				obj.romsData.(vars{i}).data = obj.romsData.(vars{i}).data .* obj.region.mask_rho;
				
				% Apply 2D flux to correct z-level to make 3D
				if strcmp(lvls{i},'sfc')
					% Get fixed variable (remove 'FG_')
					vname = vars{i}(4:end);
					tmpfg = zeros(size(obj.region.mask_rho3d));
					% Apply value into 3D grid
					tmpfg(:,:,obj.region.nz,:) = obj.romsData.(vars{i}).data .* obj.region.mask_rho;
					% Divide by z, save as 3D rate
					obj.romsData.(vname).fg = tmpfg ./ obj.budget.dzdt.dz;
				elseif strcmp(lvls{i},'sed')
					% Get fixed variable (remove 'SED_')  !!! this may need hard coding depending on budget !!!
					vname  = obj.budget.varame;
					tmpsed = zeros(size(obj.region.mask_rho3d));
					% Apply value into 3D grid
					tmpsed(:,:,1,:) = obj.romsData.(vars{i}).data .* obj.region.mask_rho;
					% Divide by z, save as 3D rate
					obj.romsData.(vname).sed = tmpsed ./ obj.budget.dzdt.dz;
				end
			end

			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				% 2D version
				obj.romsData.FG_N2O_DECOMP.data = obj.romsData.FG_N2O_SODEN.data + obj.romsData.FG_N2O_AO1.data + ...
							          obj.romsData.FG_N2O_ATM.data   + obj.romsData.FG_N2O_SIDEN.data;
				% 3D version
				obj.romsData.N2O_DECOMP.fg      = obj.romsData.N2O_SODEN.fg + obj.romsData.N2O_AO1.fg + ...
							          obj.romsData.N2O_ATM.fg   + obj.romsData.N2O_SIDEN.fg;
			end

		end % end method getFluxes

		%--------------------------------------------------------------------------------
		function obj = computeDzDt(obj,varargin)
			% -------------------
			% Compute height of grid cells from his/avg files and calculates dz/dt
			% Called in getBudg
			% End results is m
			%
			% Usage:
			% - obj = computDzDt(obj)
			% -------------------

			disp('---------------------------------');
			disp('Computing changes in z-levels')
			
			% - Grab time info from average file
			stats      = ncinfo(obj.paths.avg);
			fieldnames = {stats.Attributes.Name};
			fieldvalue = {stats.Attributes.Value};
			dt         = fieldvalue(strcmp(fieldnames,'dt'));
			dt         = double(dt{1});
			navg       = fieldvalue(strcmp(fieldnames,'navg'));
			navg       = double(navg{1});

			% - Compute dt according to output frequency
			dt         = dt*navg;

			% - Calculate dzdt
			avg_z = double(ncread(obj.paths.avg,'z_w',...
				[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1],...
				[diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf inf]));
			
			% - Save levels
			obj.budget.dzdt.dz = diff(avg_z,1,3);
			obj.budget.dzdt.dt = dt;
		end % end methods computeDzDt

		%--------------------------------------------------------------------------------
		function obj = computeDcDt(obj)
			% -------------------
			% Compute change in concentration with time
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeDcDt(obj)
			% -------------------
	
			disp('---------------------------------');
			disp('Computing dC/dt')
			
			% Load concentrations based on input from getParam
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle terms...');
				vars = {'NO3','NO2','NH4','N2','N2O'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate only...');
				vars = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite only...');
				vars = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide only...');
				vars = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide (decomp) only...');
				vars = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM'};
			elseif strcmp(obj.budget.varname,'O2')
				disp('...oxygen only...');
				vars = {'O2'};
			end

			% - Load data on either side of 'history' (first,last snapshot)
			% - Scroll through each variable and get dCdt
			for v = 1:length(vars)
				
				% - Get start/finish from history
				his1 = double(ncread(obj.paths.his,vars{v},...
				       [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1],...
                                       [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf obj.region.nt]));
				
				his2 = double(ncread(obj.paths.his,vars{v},...
				       [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 2],...
                                       [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf obj.region.nt]));

				% - Divide by dt, apply 3d mask
				obj.romsData.(vars{v}).dcdt = ((his2 - his1)/obj.budget.dzdt.dt) .* obj.region.mask_rho3d;
			end 
			
			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				obj.romsData.N2O_DECOMP.dcdt = obj.romsData.N2O_SODEN.dcdt + obj.romsData.N2O_AO1.dcdt + ...
						               obj.romsData.N2O_ATM.dcdt   + obj.romsData.N2O_SIDEN.dcdt;
			end
		end % end method computeDcDt

		%--------------------------------------------------------------------------------
		function obj = computeXYZflux(obj)
			% -------------------
			% Compute HorXAdvFlux, HorYAdvFlux, and top/bottom advection
			% Also get diffusion, if it is available
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeXYZflux(obj)
			% -------------------
			
			disp('---------------------------------');
			disp('Get advective/diffusion terms')
			
			% - Get list of N-Cycle terms (can update this for other cycles)
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle advection/diffusion terms...')
				vars = {'NO3','NH4','DON','DONR','NO2','N2','N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM','N2_SED','N2O_NEV'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate advection/diffusion only...')
				vars = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite advection/diffusion only...')
				vars = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide advection/diffusion only...')
				vars = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide advection/diffusion (decomp) only...')
				vars = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM'};
			elseif strcmp(obj.budget.varname,'O2')
				disp('...oxygen advection/diffusion only...');
				vars = {'O2'};
			end

			% - Cycle through and load XYfluxes
			for f = 1:length(vars)

				% ... for each output (daily, monthly)			
				for t = 1:obj.region.nt

					% - Load XYZ advection
					temp.X  = double(ncread(obj.paths.flux,['HorXAdvFlux_',vars{f}],...
						  [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 t  ],...
					          [diff(obj.region.lon_lim)+2 diff(obj.region.lat_lim)+1 inf 1]));
					temp.Y  = double(ncread(obj.paths.flux,['HorYAdvFlux_',vars{f}],...
						  [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 t  ],...
						  [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+2 inf 1]));
					temp.Z  = double(ncread(obj.paths.flux,['VertAdvFlux_',vars{f}],...
						  [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 t  ],...
					          [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf 1]));
					temp.Zd = double(ncread(obj.paths.flux,['VertDiffFlux_',vars{f}],...
						  [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 t  ],...
						  [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf 1]));
				
					% - Initiate matrices-to-fill, compute adv X and Y
					% - Converts from mmol/s to mmol/m3/s by dividing by grid area and height of cell (volume)
					% Get dimensions
					nx = obj.region.nx;
					ny = obj.region.ny;
					nz = obj.region.nz;
					
					% - X advection
					adx(1:nx,1:ny,1:nz) = NaN;
					adx(1:nx,:,:)       = (temp.X(1:nx,:,:) - temp.X(2:nx+1,:,:)) ./ ...
							       obj.region.area3d(1:nx,:,:) ./ obj.budget.dzdt.dz(1:nx,:,:,t);
					
					% - Y advection
					ady(1:nx,1:ny,1:nz) = NaN; 		
					ady(:,1:ny,:)       = (temp.Y(:,1:ny,:) - temp.Y(:,2:ny+1,:)) ./ ...
							       obj.region.area3d(:,1:ny,:) ./ obj.budget.dzdt.dz(:,1:ny,:,t);
					
					% - Z advection
					adz(:,:,:) = (temp.Z(:,:,1:nz) - temp.Z(:,:,2:nz+1)) ./ obj.budget.dzdt.dz(:,:,:,t);
					
					% - Z diffusion
					dfz(:,:,:) = (temp.Zd(:,:,1:nz) - temp.Zd(:,:,2:nz+1)) ./ obj.budget.dzdt.dz(:,:,:,t);
					
					% - Apply 3D mask
					adx = adx .* obj.region.mask_rho3d(:,:,:,t);
					ady = ady .* obj.region.mask_rho3d(:,:,:,t);
					adz = adz .* obj.region.mask_rho3d(:,:,:,t);

					% - Save results
					tmpout.(vars{f}).adx{t} = adx;
					tmpout.(vars{f}).ady{t} = ady;
					tmpout.(vars{f}).adz{t} = adz;
					tmpout.(vars{f}).dfz{t} = dfz;				

				end
				obj.romsData.(vars{f}).adx   = cat(4,tmpout.(vars{f}).adx{:});
				obj.romsData.(vars{f}).ady   = cat(4,tmpout.(vars{f}).ady{:});
				obj.romsData.(vars{f}).adz   = cat(4,tmpout.(vars{f}).adz{:});
				obj.romsData.(vars{f}).dfz   = cat(4,tmpout.(vars{f}).dfz{:});
			end

			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				obj.romsData.N2O_DECOMP.adx  = obj.romsData.N2O_SODEN.adx + obj.romsData.N2O_AO1.adx + ...
							       obj.romsData.N2O_ATM.adx   + obj.romsData.N2O_SIDEN.adx;
				obj.romsData.N2O_DECOMP.ady  = obj.romsData.N2O_SODEN.ady + obj.romsData.N2O_AO1.ady + ...
							       obj.romsData.N2O_ATM.ady   + obj.romsData.N2O_SIDEN.ady;
				obj.romsData.N2O_DECOMP.adz  = obj.romsData.N2O_SODEN.adz + obj.romsData.N2O_AO1.adz + ...
							       obj.romsData.N2O_ATM.adz   + obj.romsData.N2O_SIDEN.adz;
				obj.romsData.N2O_DECOMP.dfz  = obj.romsData.N2O_SODEN.dfz + obj.romsData.N2O_AO1.dfz + ...
							       obj.romsData.N2O_ATM.dfz   + obj.romsData.N2O_SIDEN.dfz;
			end

		end % end methods computeXYZFlux

		%--------------------------------------------------------------------------------
		function obj = computeSMS(obj)
			% ---------------------
			% Gathers sources and sinks
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeSMS(obj)
			% ---------------------

			disp('---------------------------------');
			disp('Computing sources-minus-sinks (SMS)');

			% Get budget to close from getParam
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle SMS terms...')
				vars = {'NO3','NH4','NO2','N2','N2O'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate SMS only...')
				vars = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite SMS only...')
				vars = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O');
				disp('...nitrous oxide SMS only...')
				vars = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp');
				disp('...nitrous oxide SMS (decomp) only...')
				vars = {'N2O_decomp'};
			elseif strcmp(obj.budget.varname,'O2');
				disp('...oxygen SMS only...');
				vars = {'O2'};
			end
			
			% Scroll through each variable and get SMS
			for v = 1:length(vars)
		
				% Nitrate
				% SMS = NITROX - DENITRIF1 - SED_DENITRIF - PHOTO_NO3_UPTAKE
				if strcmp('NO3',vars{v})
					% Get SMS 
					obj.romsData.NO3.sms  = obj.romsData.NITROX.data    - ...                                  
							        obj.romsData.DENITRIF1.data - obj.romsData.photo_NO3.data - tmpsed;
					% Apply 3D mask
					obj.romsData.NO3.sms  = obj.romsData.NO3.sms .* obj.region.mask_rho3d;
				end

				% Nitrite
				% SMS = (AMMOX - 2*N2OAMMOX) - NITROX + DENITRIF1 - DENITRIF2 - ANAMMOX - PHOTO_NO2_UPTAKE
				if strcmp('NO2',vars{v})
					% Get SMS
					obj.romsData.NO2.sms  = (obj.romsData.AMMOX.data - 2*obj.romsData.N2OAMMOX.data)   + obj.romsData.DENITRIF1.data - ...
						                 obj.romsData.NITROX.data  - obj.romsData.DENITRIF2.data - ...
						                 obj.romsData.ANAMMOX.data - obj.romsData.photo_NO2.data; 
					% Apply 3D mask
					obj.romsData.NO2.sms  = obj.romsData.NO2.sms .* obj.region.mask_rho3d;
				end

				% Nitrous Oxide (basic)
				% SMS = N2OAMMOX + DENITRIF2 - DENITRIF3
				if strcmp('N2O',vars{v})
					% Get SMS
					obj.romsData.N2O.sms  = obj.romsData.N2OAMMOX.data + obj.romsData.DENITRIF2.data./2  - ...
							        obj.romsData.DENITRIF3.data;  				       
					% Apply 3D mask
					obj.romsData.N2O.sms  = obj.romsData.N2O.sms .* obj.region.mask_rho3d;
				end

				% Nitrous Oxide (decomp)
				if strcmp('N2O_decomp',vars{v})
					% Get SMS, apply mask
					obj.romsData.N2O.sms        = [obj.romsData.N2OAMMOX.data + obj.romsData.DENITRIF2.data./2  - ...
							     	       obj.romsData.DENITRIF3.data] .* obj.region.mask_rho3d;  		
					obj.romsData.N2O_AO1.sms    = [obj.romsData.N2OAMMOX.data  - ...   	
						             	       obj.romsData.N2OAO1_CONS.data] .* obj.region.mask_rho3d;    
					obj.romsData.N2O_ATM.sms    = [-obj.romsData.N2OATM_CONS.data] .* obj.region.mask_rho3d;
					obj.romsData.N2O_SIDEN.sms  = [-obj.romsData.N2OSIDEN_CONS.data] .* obj.region.mask_rho3d;
					obj.romsData.N2O_SODEN.sms  = [obj.romsData.DENITRIF2.data./2  - ...                      
							     	       obj.romsData.N2OSODEN_CONS.data] .* obj.region.mask_rho3d;					 	
					obj.romsData.N2O_DECOMP.sms = [obj.romsData.N2OAMMOX.data + obj.romsData.DENITRIF2.data./2  - ...
							               obj.romsData.DENITRIF3_DECOMP.data] .* obj.region.mask_rho3d;  			         	 
				end
	
				% Oxygen
				% SMS = O2_PRODUCTION - O2_CONSUMPTION
				if strcmp('O2',vars{v})
					% Get SMS
					obj.romsData.O2.sms = obj.romsData.O2_PRODUCTION - obj.romsData.O2_CONSUMPTION .* obj.region.mask_rho3d;
				end
			end

		end % end method getSMS

		%--------------------------------------------------------------------------------
		function obj = computeNet(obj)
			% ---------------------
			% Computes remainder (net) from budget equation
			% dC/dt = adv + dfz + sms + fg
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeNet(obj)
			% ---------------------

			% - Get list of N-Cycle terms (can update this for other cycles)
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle advection/diffusion terms...')
				vars = {'NO3','NH4','DON','DONR','NO2','N2','N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM','N2_SED','N2O_NEV'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate advection/diffusion only...')
				vars = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite advection/diffusion only...')
				vars = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide advection/diffusion only...')
				vars = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide advection/diffusion (decomp) only...')
				vars = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM'};
			elseif strcmp(obj.budget.varname,'O2')
				disp('...oxygen advection/diffusion only...');
				vars = {'O2'};
			end

			% Calculate remainder (net)
			for i = 1:length(vars)
				obj.romsData.(vars{i}).net = obj.romsData.(vars{i}).dcdt - (obj.romsData.(vars{i}).adx + ...
							     obj.romsData.(vars{i}).ady  +  obj.romsData.(vars{i}).adz + ...
							     obj.romsData.(vars{i}).dfz  +  obj.romsData.(vars{i}).sms + ...
							     obj.romsData.(vars{i}).fg);
			end

			% Optional decomp
			if strcmp(obj.budget.varname,'N2O_decomp');
				obj.romsData.N2O_DECOMP.net  = obj.romsData.N2O_SODEN.net + obj.romsData.N2O_AO1.net + ...
							       obj.romsData.N2O_ATM.net   + obj.romsData.N2O_SIDEN.net;
			end

		end % end method computeNet
	
		%--------------------------------------------------------------------------------
		function obj = intConc(obj)
                        % ------------------
                        % Vertically integrate concentrations
			% Called in getBudg
			% Terms are all in mmol/m3, get it in mmol/m2 by using dz
                        %
			% Usage:
			% - obj = Conc(obj)
                        % ------------------

			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle integration...')
				vars  = {'NO3','NH4','NO2','N2','N2O'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate integration only...')
				vars   = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite integration only...')
				vars   = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide integration only...')
				vars   = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide integration (decomp) only...')
				vars = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM','N2O_DECOMP'};
			end
		
			% Go through each 3D rate, integrate vertically
			for i = 1:length(vars)		
				tmpdata           	   = obj.romsData.(vars{i}).data .* obj.budget.dzdt.dz .* obj.region.mask_rho3d; % mmol/m3/s --> mmol/m2/s
				tmpdata          	   = squeeze(nansum(tmpdata,3));
				tmpdata                    = tmpdata .* obj.region.mask_rho;
				obj.romsData.(vars{i}).int = tmpdata; 
			end
		end % end method intRates
	
		%--------------------------------------------------------------------------------
		function obj = intRates(obj)
                        % ------------------
                        % Vertically integrate rates
			% Called in getBudg
			% Terms are all in mmol/m3/s, get it in mmol/m2/s by using dz
                        %
			% Usage:
			% - obj = intRates(obj)
                        % ------------------

			% Load rates based on inputs from getParam	
			if strcmp(obj.budget.varname,'N');
				disp('...full nitrogen cycle rates...');
				vars = {'AMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2','DENITRIF3',...
					'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE',...	
					'SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE',...	
					'SP_NH4_UPTAKE','DIAT_NH4_UPTAKE','DIAZ_NH4_UPTAKE',...	
					'N2OAMMOX','N2OSODEN_CONS','N2OAO1_CONS','N2OATM_CONS',...
					'N2OSIDEN_CONS','N2O_PROD_NEV','N2O_CONS_NEV','POC_REMIN',...
					'TOT_PROD','DOC_REMIN','O2_CONSUMPTION','O2_PRODUCTION',...
					'GRAZE_SP','GRAZE_DIAT','GRAZE_DIAZ','SP_LOSS','DIAT_LOSS',...
					'ZOO_LOSS','DIAZ_LOSS','DON_PROD','DON_REMIN'};
			elseif strcmp(obj.budget.varname,'NO3');
				disp('...nitrate rates only...');
				vars = {'NITROX','DENITRIF1','photo_NO3'};
			elseif strcmp(obj.budget.varname,'NO2');
				disp('...nitrite rates only...');
				vars = {'NITROX','AMMOX','N2OAMMOX','ANAMMOX','photo_NO2',...
					'DENITRIF1','DENITRIF2'};
			elseif strcmp(obj.budget.varname,'N2O');
				disp('...nitrous oxide rates only...');
				vars = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
			elseif strcmp(obj.budget.varname,'N2O_decomp');
				disp('...nitrous oxide rates (decomp) only...');
				vars = {'DENITRIF2','N2OAMMOX','DENITRIF3','N2OSODEN_CONS','N2OSIDEN_CONS',...
					'N2OAO1_CONS','N2OATM_CONS','DENITRIF3_DECOMP'};
			elseif strcmp(obj.budget.varname,'NH4');
				disp('...ammonium rates only...');
				return
			elseif strcmp(obj.budget.varname,'O2');
				disp('...oxygen rates only...');
				vars = {'O2_CONSUMPTION','O2_PRODUCTION'};
			end		

			% Go through each 3D rate, integrate vertically
			for i = 1:length(vars)		
				tmpdata           	   = obj.romsData.(vars{i}).data .* obj.budget.dzdt.dz .* obj.region.mask_rho3d; % mmol/m3/s --> mmol/m2/s
				tmpdata           	   = squeeze(nansum(tmpdata,3));
				tmpdata                    = tmpdata .* obj.region.mask_rho;
				obj.romsData.(vars{i}).int = tmpdata; 
			end
		end % end method intRates

		%--------------------------------------------------------------------------------
		function obj = intBudg(obj)
                        % ------------------
                        % Vertically integrate budget terms
			% Called in getBudg
			% Terms are all in mmol/m3/s, get it in mmol/m2/s by using dz
                        %
                        % Usage:
                        % - obj = intBudg(obj)
                        % ------------------

			disp('---------------------------------');
			disp('Computing vertical integration...');
			disp('...and getting grid area fluxes');
			
			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle integration...')
				vars  = {'NO3','NH4','NO2','N2','N2O'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate integration only...')
				vars   = {'NO3'};
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite integration only...')
				vars   = {'NO2'};
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide integration only...')
				vars   = {'N2O'};
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide integration (decomp) only...')
				vars = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM','N2O_DECOMP'};
			end
		
			% Go through each variables
			for i = 1:length(vars)
				
				% Grab terms, including remainder
				dcdt = obj.romsData.(vars{i}).dcdt.* obj.region.mask_rho3d;
				adx  = obj.romsData.(vars{i}).adx .* obj.region.mask_rho3d;
				ady  = obj.romsData.(vars{i}).ady .* obj.region.mask_rho3d;
				adz  = obj.romsData.(vars{i}).adz .* obj.region.mask_rho3d;
				dfz  = obj.romsData.(vars{i}).dfz .* obj.region.mask_rho3d;
				sms  = obj.romsData.(vars{i}).sms .* obj.region.mask_rho3d;
				fg   = obj.romsData.(vars{i}).fg  .* obj.region.mask_rho3d;
				net  = obj.romsData.(vars{i}).net .* obj.region.mask_rho3d;
				
				% Integrate vertically (fg term should match real flux)
				% ...mmol/m3/s to mmol/m2/s
				obj.romsData.(vars{i}).intdcdt = squeeze(nansum(dcdt.*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intadx  = squeeze(nansum(adx .*obj.budget.dzdt.dz,3)); 
				obj.romsData.(vars{i}).intady  = squeeze(nansum(ady .*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intadz  = squeeze(nansum(adz .*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intdfz  = squeeze(nansum(dfz .*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intsms  = squeeze(nansum(sms .*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intfg   = squeeze(nansum(fg  .*obj.budget.dzdt.dz,3));
				obj.romsData.(vars{i}).intnet  = squeeze(nansum(net .*obj.budget.dzdt.dz,3));

				% ...mmol/m2/s to mmol/s
				obj.romsData.(vars{i}).totdcdt  = nansum(obj.romsData.(vars{i}).intdcdt.*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totadx   = nansum(obj.romsData.(vars{i}).intadx .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totady   = nansum(obj.romsData.(vars{i}).intady .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totadz   = nansum(obj.romsData.(vars{i}).intadz .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totdfz   = nansum(obj.romsData.(vars{i}).intdfz .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totsms   = nansum(obj.romsData.(vars{i}).intsms .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totfg    = nansum(obj.romsData.(vars{i}).intfg  .*obj.region.grid_area,'all'); 
				obj.romsData.(vars{i}).totnet   = nansum(obj.romsData.(vars{i}).intnet .*obj.region.grid_area,'all'); 
			end
		end

		%--------------------------------------------------------------------------------
		function plotConc(obj,varargin)
			% ------------------
			% Plot the integrated concentrations
			% Called in getBudg
			%
			% Inputs:
			% - time = time record to plot (default--> all)
			% - prc  = percentile to limit colorbar
			%	   ...if prc == 5, then caxis = [5th %, 95th %]
			%
			% Example:
			% - plotConc(obj,'time',10,'prc',5)
			% ------------------
	
			close all

			% - defaults for optional  arguments
			A.time  = [];
			A.prc   = [2];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

 			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...plotting full nitrogen cycle concentrations...');
				return
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...plotting NO3 concentrations...');
				vars  = {'NO3'};
				units = 'mmol N m^{-2}'; 
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...plotting NO2 concentrations...');
				vars  = {'NO2'};
				units = 'mmol N m^{-2}'; 
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...plotting N2O concentrations...');
				disp('...nitrous oxide integration only...')
				vars = {'N2O'};
				units = 'mmol N_2O m^{-2}'; 
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...plotting N2O (decomp) concentrations...');
				vars = {'N2O','N2O_AO1','N2O_SIDEN', 'N2O_SODEN', 'N2O_ATM','N2O_DECOMP'};
				units = 'mmol N_2O m^{-2}'; 
			end

			% Process inputs
			if isempty(A.time)
				A.time = 1:1:obj.region.nt;
			end

			% First generate uniform color bars for each axis	
			% Go through each variables
			for i = 1:length(vars)
				cbar.tmpdata = [];
				cbar_tmpdata = [];		
				
				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	
					
					% Gather data
					tmpdata = obj.romsData.(vars{i}).int(:,:,tt); % mmol/m2/s
					
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;

					% Get colobars
					cbar_tmpdata = prclims(tmpdata,'prc',A.prc,'bal',0);
					if t == 1
						cbar.tmpdata = cbar_tmpdata;
					end
					
					% Update colorbars?
					if max(cbar_tmpdata) > max(cbar.tmpdata)
						cbar.tmpdata = cbar_tmpdata;
					end
				end
                       
				% Restart
				for t = 1:length(A.time) 
					tt = A.time(t);

					% Gather data
					tmpdata = obj.romsData.(vars{i}).int(:,:,tt); % mmol/m2/s
					
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;
					
					% Plot
					fig     = piofigs('mfig',1); set(0,'currentfigure',fig);
					if strcmp(obj.info.time_string,'Monthly')
						fname = [vars{i},'_M',num2str(tt)];
					elseif strcmp(obj.info.time_string,'Daily');
						fname = [vars{i},'_D',num2str(tt)];
					end
					clevs   = cbar.tmpdata;
					clevs   = clevs(1):(diff(clevs)/31):clevs(2);
					tmpdata(tmpdata<clevs(1))   = clevs(1);
					tmpdata(tmpdata>clevs(end)) = clevs(end);		
					[ax] = map_plot(fig(1),obj.region.lon_rho,obj.region.lat_rho);
					m_contourf(obj.region.lon_rho,obj.region.lat_rho,tmpdata,clevs,'LineStyle','none');
					cb   = colorbar;
					title(['Integrated ',vars{i}],'Interpreter','none');
					ylabel(cb,units)
					caxis([clevs(1) clevs(end)]);
					colormap(gca,cmocean('amp'));
					ax.FontSize = 10;
					cb.FontSize = 10;
					export_fig('-jpg',[obj.paths.plots.budget,fname]);
					close(fig)
				end
			end			
		end % end methods plotConc

		%--------------------------------------------------------------------------------
		function plotFluxes(obj,varargin)
			% ------------------
			% Plot the air-sea fluxes
			% Called in getBudg
			%
			% Inputs:
			% - time = time record to plot (default--> all)
			% - prc  = percentile to limit colorbar
			%	   ...if prc == 5, then caxis = [5th %, 95th %]
			%
			% Example:
			% - plotFluxes(obj,'time',10,'prc',5)
			% ------------------
	
			close all

			% - defaults for optional  arguments
			A.time  = [];
			A.prc   = [2];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

 			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle integration...')
				return
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...no nitrate fluxes...');
				return
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...no nitrite fluxes...');
				return
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide integration only...')
				vars   = {'FG_N2O'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...plotting N2O decomp fluxes...');
				vars = {'FG_N2O','FG_N2O_AO1','FG_N2O_SIDEN', 'FG_N2O_SODEN', 'FG_N2O_ATM','FG_N2O_DECOMP'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			end

			% Process inputs
			if isempty(A.time)
				A.time = 1:1:obj.region.nt;
			end

			% Units
			
			% First generate uniform color bars for each axis	
			% Go through each variables
			for i = 1:length(vars)
				cbar.tmpdata = [];
				cbar_tmpdata = [];		
				
				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	
					
					% Gather data
					tmpdata = obj.romsData.(vars{i}).data(:,:,tt); % mmol/m2/s
					
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;

					% Get colobars
					cbar_tmpdata = prclims(tmpdata,'prc',A.prc,'bal',1);
					if t == 1
						cbar.tmpdata = cbar_tmpdata;
					end
					
					% Update colorbars?
					if max(cbar_tmpdata) > max(cbar.tmpdata)
						cbar.tmpdata = cbar_tmpdata;
					end
				end
                       
				% Restart
				for t = 1:length(A.time) 
					tt = A.time(t);

					% Gather data
					tmpdata = obj.romsData.(vars{i}).data(:,:,tt); % mmol/m2/s
					
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;
					
					% Plot
					fig     = piofigs('mfig',1); set(0,'currentfigure',fig);
					if strcmp(obj.info.time_string,'Monthly')
						fname = [vars{i},'_M',num2str(tt)];
					elseif strcmp(obj.info.time_string,'Daily');
						fname = [vars{i},'_D',num2str(tt)];
					end
					clevs   = cbar.tmpdata;
					clevs   = clevs(1):(diff(clevs)/31):clevs(2);
					tmpdata(tmpdata<clevs(1))   = clevs(1);
					tmpdata(tmpdata>clevs(end)) = clevs(end);		
					[ax] = map_plot(fig(1),obj.region.lon_rho,obj.region.lat_rho);
					m_contourf(obj.region.lon_rho,obj.region.lat_rho,tmpdata,clevs,'LineStyle','none');
					cb   = colorbar;
					title(['Air-sea Flux: ',vars{i}],'Interpreter','none');
					ylabel(cb,units)
					caxis([clevs(1) clevs(end)]);
					colormap(gca,cmocean('balance'));
					ax.FontSize = 10;
					cb.FontSize = 10;
					export_fig('-jpg',[obj.paths.plots.budget,fname]);
					close(fig)
				end
			end			
		end % end methods plotFluxes

		%--------------------------------------------------------------------------------
		function plotRates(obj,varargin)
			% ------------------
			% Plot the integrated rates
			% Called in getBudg
			%
			% Inputs:
			% - time = time record to plot (default--> all)
			% - prc  = percentile to limit colorbar
			%	   ...if prc == 5, then caxis = [5th %, 95th %]
			%
			% Example:
			% - plotRates(obj,'time',10,'prc',5)
			% ------------------
	
			close all

			% - defaults for optional  arguments
			A.time  = [];
			A.prc   = [2];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

 			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle integration...')
				return
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate integration only...')
				vars = {'NITROX','DENITRIF1','photo_NO3'};
				units = 'mmol N m^{-2} s^{-1}';
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite integration only...')
				vars  = {'NITROX','AMMOX','DENITRIF1','DENITRIF2','ANAMMOX','photo_NO2','N2OAMMOX'};
				units = 'mmol N m^{-2} s^{-1}';
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide integration only...')
				vars  = {'DENITRIF2','DENITRIF3','N2OAMMOX'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...plotting N2O decomp fluxes...');
                        	vars = {'DENITRIF2','DENITRIF3','N2OAMMOX','N2OSODEN_CONS','N2OAO1_CONS', ...
                                        'N2OATM_CONS', 'N2OSIDEN_CONS','DENITRIF3_DECOMP'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			end

			% Process inputs
			if isempty(A.time)
				A.time = 1:1:obj.region.nt;
			end

			% First generate uniform color bars for each axis	
			% Go through each variables
			for i = 1:length(vars)
				cbar.tmpdata = [];
				cbar_tmpdata = [];		
				
				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	
					
					% Gather data
					if strcmp(vars{i},'DENITRIF2') & strcmp(obj.budget.varname,'N2O')
						tmpdata = obj.romsData.(vars{i}).int ./ 2;  
						tmpdata = nansum(tmpdata,3);
					elseif strcmp(vars{i},'DENITRIF2') & strcmp(obj.budget.varname,'N2O_decomp')
						tmpdata = obj.romsData.(vars{i}).int ./ 2;  
						tmpdata = nansum(tmpdata,3);
					else
						tmpdata = obj.romsData.(vars{i}).int; 
						tmpdata = nansum(tmpdata,3);
					end
					
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;

					% Get colobars
					cbar_tmpdata = prclims(tmpdata,'prc',A.prc,'bal',1);
					if t == 1
						cbar.tmpdata = cbar_tmpdata;
					end
					
					% Update colorbars?
					if max(cbar_tmpdata) > max(cbar.tmpdata)
						cbar.tmpdata = cbar_tmpdata;
					end
				end
                       
				% Restart
				for t = 1:length(A.time) 
					tt = A.time(t);

					% Gather data
					if strcmp(vars{i},'DENITRIF2') & strcmp(obj.budget.varname,'N2O')
						tmpdata = obj.int.rates.(vars{i}) ./ 2;  
						tmpdata = nansum(tmpdata,3);
					elseif strcmp(vars{i},'DENITRIF2') & strcmp(obj.budget.varname,'N2O_decomp')
						tmpdata = obj.int.rates.(vars{i}) ./ 2;  
						tmpdata = nansum(tmpdata,3);
					else
						tmpdata = obj.int.rates.(vars{i}); 
						tmpdata = nansum(tmpdata,3);
					end
       
					% Blank land
					tmpdata = tmpdata .* obj.region.mask_rho;
					
					% Plot
					fig     = piofigs('mfig',1); set(0,'currentfigure',fig);
					if strcmp(obj.info.time_string,'Monthly')
						fname = [vars{i},'_M',num2str(tt)];
					elseif strcmp(obj.info.time_string,'Daily');
						fname = [vars{i},'_D',num2str(tt)];
					end
					clevs   = cbar.tmpdata;
					clevs   = clevs(1):(diff(clevs)/31):clevs(2);
					tmpdata(tmpdata<clevs(1))   = clevs(1);
					tmpdata(tmpdata>clevs(end)) = clevs(end);		
					[ax] = map_plot(fig(1),obj.region.lon_rho,obj.region.lat_rho);
					m_contourf(obj.region.lon_rho,obj.region.lat_rho,tmpdata,clevs,'LineStyle','none');
					cb   = colorbar;
					title(['Integrated ',vars{i}],'Interpreter','none');
					ylabel(cb,units)
					caxis([clevs(1) clevs(end)]);
					colormap(gca,cmocean('balance'));
					ax.FontSize = 10;
					cb.FontSize = 10;
					export_fig('-jpg',[obj.paths.plots.budget,fname]);
					close(fig)
				end
			end			
		end % end methods plotRates

		%--------------------------------------------------------------------------------
		function plotBudg(obj,varargin);
			% ------------------
			% Plot the intergrated budget terms
			% Called in getBudg
			%
			% Inputs:
			% - time = time record to plot (default--> all)
			% - prc  = percentile to limit colorbar
			%	   ...if prc == 5, then caxis = [5th %, 95th %]
			%
			% Example:
			% - plotBudg(obj,'time',10,'prc',5)
			% ------------------

                        % - defaults for optional  arguments
                        A.time  = [];
                        A.prc   = [2];
                        A       = parse_pv_pairs(A,varargin); % parse method arguments to A
			
			% Get variables to integrate
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle plots...')
				vars  = {'NO3','NH4','NO2','N2','N2O'};
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate plots only...')
				vars   = {'NO3'};
				units  = 'mmol N/m2/s'; 
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite plots only...')
				vars   = {'NO2'};
				units  = 'mmol N/m2/s'; 
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide plots only...')
				vars   = {'N2O'};
				units  = 'mmol N/m2/s'; 
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...nitrous oxide plots (decomp) only...')
				vars   = {'N2O','N2O_AO1','N2O_SIDEN','N2O_SODEN','N2O_ATM','N2O_DECOMP'};
				units  = 'mmol N/m2/s'; 
			elseif strcmp(obj.budget.varname,'O2')
				disp('...oxygen plots only...');
				vars   = {'O2'};
				units  = 'mmol O2/m2/s';
			end
	
			% Process inputs
			if isempty(A.time)
				A.time = 1:1:obj.region.nt;
			end

			% First generate uniform color bars for each axis	
			% Go through each variables
			for i = 1:length(vars)
				
				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	
					
					% Gather terms
					dcdt = obj.romsData.(vars{i}).intdcdt(:,:,tt);
					adv  = obj.romsData.(vars{i}).intadx(:,:,tt) + ...
					       obj.romsData.(vars{i}).intady(:,:,tt) + ...
					       obj.romsData.(vars{i}).intadz(:,:,tt);
					dfz  = obj.romsData.(vars{i}).intdfz(:,:,tt);
					sms  = obj.romsData.(vars{i}).intsms(:,:,tt);
					fg   = obj.romsData.(vars{i}).intfg(:,:,tt);
					net  = obj.romsData.(vars{i}).intnet(:,:,tt);
					
					% Blank land
					dcdt = dcdt .* obj.region.mask_rho;
					adv  = adv  .* obj.region.mask_rho;
					dfz  = dfz  .* obj.region.mask_rho;
					sms  = sms  .* obj.region.mask_rho;
					fg   = fg   .* obj.region.mask_rho;
					net  = net  .* obj.region.mask_rho;

					% Get colobars
					cbar_dcdt = prclims(dcdt,'prc',A.prc);
					cbar_adv  = prclims(adv,'prc',A.prc);
					cbar_dfz  = prclims(dfz,'prc',A.prc);
					cbar_sms  = prclims(sms,'prc',A.prc);
					cbar_fg   = prclims(fg,'prc',A.prc);
					cbar_net  = prclims(net,'prc',A.prc);
					if t == 1
						cbar.dcdt = cbar_dcdt;
						cbar.adv  = cbar_adv;
						cbar.dfz  = cbar_dfz;
						cbar.sms  = cbar_sms;
						cbar.fg   = cbar_fg;
						cbar.net  = cbar_net;
					end
					
					% Update colorbars?
					if max(cbar_dcdt) > max(cbar.dcdt)
						cbar.dcdt = cbar_dcdt;
					end
					if max(cbar_adv) > max(cbar.adv)
						cbar.adv = cbar_adv;
					end
					if max(cbar_dfz) > max(cbar.dfz);
						cbar.dfz = cbar_dfz;
					end
					if max(cbar_sms) > max(cbar.sms);
						cbar.sms = cbar_sms;
					end
					if max(cbar_fg)  > max(cbar.fg);
						cbar.fg  = cbar_fg;
					end
					if max(cbar_net) > max(cbar.net)
						cbar.net = cbar_net;
					end
				end

				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	

					% Gather terms
					dcdt = obj.romsData.(vars{i}).intdcdt(:,:,tt);
					adv  = obj.romsData.(vars{i}).intadx(:,:,tt) + ...
					       obj.romsData.(vars{i}).intady(:,:,tt) + ...
					       obj.romsData.(vars{i}).intadz(:,:,tt);
					dfz  = obj.romsData.(vars{i}).intdfz(:,:,tt);
					sms  = obj.romsData.(vars{i}).intsms(:,:,tt);
					fg   = obj.romsData.(vars{i}).intfg(:,:,tt);
					net  = obj.romsData.(vars{i}).intnet(:,:,tt);

                                        % Blank land
					dcdt = dcdt .* obj.region.mask_rho;
					adv  = adv  .* obj.region.mask_rho;
					dfz  = dfz  .* obj.region.mask_rho;
					sms  = sms  .* obj.region.mask_rho;
					fg   = fg   .* obj.region.mask_rho;
					net  = net  .* obj.region.mask_rho;
					
					% Start plots
					for j = 1:6
				
						% Initiate map figure
						fig = piofigs('mfig',1);
						if j == 1
							dat    = dcdt;
							titstr = ['Integrated ',vars{i},': change over time'];
							fstr   = ['dcdt'];
							clevs  = cbar.dcdt;
						elseif j == 2
							dat    = adv;
							titstr = ['Integrated ',vars{i},': advection'];
							fstr   = ['adv'];
							clevs  = cbar.adv;
						elseif j == 3
							dat    = dfz;
							titstr = ['Integrated ',vars{i},': diffusion'];
							fstr   = ['diff'];
							clevs  = cbar.dfz;
						elseif j == 4
							dat    = sms;
							titstr = ['Integrated ',vars{i},': sources-minus-sinks'];
							fstr   = ['sms'];
							clevs  = cbar.sms;
						elseif j == 5
							dat    = fg;
							titstr = ['Integrated ',vars{i},': air-sea flux'];
							fstr   = ['fg'];
							clevs  = cbar.fg;
						elseif j == 6
							dat    = net;
							titstr = ['Integrated ',vars{i},': net remainder'];
							fstr   = ['net'];
							clevs  = cbar.net;
						end
						
						% Plot results
						[ax] = map_plot(fig,obj.region.lon_rho,obj.region.lat_rho);
						clevs        	    = clevs(1):(diff(clevs)/31):clevs(2);
						dat(dat<clevs(1))   = clevs(1);
						dat(dat>clevs(end)) = clevs(end);
						[tmp hc]            = m_contourf(obj.region.lon_rho,...
								      obj.region.lat_rho,dat,clevs,'LineStyle','none');
						cb = colorbar;
						title(titstr,'Interpreter', 'none');
						ylabel(cb,units,'FontSize',6)
						caxis([clevs(1) clevs(end)]);
						colormap(gca,cmocean('balance'));
						ax.FontSize = 10;
						cb.FontSize = 10;
					
						% Print	
						if strcmp(obj.info.time_string,'Monthly')
							fname = [vars{i},'_',fstr,'_M',num2str(tt)];
						elseif strcmp(obj.info.time_string,'Daily');
							fname = [vars{i},'_',fstr,'_D',num2str(tt)];
						end
						export_fig('-jpg',[obj.paths.plots.budget,fname]);
						close(fig)	
					end
					
				end
			end
		end % end method plotBudg
		
		%--------------------------------------------------------------------------------
		function plot1D(obj,varargin)
			% ------------------
			% - Plot the relevant rates versus depth at a lat,lon point
			%
			% - Inputs:
			%   time - time record to plot (default--> all)
			%   lon  - longitude point (can be a vector)
			%   lat  - latitude point (can be a vector)
			%
			% - Example:
			%   plot1D(obj,'time',12,'lon',200,'lat',0)
			% ------------------

                        % Defaults for optional  arguments
                        A.time  = [];
                        A.lon   = [];
			A.lat   = [];
			A.zlim  = [];
                        A       = parse_pv_pairs(A,varargin); % parse method arguments to A

			% Process inputs
			if isempty(A.time)
				A.time = 1:obj.region.nt;
			end 	
			if isempty(A.lon) | isempty(A.lat)
				disp('No lat/lon point selected, killing');
				return
			end
			
			% Get index of nearest lat/lon grid cell
			lon_idx = [];
			lat_idx = [];
			for i = 1:length(A.lon)
				all_dist   = distance(A.lat(i),A.lon(i),obj.region.lat_rho,obj.region.lon_rho);
				[lon_idx(i),lat_idx(i)]  = find(all_dist == min(all_dist(:)));
			end

			% Get colormix
			clrs = colormix(length(A.time),'w');
			for i = 1:length(A.time)
				if clrs(i,:) == [0 0 0];
					clrs(i,:) = [0.3 0.3 0.3];
				end
			end

			% Process time
			if obj.region.nt == 12
				tstr = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
			elseif obj.region.nt == 12 & length(A.time) < 12
				ttstr = [];
				for i = 1:length(A.time)
					if i == 1
						ttstr = [tstr{A.time(i)}];
					else
						ttstr = [ttstr,', ',tstr{A.time(i)}];
					end
				end
			end

 			% Get rates
			if strcmp(obj.budget.varname,'N')	
				disp('...full nitrogen cycle integration...')
				return
			elseif strcmp(obj.budget.varname,'NO3')
				disp('...nitrate integration only...')
				vars = {'NITROX','DENITRIF1','photo_NO3'};
				units = 'mmol N m^{-2} s^{-1}';
			elseif strcmp(obj.budget.varname,'NO2')
				disp('...nitrite integration only...')
				vars  = {'NITROX','DENITRIF1','DENITRIF2','ANAMMOX','photo_NO2','NO2AMMOX'};
				units = 'mmol N m^{-2} s^{-1}';
			elseif strcmp(obj.budget.varname,'N2O')
				disp('...nitrous oxide integration only...')
				vars  = {'DENITRIF2','DENITRIF3','N2OAMMOX'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			elseif strcmp(obj.budget.varname,'N2O_decomp')
				disp('...plotting N2O decomp fluxes...');
                        	vars = {'DENITRIF2','DENITRIF3','N2OAMMOX','N2OSODEN_CONS','N2OAO1_CONS', ...
                                        'N2OATM_CONS', 'N2OSIDEN_CONS','DENITRIF3_DECOMP'};
				units = 'mmol N_2O m^{-2} s^{-1}'; 
			end

			% Plot 1D rates
			for i = 1:length(A.lon)
				for j = 1:length(vars)
					for k = 1:length(A.time)
				
						% Gather data
						if strcmp(obj.budget.varname,'NO2') & strcmp(vars{j},'NO2AMMOX')
							tmpdata = obj.budget.rates.AMMOX - 2*obj.budget.rates.N2OAMMOX;
							tmpdata = squeeze(tmpdata(lon_idx(i),lat_idx(i),:,A.time(k)));
						else
							tmpdata = squeeze(obj.budget.rates.(vars{j})(lon_idx(i),lat_idx(i),:,A.time(k)));
						end

						% Gather depth
						depth = -squeeze(obj.region.z_r(lon_idx(i),lat_idx(i),:,A.time(k)));

						% Make plot(s)
						if k == 1
							fig = piofigs('lfig',1.5);
							set(gca,'YDir','Reverse');
							xlabel(units);
							ylabel('Depth (m)');
							hold on; grid on
						end
						if length(A.time) == 1
							title({[vars{j},' vs. Depth: ',tstr{A.time(k)},' Average'],...
						       	       ['Lon=',num2str(A.lon(i)),', Lat=',num2str(A.lat(i))]});
							legend(tstr{A.time(k)});
						end
						plot(tmpdata,depth,'LineWidth',1.5,'Color',clrs(k,:));
						hold on;
						if ~isempty(A.zlim) & length(A.zlim)==1
							ylim([0 A.zlim]);
						elseif ~isempty(A.zlim) & length(A.zlim)==2
							ylim([A.zlim]);
						end
						if length(A.time) < obj.region.nt
							fname = [vars{j},'_vs_z_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_',tstr{A.time(k)}];
							print('-djpeg',[obj.paths.plots.budget,fname])
							close all
						end
					end
					
					% Add averages if many months plotted
					if length(A.time) == obj.region.nt
						% Gather data
						if strcmp(obj.budget.varname,'NO2') & strcmp(vars{j},'NO2AMMOX')
							tmpdata = obj.budget.rates.AMMOX - 2*obj.budget.rates.N2OAMMOX;
							tmpdata = squeeze(tmpdata(lon_idx(i),lat_idx(i),:,:));
						else
							tmpdata = squeeze(obj.budget.rates.(vars{j})(lon_idx(i),lat_idx(i),:,:));
						end
						tmpdata = nanmean(tmpdata,2);
						plot(tmpdata,depth,'--k','linewidth',4);
						lstr = [tstr];
						lstr(13) = {'Avg'};
						legend(lstr,'Location','Southeast');
						title({[vars{j},' vs. Depth'],...
						       ['Lon=',num2str(A.lon(i)),', Lat=',num2str(A.lat(i))]});
						fname = [vars{j},'_vs_z_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_annual'];
						print('-djpeg',[obj.paths.plots.budget,fname])
						close all
					end
				end

				% Plot locations
				if length(A.time) == obj.region.nt
					pltdata = obj.int.conc.(obj.budget.varname);
				else
					pltdata = obj.int.conc.(obj.budget.varname)(:,:,A.time);
				end
				if length(A.time) == obj.region.nt
					avgdata = nanmean(pltdata,3);
					[fig,ax,cb] = map_ctf(obj.region.lon_rho,obj.region.lat_rho,avgdata,'bal',0,'cmap',cmocean('amp'));
					m_plot(obj.region.lon_rho(lon_idx,lat_idx),obj.region.lat_rho(lon_idx,lat_idx),'k*','markersize',10);
					m_grid('box','on','linestyle','none','backgroundcolor',rgb('LightGray'));			
					fname = [obj.budget.varname,'int_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_annual'];
					ylabel(cb,units(1:end-7))
					title({['Integrated ',obj.budget.varname],['Annual Average']})
					print('-djpeg',[obj.paths.plots.budget,fname]);
					close all
				else
					for t = 1:length(A.time)
						[fig,ax,cb] = map_ctf(obj.region.lon_rho,obj.region.lat_rho,...
							      pltdata(:,:,t),'bal',0,'cmap',cmocean('amp'));
						m_plot(obj.region.lon_rho(lon_idx,lat_idx),obj.region.lat_rho(lon_idx,lat_idx),'k*','markersize',10);
						m_grid('box','on','linestyle','none','backgroundcolor',rgb('LightGray'));			
						fname = [obj.budget.varname,'int_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_',tstr{A.time(t)}];
						ylabel(cb,units(1:end-7))
						title({['Integrated ',obj.budget.varname],[tstr{A.time(t)},' Average']})
						print('-djpeg',[obj.paths.plots.budget,fname]);
						close all
					end
				end
			end
		end % end method plot1D
	
		%--------------------------------------------------------------------------------
		function [obj] = initDiag(obj);	
        		% -------------------
        		% Initialize diagnostic plots: set metadata, save paths, get axis limits 
        		% colorbar bounds, contour levels, titles, and more
        		% 
			% Usage:
			% - obj = initDiag(obj)
			% -------------------
			
			% - check that object is initialized
			try; x = obj.info.var3d; catch; obj = init(obj); end

        	 	% - set paths of data for validation
			% - 0.25 and 1.00 WOA masks
         		obj.paths.diag.woamasks = '/data/project3/data/woa18/masks/WOA_bmsks.mat';	

			% - temperature
			obj.paths.diag.temp.file  = '/data/project3/data/woa18/temperature/0p25/temp_woa18_clim.nc';
			obj.paths.diag.temp.type  = 'nc';
         		obj.paths.diag.temp.var   = 'temp';
			obj.paths.diag.temp.zvar  = 'depth';
			obj.paths.diag.temp.dim   = 'xyzt';
			obj.paths.diag.temp.lon   = 'lon';
			obj.paths.diag.temp.lat   = 'lat';

			% - salinity
			obj.paths.diag.salt.file  = '/data/project3/data/woa18/salinity/0p25/salt_woa18_clim.nc';
			obj.paths.diag.salt.type  = 'nc';
			obj.paths.diag.salt.var   = 'salt';
			obj.paths.diag.salt.zvar  = 'depth';
			obj.paths.diag.salt.dim   = 'xyzt';
			obj.paths.diag.salt.lon   = 'lon';
			obj.paths.diag.salt.lat   = 'lat';

			% - sea surface height
			obj.paths.diag.ssh.file = '/data/project1/demccoy/ROMS/validation/AVISO/monthly_AVISO.mat';
			obj.paths.diag.ssh.type = 'mat';
			obj.paths.diag.ssh.var  = 'adt_month_av';
			obj.paths.diag.ssh.zvar = [];
			obj.paths.diag.ssh.dim  = 'xyt';
			obj.paths.diag.ssh.lon  = 'lon_aviso';
			obj.paths.diag.ssh.lat  = 'lat_aviso';

			% - zonal wind stress
			obj.paths.diag.zws.file = '/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat';
			obj.paths.diag.zws.type = 'mat';
			obj.paths.diag.zws.var  = 'u';
			obj.paths.diag.zws.zvar = [];
			obj.paths.diag.zws.dim  = 'xyt';
			obj.paths.diag.zws.lon  = 'lon';
			obj.paths.diag.zws.lat  = 'lat';

			% - meridional wind stress
			obj.paths.diag.mws.file = '/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat';
			obj.paths.diag.mws.type = 'mat';
			obj.paths.diag.mws.var  = 'v';
			obj.paths.diag.mws.zvar = [];
			obj.paths.diag.mws.dim  = 'xyt';
			obj.paths.diag.mws.lon  = 'lon';
			obj.paths.diag.mws.lat  = 'lat';
			
			% - wind stress
			obj.paths.diag.ws.file = '/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat';
			obj.paths.diag.ws.type = 'mat';
			obj.paths.diag.ws.var  = 'ws';
			obj.paths.diag.ws.zvar = [];
			obj.paths.diag.ws.dim  = 'xyt';
			obj.paths.diag.ws.lon  = 'lon';
			obj.paths.diag.ws.lat  = 'lat';
		
			% - wind stress curl
			obj.paths.diag.wsc.file = '/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat';
			obj.paths.diag.wsc.type = 'mat';
			obj.paths.diag.wsc.var  = 'wsc';
			obj.paths.diag.wsc.zvar = [];
			obj.paths.diag.wsc.dim  = 'xyt';
			obj.paths.diag.wsc.lon  = 'lon';
			obj.paths.diag.wsc.lat  = 'lat';
	
			% - oxygen (only 1.00 available)
			obj.paths.diag.O2.file  = '/data/project3/data/woa18/oxygen/1p0/o2_woa18_clim.nc';
         		obj.paths.diag.O2.type  = 'nc';
			obj.paths.diag.O2.var   = 'o2';
			obj.paths.diag.O2.zvar  = 'depth';
         		obj.paths.diag.O2.dim   = 'xyzt';
         		obj.paths.diag.O2.lon   = 'lon';
         		obj.paths.diag.O2.lat   = 'lat';

			% - nitrate (only 1.00 available)
         		obj.paths.diag.NO3.file = '/data/project3/data/woa18/nitrate/1p0/no3_woa18_clim.nc';
         		obj.paths.diag.NO3.type = 'nc';
         		obj.paths.diag.NO3.var  = 'no3';
			obj.paths.diag.NO3.zvar = 'depth';
         		obj.paths.diag.NO3.dim  = 'xyzt';
         		obj.paths.diag.NO3.lon  = 'lon';
         		obj.paths.diag.NO3.lat  = 'lat';

			% - phosphate (only 1.00 available)
         		obj.paths.diag.PO4.file = '/data/project3/data/woa18/phosphate/1p0/po4_woa18_clim.nc';
         		obj.paths.diag.PO4.type = 'nc';
         		obj.paths.diag.PO4.var  = 'po4';
			obj.paths.diag.PO4.zvar = 'depth';
         		obj.paths.diag.PO4.dim  = 'xyzt';
         		obj.paths.diag.PO4.lon  = 'lon';
			obj.paths.diag.PO4.lat  = 'lat';
        	 
			% - silicate (only 1.00 available)
			obj.paths.diag.SiO3.file = '/data/project3/data/woa18/silicate/1p0/si_woa18_clim.nc';
			obj.paths.diag.SiO3.type = 'nc';
			obj.paths.diag.SiO3.var  = 'si';
			obj.paths.diag.SiO3.zvar = 'depth';
			obj.paths.diag.SiO3.dim  = 'xyzt';
			obj.paths.diag.SiO3.lon  = 'lon';
			obj.paths.diag.SiO3.lat  = 'lat';

			% - N2O (only 1.00 available, reconstructed)
			obj.paths.diag.N2O.file = '/data/project2/yangsi/analysis/n2oInterior/processed/n2opredRF_05-27-2020.nc';
			obj.paths.diag.N2O.type = 'nc';
			obj.paths.diag.N2O.var  = 'n2o';
			obj.paths.diag.N2O.zvar = 'depth';
			obj.paths.diag.N2O.dim  = 'xyzt';
			obj.paths.diag.N2O.lon  = 'lon';
			obj.paths.diag.N2O.lat  = 'lat';

			% - NO2 (only 1.00 available, reconstructed)
			obj.paths.diag.NO2.file = '/data/project2/yangsi/analysis/no2Interior/processed/no2predRF_05-27-2020.nc';
			obj.paths.diag.NO2.type = 'nc';
			obj.paths.diag.NO2.var  = 'no2';
			obj.paths.diag.NO2.zvar = 'depth';
			obj.paths.diag.NO2.dim  = 'xyzt';
			obj.paths.diag.NO2.lon  = 'lon';
			obj.paths.diag.NO2.lat  = 'lat';

			% - NH4 (not available)
			obj.paths.diag.NH4.file = [];
			obj.paths.diag.NH4.type = [];
			obj.paths.diag.NH4.var  = [];
			obj.paths.diag.NH4.zvar = [];
			obj.paths.diag.NH4.dim  = [];
			obj.paths.diag.NH4.lon  = [];
			obj.paths.diag.NH4.lat  = [];

			% - Chl 
            		obj.paths.diag.chl.file = '/data/project1/data/MODIS-Aqua/backup/res0p25/chl_clim_0p25.nc';
            		obj.paths.diag.chl.type = 'nc';
            		obj.paths.diag.chl.dim  = 'xyt';
			obj.paths.diag.chl.var  = 'chl';
            		obj.paths.diag.chl.zvar = [];
			obj.paths.diag.chl.lon  = 'lon';
            		obj.paths.diag.chl.lat  = 'lat';
         		
            		% - NPP-VGPM	
			obj.paths.diag.nppvgpm.file = '/data/project2/yangsi/analysis/NPPcode/std-VGPM/std-VGPMnpp_MODIS_clim2002-2018.nc';
            		obj.paths.diag.nppvgpm.type = 'nc';
            		obj.paths.diag.nppvgpm.var  = 'npp';
			obj.paths.diag.nppvgpm.zvar = [];
            		obj.paths.diag.nppvgpm.dim  = 'xyt';
            		obj.paths.diag.nppvgpm.lon  = 'lon';
            		obj.paths.diag.nppvgpm.lat  = 'lat';
            		
			% - NPP-CBPM
			obj.paths.diag.nppcbpm.file = '/data/project2/yangsi/analysis/NPPcode/CbPM2/CbPM2npp_MODIS_clim2002-2018.nc';
            		obj.paths.diag.nppcbpm.type = 'nc';
            		obj.paths.diag.nppcbpm.var  = 'npp';
			obj.paths.diag.nppcbpm.zvar = [];
            		obj.paths.diag.nppcbpm.dim  = 'xyt';
            		obj.paths.diag.nppcbpm.lon  = 'lon';
            		obj.paths.diag.nppcbpm.lat  = 'lat';

			% - NPP-CAFE
			obj.paths.diag.nppcafe.file = '/data/project1/data/MODIS-Aqua/CAFE_NPP/climatology/nppclim_CAFE_MODIS_9km.nc';
			obj.paths.diag.nppcafe.type = 'nc';
			obj.paths.diag.nppcafe.var  = 'npp';
			obj.paths.diag.nppcafe.zvar = [];
			obj.paths.diag.nppcafe.dim  = 'xyt';
			obj.paths.diag.nppcafe.lon  = 'lon';
			obj.paths.diag.nppcafe.lat  = 'lat';
		
			% - Set pathways for diagnostic plots	
			% - get colormap, determine surface level for WOA, get number of levels
	        	cmap                = load('/data/project1/yangsi/MATLAB/functions/colormaps/cmap.mat');
       			obj.plots.diag.cmap_diff = cmap.cmap_ferretdiff;
         		obj.plots.diag.surfAvg   = 50;
        		obj.plots.diag.surfLevs  = find(obj.grid.woa1p0.depth <= obj.plots.diag.surfAvg);
       			obj.plots.diag.surfNlevs = length(obj.plots.diag.surfLevs);

         		% O2
			obj.plots.diag.O2.roms.var        = 'O2';
			obj.plots.diag.O2.roms.dim        = 'xyzt';
			obj.plots.diag.O2.data.var        = 'o2';
			obj.plots.diag.O2.data.varsource  = 'WOA18';
			obj.plots.diag.O2.units           = 'mMol/m^3';
			obj.plots.diag.O2.surf.absLevs    = [160:20:340];
			obj.plots.diag.O2.surf.absCaxis   = [160,340];
			obj.plots.diag.O2.surf.diffLevs   = [-50:10:50];
			obj.plots.diag.O2.surf.diffCaxis  = [-50,50];
			obj.plots.diag.O2.surf.logplt     = 0;
			
			% PO4
			obj.plots.diag.PO4.roms.var       = 'PO4';
         		obj.plots.diag.PO4.roms.dim       = 'xyzt';
         		obj.plots.diag.PO4.data.var       = 'po4';
         		obj.plots.diag.PO4.data.varsource = 'WOA18';
         		obj.plots.diag.PO4.units          = 'mMol/m^3';
         		obj.plots.diag.PO4.surf.absLevs   = (0:0.2:1.8);
         		obj.plots.diag.PO4.surf.absCaxis  = [0,1.8];
         		obj.plots.diag.PO4.surf.diffLevs  = (-1.0:0.25: 1.0);
         		obj.plots.diag.PO4.surf.diffCaxis = [-1,1];
         		obj.plots.diag.PO4.surf.logplt    = 0;

         		% NO3
         		obj.plots.diag.NO3.roms.var       = 'NO3';
         		obj.plots.diag.NO3.roms.dim       = 'xyzt';
         		obj.plots.diag.NO3.data.var       = 'no3';
         		obj.plots.diag.NO3.data.varsource = 'WOA18';
         		obj.plots.diag.NO3.units          = 'mMol/m^3';
         		obj.plots.diag.NO3.surf.absLevs   = (0:1.0:24);
         		obj.plots.diag.NO3.surf.absCaxis  = [0,24];
         		obj.plots.diag.NO3.surf.diffLevs  = (-12:1:12);
         		obj.plots.diag.NO3.surf.diffCaxis = [-12,12];
         		obj.plots.diag.NO3.surf.logplt    = 0;

			% SiO3
			obj.plots.diag.SiO3.roms.var       = 'SiO3';
			obj.plots.diag.SiO3.roms.dim       = 'xyzt';
			obj.plots.diag.SiO3.data.var       = 'si';
			obj.plots.diag.SiO3.data.varsource = 'WOA18';
			obj.plots.diag.SiO3.units          = 'mMol/m^3';
			obj.plots.diag.SiO3.surf.absLevs   = (0:20:150);
			obj.plots.diag.SiO3.surf.absCaxis  = [0,150];
			obj.plots.diag.SiO3.surf.diffLevs  = (-50:10:50);
			obj.plots.diag.SiO3.surf.diffCaxis = [-50,50];
			obj.plots.diag.SiO3.surf.logplt    = 0;

			% N2O
			obj.plots.diag.N2O.roms.var       = 'N2O';
			obj.plots.diag.N2O.roms.dim       = 'xyzt';
			obj.plots.diag.N2O.data.var       = 'n2o';
			obj.plots.diag.N2O.data.varsource = 'GLODAPv2';
			obj.plots.diag.N2O.units          = 'mMol/m^3';
			obj.plots.diag.N2O.surf.absLevs   = [0:0.005:0.05];
			obj.plots.diag.N2O.surf.absCaxis  = [0,0.05];
			obj.plots.diag.N2O.surf.diffLevs  = [-0.025:0.005:0.025];
			obj.plots.diag.N2O.surf.diffCaxis = [-0.025,0.025];
			obj.plots.diag.N2O.surf.logplt    = 0;
         		
			% NO2
			obj.plots.diag.NO2.roms.var       = 'NO2';
			obj.plots.diag.NO2.roms.dim       = 'xyzt';
			obj.plots.diag.NO2.data.var       = 'no2';
			obj.plots.diag.NO2.data.varsource = 'GLODAPv2';
			obj.plots.diag.NO2.units          = 'mMol/m^3';
			obj.plots.diag.NO2.surf.absLevs   = (0:0.3:3);
			obj.plots.diag.NO2.surf.absCaxis  = [0,3];
			obj.plots.diag.NO2.surf.diffLevs  = (-3:0.6:3);
			obj.plots.diag.NO2.surf.diffCaxis = [-3,3];
			obj.plots.diag.NO2.surf.logplt    = 0;
	
			% NH4
			obj.plots.diag.NH4.roms.var       = 'NH4';
			obj.plots.diag.NH4.roms.dim       = 'xyzt';
			obj.plots.diag.NH4.data.var       = [];
			obj.plots.diag.NH4.data.varsource = [];
			obj.plots.diag.NH4.units          = 'mMol/m^3';
			obj.plots.diag.NH4.surf.absLevs   = (0:0.05:0.2);
			obj.plots.diag.NH4.surf.absCaxis  = [0,0.2];
			obj.plots.diag.NH4.surf.diffLevs  = (-1:0.5:1);
			obj.plots.diag.NH4.surf.diffCaxis = [-1,1];
			obj.plots.diag.NH4.surf.logplt    = 0;	

			% N*
			obj.plots.diag.nstar.roms.var       = 'nstar';
         		obj.plots.diag.nstar.roms.dim       = 'xyzt';
		        obj.plots.diag.nstar.data.var       = 'nstar';
         		obj.plots.diag.nstar.data.varsource = 'WOA18';
         		obj.plots.diag.nstar.units          = 'mMol/m^3';
         		obj.plots.diag.nstar.surf.absLevs   = (-10:0.5:10);
         		obj.plots.diag.nstar.surf.absCaxis  = [-10,10];
         		obj.plots.diag.nstar.surf.diffLevs  = (-10:0.5: 10);
         		obj.plots.diag.nstar.surf.diffCaxis = [-10,10];
         		obj.plots.diag.nstar.surf.logplt    = 0;

        		% chlA
         		obj.plots.diag.chl.roms.var       = 'TOT_CHL';
		        obj.plots.diag.chl.roms.dim       = 'xyt';
		        obj.plots.diag.chl.data.var       = 'chl';
		        obj.plots.diag.chl.data.varsource = 'Modis-Aqua';
         		obj.plots.diag.chl.units          = 'mg/m^3';
        		obj.plots.diag.chl.surf.absLevs   = [0.01 0.02 0.05  0.1 0.2 0.4  0.8 1.5 3.0 6.0 10.0];
         		obj.plots.diag.chl.surf.absCaxis  = [0.01,10];
         		obj.plots.diag.chl.surf.diffLevs  = [0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0];
         		obj.plots.diag.chl.surf.diffLevs  = [-fliplr(obj.plots.diag.chl.surf.diffLevs) obj.plots.diag.chl.surf.diffLevs];
         		obj.plots.diag.chl.surf.diffCaxis = [-10,10];
        		obj.plots.diag.chl.surf.logplt    = 1;

		        % npp vgpm
        		obj.plots.diag.nppvgpm.roms.var       = 'NPP';
       			obj.plots.diag.nppvgpm.roms.dim       = 'xyt';
       			obj.plots.diag.nppvgpm.data.var       = 'nppvgpm';
       			obj.plots.diag.nppvgpm.data.varsource = 'Modis-Aqua';
       			obj.plots.diag.nppvgpm.units          = 'mgC/m^{2}/d';
       			obj.plots.diag.nppvgpm.surf.absLevs   = (0:50:1500); 
       			obj.plots.diag.nppvgpm.surf.absCaxis  = [0,1500];
       			obj.plots.diag.nppvgpm.surf.diffLevs  = (-1000:50:1000);
       			obj.plots.diag.nppvgpm.surf.diffCaxis = [-1000,1000];
        		obj.plots.diag.nppvgpm.surf.logplt    = 0;

      		 	% npp cbpm
       			obj.plots.diag.nppcbpm.roms.var       = 'NPP';
       			obj.plots.diag.nppcbpm.roms.dim       = 'xyt';
       			obj.plots.diag.nppcbpm.data.var       = 'nppcbpm';
       			obj.plots.diag.nppcbpm.data.varsource = 'Modis-Aqua';
       			obj.plots.diag.nppcbpm.units          = 'mgC/m^{2}/d';
       			obj.plots.diag.nppcbpm.surf.absLevs   = (0:50:1500) ;
       			obj.plots.diag.nppcbpm.surf.absCaxis  = [0,1500];
       			obj.plots.diag.nppcbpm.surf.diffLevs  = (-1000:50:1000);
       			obj.plots.diag.nppcbpm.surf.diffCaxis = [-1000,1000];
       			obj.plots.diag.nppcbpm.surf.logplt    = 0;

			% npp cafe
			obj.plots.diag.nppcafe.roms.var       = 'NPP';
			obj.plots.diag.nppcafe.roms.dim       = 'xyt';
			obj.plots.diag.nppcafe.data.var       = 'nppcafe';
			obj.plots.diag.nppcafe.data.varsource = 'Modis-Aqua';
			obj.plots.diag.nppcafe.units          = 'mgC/m^{2}/d';
			obj.plots.diag.nppcafe.surf.absLevs   = (0:50:1500);
			obj.plots.diag.nppcafe.surf.absCaxis  = [0,1500];
			obj.plots.diag.nppcafe.surf.diffLevs  = (-1000:50:1000);
			obj.plots.diag.nppcafe.surf.diffCaxis = [-1000,1000];
			obj.plots.diag.nppcafe.surf.logplt    = 0;

			% temperature
			obj.plots.diag.temp.roms.var       = 'temp';
			obj.plots.diag.temp.roms.dim       = 'xyzt';
			obj.plots.diag.temp.data.var       = 'temp';
			obj.plots.diag.temp.data.varsource = 'WOA18';
			obj.plots.diag.temp.units          = '^oC';
			obj.plots.diag.temp.surf.absLevs   = [0:2:30];
			obj.plots.diag.temp.surf.absCaxis  = [0,30];
			obj.plots.diag.temp.surf.diffLevs  = [-3:1:3];
			obj.plots.diag.temp.surf.diffCaxis = [-3,3];
			obj.plots.diag.temp.surf.logplt    = 0;

			% salinity
			obj.plots.diag.salt.roms.var       = 'salt';
			obj.plots.diag.salt.roms.dim       = 'xyzt';
			obj.plots.diag.salt.data.var       = 'salt';
			obj.plots.diag.salt.data.varsource = 'WOA18';
			obj.plots.diag.salt.units          = 'PSU';
			obj.plots.diag.salt.surf.absLevs   = [32:0.25:38];
			obj.plots.diag.salt.surf.absCaxis  = [32,38];
			obj.plots.diag.salt.surf.diffLevs  = [-2:1:2];
			obj.plots.diag.salt.surf.diffCaxis = [-2,2];
			obj.plots.diag.salt.surf.logplt    = 0;

			% sea surface height
			obj.plots.diag.ssh.roms.var       = 'zeta';
			obj.plots.diag.ssh.roms.dim       = 'xyt';
			obj.plots.diag.ssh.data.var       = 'adt_month_av';
			obj.plots.diag.ssh.data.varsource = 'AVISO';
			obj.plots.diag.ssh.units          = 'meters';
			obj.plots.diag.ssh.surf.absLevs   = [-1.5:0.1:1.5];
			obj.plots.diag.ssh.surf.absCaxis  = [-1.5,1.5];
			obj.plots.diag.ssh.surf.diffLevs  = [-0.5:0.1:0.5];
			obj.plots.diag.ssh.surf.diffCaxis = [-0.5,0.5];
			obj.plots.diag.ssh.surf.logplt    = 0;

			% zonal wind stress
			obj.plots.diag.zws.roms.var       = 'sustr';
			obj.plots.diag.zws.roms.dim       = 'xyt';
			obj.plots.diag.zws.data.var 	  = 'u';
			obj.plots.diag.zws.data.varsource = 'SCOW 2010';
			obj.plots.diag.zws.units          = 'N m^{-2}';
			obj.plots.diag.zws.surf.absLevs   = [0:0.01:0.15];
			obj.plots.diag.zws.surf.absCaxis  = [0,0.15];
			obj.plots.diag.zws.surf.diffLevs  = [-0.05:0.01:0.05];
			obj.plots.diag.zws.surf.diffCaxis = [-0.05,0.05];
			obj.plots.diag.zws.surf.logplt    = 0;

			% meridional wind stress
			obj.plots.diag.mws.roms.var       = 'svstr';
			obj.plots.diag.mws.roms.dim       = 'xyt';
			obj.plots.diag.mws.data.var       = 'v';
			obj.plots.diag.mws.data.varsource = 'SCOW 2010';
			obj.plots.diag.mws.units          = 'N m^{-2}';
			obj.plots.diag.mws.surf.absLevs   = [0:0.01:0.15];
			obj.plots.diag.mws.surf.absCaxis  = [0,0.15];
			obj.plots.diag.mws.surf.diffLevs  = [-0.05:0.01:0.05];
			obj.plots.diag.mws.surf.diffCaxis = [-0.05,0.05];
			obj.plots.diag.mws.surf.logplt    = 0;

			% wind stress
			obj.plots.diag.ws.roms.var       = 'ws';
			obj.plots.diag.ws.roms.dim       = 'xyt';
			obj.plots.diag.ws.data.var       = 'ws';
			obj.plots.diag.ws.data.varsource = 'SCOW 2010';
			obj.plots.diag.ws.units          = 'N m^{-2}';
			obj.plots.diag.ws.surf.absLevs   = [0:0.01:0.15];
			obj.plots.diag.ws.surf.absCaxis  = [0,0.15];
			obj.plots.diag.ws.surf.diffLevs  = [-0.05:0.01:0.05];
			obj.plots.diag.ws.surf.diffCaxis = [-0.05,0.05];
			obj.plots.diag.ws.surf.logplt    = 0;
	
			% wind stress curl
			obj.plots.diag.wsc.roms.var       = 'wsc';
			obj.plots.diag.wsc.roms.dim       = 'xyt';
			obj.plots.diag.wsc.data.var       = 'wsc';
			obj.plots.diag.wsc.data.varsource = 'SCOW 2010';
			obj.plots.diag.wsc.units          = 'N m^{-2}';
			obj.plots.diag.wsc.surf.absLevs   = [-2e-7:5e-8:2e-7];
			obj.plots.diag.wsc.surf.absCaxis  = [-2e-7,2e-7];
			obj.plots.diag.wsc.surf.diffLevs  = [-2e-7:5e-8:2e-7];
			obj.plots.diag.wsc.surf.diffCaxis = [-2e-7,2e-7];
			obj.plots.diag.wsc.surf.logplt    = 0;

			% - MLD from ifremer density threshold
			obj.plots.diag.mld_ifremer.roms.var       = 'mld';
			obj.plots.diag.mld_ifremer.roms.dim       = 'xyt';
			obj.plots.diag.mld_ifremer.data.var       = 'mld_ifremer';
			obj.plots.diag.mld_ifremer.data.varsource = 'WOCE/NODC/ARGO';
			obj.plots.diag.mld_ifremer.units          = 'meters';
			obj.plots.diag.mld_ifremer.surf.absLevs   = [0:10:150];
			obj.plots.diag.mld_ifremer.surf.absCaxis  = [0,150];
			obj.plots.diag.mld_ifremer.surf.diffLevs  = [-50:10:50];
			obj.plots.diag.mld_ifremer.surf.diffCaxis = [-50,50];
			obj.plots.diag.mld_ifremer.surf.logplt    = 0;

			% - MLD from argo climatology
			obj.plots.diag.mld_argo.roms.var       = 'mld';
			obj.plots.diag.mld_argo.roms.dim       = 'xyt';
			obj.plots.diag.mld_argo.data.var       = 'mld_argo';
			obj.plots.diag.mld_argo.data.varsource = 'ARGO';
			obj.plots.diag.mld_argo.units          = 'meters';
			obj.plots.diag.mld_argo.surf.absLevs   = [0:10:150];
			obj.plots.diag.mld_argo.surf.absCaxis  = [0,150];
			obj.plots.diag.mld_argo.surf.diffLevs  = [-50:10:50];
			obj.plots.diag.mld_argo.surf.diffCaxis = [-50,50];
			obj.plots.diag.mld_argo.surf.logplt    = 0;
		
			% - u velocity
			obj.plots.diag.u.roms.var       = 'u';
			obj.plots.diag.u.roms.dim	= 'xyzt';
			obj.plots.diag.u.data.var       = 'u';
			obj.plots.diag.u.data.varsource = 'Cravatte (2017)';
			obj.plots.diag.u.units          = 'm/s';
			obj.plots.diag.u.surf.absLevs   = [-1:0.2:1];
			obj.plots.diag.u.surf.absCaxis  = [-1,1];
			obj.plots.diag.u.surf.diffLevs  = [-0.5:0.1:0.5];
			obj.plots.diag.u.surf.diffCaxis = [-0.5,0.5];
			obj.plots.diag.u.surf.logplt    = 0;
		end % end methods initPlots

		%--------------------------------------------------------------------------------
		function obj = loadData(obj,varargin)
        		% -------------------
			% Method to Load ROMS data. Load either raw data,
			% or data interpolated to standard WOCE depths (z_avg)
			%
			% Usage:
			% - obj = loadData(obj,varargin)
			%
			% Inputs (varargin):
                        % - (empty) = select ROMS variables			
 			% - vars    = ROMS variable to load, as a cell array
			% - type    = 'z_avg' or 'raw'
			% 
			% Examples:
			% - obj = loadData(obj)
			% - obj = loadData(obj,'vars',{'temp','salt'},'type','z_avg');
			% -------------------
			
			% - check that object is initialized
			try; x = obj.info.var3d; catch; obj = init(obj); end
			
			% - defaults for optional  arguments
       			A.vars = []; % input variable
			A.type = []; % input variable type
        		A      = parse_pv_pairs(A,varargin); % parse method arguments to A
		
			% - check that region has been defined
			if isempty(obj.region)
				obj = defineRegion(obj);
			end
	
			% - check inputs
			if ~isempty(A.type) % type
				if ~strcmp(A.type,'raw') & ~strcmp(A.type,'z_avg');
					disp(' ');
					disp('type must be: raw or z_avg');
					disp(' '); return
				end
			end
			if ~isempty(A.vars) % vars
				for i = 1:length(A.vars)
					if ~iscell(A.vars(i));
						disp(' ');
						disp('vars must be a cell array');
						disp(' '); 
						return
					end
					if ~ismember(A.vars{i},obj.info.var3d) & ~ismember(A.vars{i},obj.info.var2d)
						disp(' ');
						disp([A.vars{i},' does not exist']);
						disp(' '); return
					end
				end
			end			

			% - if type is defined, select correct choices
			% - if not, then define type
			if strcmp(A.type,'raw');
				q0 = 1;
			elseif strcmp(A.type,'z_avg');
				q0 = 2;
			elseif isempty(A.type)
                        	q0 = input('---------------------------\nraw(1) or z_avg(2)\n---------------------------\n>> ');
			end

			% - load variables 
			if ~isempty(A.vars) % - user provided inputs, load variables
				for i = 1:length(A.vars)
					idx = find(ismember(obj.info.var2d,A.vars{i})==1);
					idx = obj.info.idx2d(idx);
					typ = 2;
					if isempty(idx)
						idx = find(ismember(obj.info.var3d,A.vars{i})==1);
						idx = obj.info.idx3d(idx);
						typ = 3;
					end
					disp('-------------------------');
					disp(['Loading ',A.vars{i},' data...']);
					if q0 == 1
						tmpdata = ncread(obj.paths.avg,A.vars{i});
					elseif q0 == 2
						tmpdata = ncread(obj.paths.zavg,A.vars{i});
					end
					% Reduce to region
					if typ == 2
						tmpdata 		      = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
                                                                 			obj.region.lat_lim(1):obj.region.lat_lim(2),:);
						tmpdata                       = tmpdata .* obj.region.mask_rho;
						obj.romsData.(A.vars{i}).data = tmpdata;
					elseif typ == 3
						tmpdata 		      = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
                                                                 		       	obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
						if q0 == 1 
							tmpdata = tmpdata .* obj.region.mask_rho3d;
						elseif q0 == 2
							tmpdata = tmpdata .* obj.region.mask_rhoz3d;
						end
						obj.romsData.(A.vars{i}).data = tmpdata;
					end
					disp(['...done!']);
					disp('-------------------------');
					obj.romsData.(A.vars{i}).name  = obj.info.Variables(idx).Attributes(1).Value;
					obj.romsData.(A.vars{i}).units = obj.info.Variables(idx).Attributes(2).Value;
				end
			else % - no input, show available fields
				if q0 == 2
					q1 = 1; % only 3D available in z_avg file
				else
                        		q1 = input('---------------------------\n3D(1) or 2D(2) Variable:\n---------------------------\n>> ');
				end
				if q1<1 | q1 >2
					disp('Bad choice, only 1 or 2  allowed');
					return
				elseif q1 == 1
					disp('-------------------------');
					disp('3-D Variables:');
					disp('-------------------------');
					for i = 1:length(obj.info.var3d);
						disp([num2str(i),'==',obj.info.var3d{i}]);	
					end
                        		q2 = input('---------------------------\nChoose from above:\n---------------------------\n>> ');
					if ~ismember(q2,[1:length(obj.info.var3d)]);
						disp('Bad choice, number out of range');
						return
					else
						for j = 1:length(q2)
							disp('-------------------------');
							disp(['Loading ',obj.info.var3d{q2(j)},' data...']);
							if q0 == 1
								tmpdata = ncread(obj.paths.avg,obj.info.var3d{q2(j)});

							elseif q0 == 2
								tmpdata = ncread(obj.paths.zavg,obj.info.var3d{q2(j)});
							end
							% Reduce to region
							tmpdata = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
                                                                         obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
							if q0 == 1 
								tmpdata = tmpdata .* obj.region.mask_rho3d;
							elseif q0 == 2
								tmpdata = tmpdata .* obj.region.mask_rhoz3d;
							end
							obj.romsData.(obj.info.var3d{q2(j)}).data = tmpdata;
							disp(['...done!']);
							disp('-------------------------');
							obj.romsData.(obj.info.var3d{q2(j)}).name  = obj.info.name3d{q2(j)};
							obj.romsData.(obj.info.var3d{q2(j)}).units = obj.info.unit3d{q2(j)};
						end
					end
				elseif q1 == 2
					disp('-------------------------');
					disp('2-D Variables:');
					disp('-------------------------');
					for i = 1:length(obj.info.var2d);
						disp([num2str(i),'==',obj.info.var2d{i}]);	
					end
                        		q2 = input('---------------------------\nChoose from above:\n---------------------------\n>> ');
					if ~ismember(q2,[1:length(obj.info.var2d)]);
						disp('Bad choice');
						return
					else
						for j = 1:length(q2)
							disp('-------------------------');
							disp(['Loading ',obj.info.var2d{q2(j)},' data...']);
	
							if q0 == 1
								tmpdata = ncread(obj.paths.avg,obj.info.var2d{q2(j)});
							elseif q0 == 2
								tmpdata = ncread(obj.paths.zavg,obj.info.var2d{q2(j)});
							end
							% Reduce to region
							tmpata = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
                                                                         obj.region.lat_lim(1):obj.region.lat_lim(2),:);
							tmpdata = tmpdata .* obj.region.mask_rho3d;
							obj.romsData.(obj.info.var2d{q2(j)}).data = tmpdata;
							disp(['...done!']);
							disp('-------------------------');
							obj.romsData.(obj.info.var2d{q2(j)}).name  = obj.info.name2d{q2(j)};
							obj.romsData.(obj.info.var2d{q2(j)}).units = obj.info.unit2d{q2(j)};					
						end
					end
				end	
			end
		end % end methods loadData

		%--------------------------------------------------------------------------------
		function [fig,ax,cb] = plotData(obj,varargin)
        		% -------------------
			% Plot ROMS data maps based on inputs.
			% If choosing both 2D and 3D, note that depth inputs will only be applied
			% to 3D data, as expected.
			%
			% Usage:
			% - [fig,ax,cb] = plotData(obj,varargin)
			% 
			% Inputs (varargin):
			% - (empty) = plot loaded data or select ROMS variables			
			% - vars    = ROMS variable to plot, as a string or cell array
			% - depths  = Depths to plot (see obj.grid.z_avg_dep for list of available depths)
			%             Will use nearest depth when necessary
			% - levels  = Levels to plot (if using terrain-following data)
			%             Cannot be used with 'depths' input
			% - time    = cell array of times to plot
			%             'DJF','MAM','JJA','SON','ANN',or 3-letter month name (Dec, Jan, etc)
			%	      i.e. {'DJF','MAM','JJA','SON'} for seasonal plots
			% - plttype = none (default), tmpfigs, or define file type ('png' or 'pdf');
			%
			% Outputs:
			% - fig = figure handle(s)
			% - ax  = axes handle(s)
			% - cb  = colorbar handle(s)
			%
			% Examples:
			% - obj          = plotData(obj)
			% - [fig,ax,cb]  = plotData(obj,'vars',{'temp','salt'},'depths',[100 500 1000],'time',{'DJF'},'plttype','png');
			% -------------------

			close all

			% - defaults for optional  arguments
       			A.vars     = []; % - input variable
			A.plttype  = []; % - plot type, tmpfig by default
			A.depths   = []; % - depths to plot
			A.levels   = []; % - levels to plot
			A.time     = []; % - time averaging or month to plot
			A.caxis    = []; % - optional caxis input (not recommended)
        		A          = parse_pv_pairs(A,varargin); % - parse method arguments to A

			% - month string
			mstr = {'Jan','Feb','Mar','Apr','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

			% - depth and level indexes
			lidx = []; didx = [];

			% - check that object is initialized
			try; x = obj.info.var3d; catch; obj = init(obj); end

			% - check inputs
			if ~isempty(A.vars) % vars
                                for i = 1:length(A.vars)
					if ~iscell(A.vars(i));
						disp(' ');
						disp('vars must be a cell array');
						disp(' '); return
					end
                                        if ~ismember(A.vars(i),obj.info.var3d) & ~ismember(A.vars(i),obj.info.var2d)
                                                disp(' ');
                                                disp([A.vars{i},'does not exist']);
                                                disp(' '); return
                                        end
                                end
			end
			if ~isempty(A.plttype) % plttype
				if ~strcmp(A.plttype,'tmpfigs') & ~strcmp(A.plttype,'png') & ~strcmp(A.plttype,'pdf') 
					disp(' ');
					disp('plttype must be: tmpfigs, png, pdf');
					disp(' '); return
				end
			end
			if ~isempty(A.depths) & isempty(A.levels) % depths
				for i = 1:length(A.depths)
					if A.depths(i) < 0 & A.depths(i) > max(obj.grid.z_avg_dep)
						disp(' ');
						disp(['depths must be positive and less than ',num2str(max(obj.grid.z_avg_dep))]);
						disp(' '); return
					end
				end
			elseif isempty(A.depths) & ~isempty(A.levels) % levels
				for i = 1:length(A.levels)
					if A.levels(i) <0 & A.levels(i) > obj.region.nz;
						disp(' ');
						disp(['levels must be positive and less than ',num2str(obj.region.nz)]);
						disp(' '); return
					end
				end
			elseif ~isempty(A.depths) & ~isempty(A.levels)
				disp(' ');
				disp('Cannot process both depths and levels, type help obj.plotData');
				disp(' '); return
			end
			if ~isempty(A.time) & obj.region.nt == 12 % time
				for i = 1:length(A.time)
					if ismember(A.time(i),{'DJF','MAM','JJA','SON','ANN'}) | ismember(A.time(i),mstr);
					else
						disp(' ');
						disp('time must be DJF, MAM, JJA, SON, ANN, or a 3-letter month cell string');
						disp(' '); return
					end
				end
			elseif ~isempty(A.time)
				for i = 1:length(A.time)
					if A.time(i) > obj.region.nt
						disp(' ');
						disp('time input is greater than number of records');
						disp(' '); return
					end
				end
			end

			% - process depth input
			if ~isempty(A.depths)
				for i = 1:length(A.depths)
					diffd       = abs(A.depths(i) - obj.grid.z_avg_dep);
					ind         = find(diffd == min(diffd));
					didx(i)     = ind;
					dstr{i}     = [num2str(obj.grid.z_avg_dep(ind)),'m'];
				end
			end

			% - process levels input
			if ~isempty(A.levels)
				lidx = A.levels;
				for i = 1:length(lidx)
					lstr{i} = ['Level ',num2str(lidx(i))];
				end
			end
			
			% - if depth/levels empty, choose
			if isempty(A.depths) & isempty(A.levels)
				q0 = input('---------------------------\nSurface(0), Depths(1), or Levels(2):\n---------------------------\n>> ');
				if ~ismember(q0,[0 1 2]);
					disp(' ');
					disp('Value must be 0, 1, or 2');
					disp(' '); return
				elseif q0 == 0
					didx = [];
					lidx = [];
				elseif q0 == 1
					lidx = [];
			        	disp('-------------------------');
			     		disp('z_avg Depths:');
					for i = 1:length(obj.grid.z_avg_dep)
						if i < 10
							disp([num2str(i),'   == ',num2str(obj.grid.z_avg_dep(i)),'m']);
						elseif i < 100
							disp([num2str(i),'  == ', num2str(obj.grid.z_avg_dep(i)),'m']);
						else
							disp([num2str(i),' == ',  num2str(obj.grid.z_avg_dep(i)),'m']);
						end
					end
					didx = input('---------------------------\nChoose depth(s):\n---------------------------\n>> ');
                               		for i = 1:length(didx)
						if ~ismember(didx(i),[1:length(obj.grid.z_avg_dep)]);
                                        		disp(' ');
                                       	       	        disp('Bad Choice, number out of range');
                                               	        disp(' '); return		
						else
							dstr{i} = [num2str(obj.grid.z_avg_dep(didx(i))),'m'];
						end
					end	
				elseif q0 == 2
					didx = [];
			        	disp('-------------------------');
			     		disp(['Number of levels = ',num2str(obj.region.nz)]);
					lidx = input('---------------------------\nChoose level(s):\n---------------------------\n>> ');
					for i = 1:length(lidx)
						if ~ismember(tidx(i),[1:length(obj.region.nz)]);
                                       			disp(' ');
                                      		       	disp('Bad Choice, number out of range');
                                	       	     	disp(' '); return		
						else
							lstr{i} = ['Level ',num2str(lidx(i))];
						end
					end
				end
			end

			% - process variable input
			if ~isempty(A.vars) % - user provided input
				for i = 1:length(A.vars)
					if isempty(didx) & isempty(lidx)
						obj = loadData(obj,'vars',A.vars(i),'type','raw');
					elseif isempty(didx)
						obj = loadData(obj,'vars',A.vars(i),'type','raw');
					elseif isempty(lidx)
						obj = loadData(obj,'vars',A.vars(i),'type','z_avg');
					end
				end
			else % - no user input
				if isempty(didx)
					obj    = loadData(obj,'type','raw');
				else
					obj    = loadData(obj,'type','z_avg');
				end
				A.vars = fieldnames(obj.romsData);
			end

			% - process time input
			if isempty(A.time) & obj.region.nt == 12
			        disp('-------------------------');
				disp('DJF, MAM, JJA, SON, ANN, or ');
				disp('3-letter month string cell array');
				A.time = input('---------------------------\nChoose time(s):\n---------------------------\n>> ');
				for i = 1:length(A.time)
	                                if ismember(A.time(i),{'DJF','MAM','JJA','SON','ANN'}) | ismember(A.time(i),mstr);
                                	else
                                		disp(' ');
                                		disp('time must be DJF, MAM, JJA, SON, ANN, or a 3-letter month cell string');
                                		disp(' '); return
                			end
					if strcmp(A.time(i),'DJF');
						tidx{i} = [1 2 12];
						tstr{i} = 'DJF';
					elseif strcmp(A.time(i),'MAM');
						tidx{i} = [3:5];
						tstr{i} = 'MAM';
					elseif strcmp(A.time(i),'JJA');
						tidx{i} = [6:8];
						tstr{i} = 'JJA';
					elseif strcmp(A.time(i),'SON');
						tidx{i} = [9:11];
						tstr{i} = 'SON';
					elseif strcmp(A.time(i),'ANN');
						tidx{i} = [1:12];
						tstr{i} = 'Annual';
					elseif ismember(A.time(i),mstr)
						tidx{i} = find(strcmp(mstr,A.time(i))==1);
						tstr{i} = mstr{i};
					end
				end
			elseif isempty(A.time) 
				for i = 1:length(A.time)
					tidx{i} = A.time(i);
					tstr{i} = ['Record ',num2str(A.time(i))];
				end
			elseif ~isempty(A.time)
				for i = 1:length(A.time)
					if strcmp(A.time(i),'DJF');
						tidx{i} = [1 2 12];
						tstr{i} = 'DJF';
					elseif strcmp(A.time(i),'MAM');
						tidx{i} = [3:5];
						tstr{i} = 'MAM';
					elseif strcmp(A.time(i),'JJA');
						tidx{i} = [6:8];
						tstr{i} = 'JJA';
					elseif strcmp(A.time(i),'SON');
						tidx{i} = [9:11];
						tstr{i} = 'SON';
					elseif strcmp(A.time(i),'ANN');
						tidx{i} = [1:12];
						tstr{i} = 'Annual';
					elseif ismember(A.time(i),mstr)
						tidx{i} = find(strcmp(mstr,A.time(i))==1);
						tstr{i} = mstr{i};
					elseif A.time(i) <= obj.region.nt
						tidx{i} = A.time(i);
						tstr{i} = ['Record ',num2str(A.time(i))]; 
					end
				end
			end

			% - get lonbounds/latbounds		
			lonbounds = [floor(obj.region.minlon_rho) ceil(obj.region.maxlon_rho)];
			latbounds = [floor(obj.region.minlat_rho) ceil(obj.region.maxlat_rho)];
			
			% - get levels or depths
			if isempty(lidx) & ~isempty(didx)
				zidx = didx; zstr = dstr;
			elseif isempty(didx) & ~isempty(lidx)
				zidx = lidx; zstr = lstr;
			else
				zidx = [0]; zstr = {'SFC'};
			end
			
			% - make plot for each variable, depth, and time combo
			pcnt = [0];
			for i = 1:length(A.vars)
				for j = 1:length(zidx)
					for k = 1:length(A.time);
						try % 3D
							data = squeeze(obj.romsData.(A.vars{i}).data(:,:,zidx(j),tidx{k}));
						catch % 2D
							data = squeeze(obj.romsData.(A.vars{i}).data(:,:,tidx{k}));
						end
						
						% - apply averaging
						data = nanmean(data,3);
	
						% - mask coast
						data = data .* obj.region.mask_rho;

						% - colorbar limits
						if ~isempty(A.caxis)
							datarng = [A.caxis(1) A.caxis(2)];
						else
							datarng = prctile(data(:),[1 99]);
						end
						[axisLims axisTicks]   = romsMaster.oom_levs(datarng);	
						lvls                   = [axisLims(1):(diff(axisLims)/255):axisLims(2)];	
						data(data>axisLims(2)) = axisLims(2);
						data(data<axisLims(1)) = axisLims(1);
						
						% - initiate figure
						fig(i,j,k)            = piofigs('mfig',1);
						set(0,'CurrentFigure',fig(i,j,k));
						[ax(i,j,k)] = map_plot(fig(i,j,k),obj.region.lon_rho,obj.region.lat_rho);
						m_contourf(obj.region.lon_rho,obj.region.lat_rho,data,lvls,'LineStyle','none');
						cb(i,j,k)   = colorbar;
						title([tstr{k},' ',obj.romsData.(A.vars{i}).name,' @ ',zstr{j}]);
						ylabel(cb(i,j,k),obj.romsData.(A.vars{i}).units);
	
						% - axis limits
						cb(i,j,k).XLim       	     = axisLims;
						cb(i,j,k).XTick              = axisTicks;
						cb(i,j,k).XTickLabel         = axisTicks;
						caxis([axisLims]);
					
						% - save
						if isempty(A.plttype)
							pcnt = pcnt+1;
							pltjpg(pcnt);
						elseif strcmp(A.plttype,'tmpfigs')
							if i == 1 & j == 1 & k == 1
								pltshow('mode',1);
								pltshow;
							else
								pltshow;
							end
						elseif zidx > 0
							fdir  = obj.paths.plots.roms.zsurfacefigs;
							fname = [tstr{k},'_',A.vars{i},'_',zstr{j}];
							disp([fdir,fname])
							export_fig([fdir,fname],['-',A.plttype],'-p0.05'); 
						else
							fdir = obj.paths.plots.roms.surfacefigs;
							fname = [tstr{k},'_',A.vars{i},'_',zstr{j}];
							disp([fdir,fname])
							export_fig([fdir,fname],['-',A.plttype],'-p0.05'); 
						end	
						
						% - rinse and repeat
					end
				end
			end
		end % end methods plotData
	
		%--------------------------------------------------------------------------------
      		function obj = sliceROMS(obj,varargin);
			% -------------------
			% Takes 2D depth slice of ROMS along a given latitude or longitude
			% Options are set by user
			% 
			% Usage:
			% - obj = sliceROMS(obj,varargin);
			% 
			% Inputs (varargin):
			% - (empty) = Choose inputs
			% - choice = 'lat' or 'lon' depending on slice
			% - deg    = latitude or longitude to slice (0 - 360 for lon)
			% - vars   = ROMS variable to slice, as a string or cell array
			% - vali   = array size of 'vars', with 1 (yes) or 0 (no)...
			%	     ...use to also slice validation data
			% 
			% Example
			% - obj = sliceROMS(obj,'choice','lon','deg',0,'vars',{'temp','salt'},'vali',[1 0]);
			% 
			% This will slice temp and salt data at 0 degrees longitude, and will also return
			% sliced temperature (but not salinity) data from the validation dataset
        		% -------------------
		
			% - Grab user inputs
			A.choice  = [];
			A.deg     = [];
			A.vars    = [];
			A.vali    = [];
			A         = parse_pv_pairs(A,varargin);

			% - check that object initDiag has run
			try; obj.paths.diag.woamasks;
			catch; obj = initDiag(obj); 
			end
	
			% - Check inputs
			if ~isempty(A.choice) % lon or lat
				if ~ismember(A.choice,{'lon','lat'})
					disp(' ');
					disp('Choice = lat or lon ONLY')
					disp(' ');
					return
				end
			end
			if ~isempty(A.deg) % degree
				for i = 1:length(A.deg)
					if A.deg(i) <-90 | A.deg(i) > 360
						disp(' ');
						disp('Lon = 0:360, Lat = -90:90, try again');
						disp(' ');
						return
					end
				end
			end
			if ~isempty(A.vars) % variables
				for i = 1:length(A.vars)
					if ~ismember(A.vars{i},obj.info.var3d)
						disp(' ');
						disp([A.vars{i},'is not 3D or is not recognized']);
						disp(' ');
						return	
					end
				end					
			end
			if ~isempty(A.vali) % validation switches
				for i = 1:length(A.vali)
					if ~ismember(A.vali(i),[0 1]);
						disp(' ');
						disp('vali must be 1 or 0');
						return
					end
				end
				if length(A.vali) ~= length(A.vars)
					disp(' ');
					disp('vali must be the same size as vars')
					return
				end
			end

			% Process choice input
			if isempty(A.choice)
                        	A.choice = input('---------------------------\nLon(1) or Lat(2) slice:\n---------------------------\n>> ');
				if ~ismember(A.choice,[1 2])
					disp(' ');
					disp('Bad choice')
					disp(' ');
					return
				elseif A.choice == 1
					A.choice = 'lon';
				elseif A.choice == 2
					A.choice = 'lat';
				end
			end

			% Process degree input
			if isempty(A.deg)
				A.deg = input('---------------------------\nDegree of lon/lat slice (0:360 for lon), multiple sections OK:\n---------------------------\n>> ');
				for i = 1:length(A.deg)
					if A.deg(i) < -90 | A.deg(i) > 360
						disp(' ');
						disp('Number must be between 0:360 for lon, -90:90 for lat');
						disp(' ');
						return
					end
				end
			end

			% Process variable input
			if isempty(A.vars)	
				disp('-------------------------');
				disp('3-D Variables:');
				disp('-------------------------');
				for i = 1:length(obj.info.var3d);
					disp([num2str(i),'==',obj.info.var3d{i}]);	
				end
                        	q0 = input('---------------------------\nChoose from above:\n---------------------------\n>> ');
				if ~ismember(q0,[1:length(obj.info.var3d)]);
					disp('Bad choice');
					return
				else
					A.vars = obj.info.var3d([q0]);
				end	
			end

			% Process validation switch input
			if isempty(A.vali)
				q0 = input('---------------------------\nSlice validation data? (1 == Yes)\n---------------------------\n>> '); 
				if q0 == 1
					disp('-------------------------');
					disp('1 == yes, 0 == no');
					disp('-------------------------');
					for i = 1:length(A.vars)
						q1 = input([A.vars{i},' >>']);
						if ~ismember(q1,[0 1]);
							disp('1 or 0 only!')
							return
						else
							A.vali(i) = q1;
						end
					end
				else
					A.vali = zeros(1,length(A.vars));
				end
			end
			
			% - Grab longitude/latitude data
			lon  = obj.region.lon_rho;
			lat  = obj.region.lat_rho;

			% - Get mask, convert to double
			mask    = obj.region.mask_rho;
			mask(isnan(mask)) = 0;
			newmask = nan(size(mask));
			for i = 1:size(mask,1);
				for j = 1:size(mask,2);
					if mask(i,j)
						newmask(i,j) = 1;
					else
						newmask(i,j) = 0;
					end
				end
			end
			mask = double(newmask);
	
			% - Grab depth from z_avg file
			depth = ncread(obj.paths.zavg,'depth');

			% - Choose latitude or longitude
			if strcmp(A.choice,'lat')
				lonlat = lat; dmsn = obj.region.nx; nz = length(obj.grid.z_avg_dep); nl = obj.region.nt; ns = length(A.deg); dstr = ['latitude'];
			elseif strcmp(A.choice,'lon')
				lonlat = lon; dmsn = obj.region.ny; nz = length(obj.grid.z_avg_dep); nl = obj.region.nt; ns = length(A.deg); dstr = ['longitude'];
			end

			% - Dimensions to fill
			fillmat = [dmsn nz ns];	

			% - Copy mask
			mask = repmat(mask,1,1,nz);
		
			% - Start parpool
			%delete(gcp('nocreate'));
			%parpool(8);

			% - Load 3D ROMS data and slice
			for ff = 1:length(A.vars)

				% - Load z_avg variable
				data  = ncread(obj.paths.zavg,A.vars{ff},[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
				
				% - Blank bad z_avg data
				data(abs(data)>1e4) = NaN;
				
				% - Take slice of data for each month
				tmpslice = NaN(dmsn,nz,nl,ns);
				tmpdepth = NaN(dmsn,nz,nl,ns);
				tmpmask  = NaN(dmsn,nz,nl,ns);
                                disp(' '); disp(['Slicing ROMS ',A.vars{ff},' data @ ',num2str(A.deg),'deg ',dstr,'...']);disp(' ');
				%parfor rcrd = 1:nl
				for rcrd = 1:nl
					% - Display progress
					fprintf([num2str(rcrd),'...']);
                                        % - Only get data for that month
                                        tmpdata = squeeze(data(:,:,:,rcrd));
					% - Interp to lon/lat line
					for i = 1:dmsn
						for j = 1:ns
							% - Fill tmpslice
							for z = 1:nz
								if dmsn == obj.region.nx
									tmpslice(i,z,rcrd,j) = interp1(squeeze(lonlat(i,:)),squeeze(tmpdata(i,:,z)),A.deg(j));
									tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(i,:)),squeeze(mask(i,:,z)),A.deg(j));
								elseif dmsn == obj.region.ny
									tmpslice(i,z,rcrd,j) = interp1(squeeze(lonlat(:,i)),squeeze(tmpdata(:,i,z)),A.deg(j));
									tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(:,i)),squeeze(mask(:,i,z)),A.deg(j));
								end
							end
						end
                                        end
                                        tmpdepth(:,:,rcrd,:) = depth'.*ones(fillmat);
				end
                              
				% - Get lon/lat data (outside parfor)
                                tmpdeg  = NaN(dmsn,ns);
                                for i = 1:dmsn
					for j = 1:ns	
						if dmsn == obj.region.nx;
							tmpdeg(i,j)  = interp1(squeeze(lat(i,:)),squeeze(lon(i,:)),A.deg(j));
						elseif dmsn == obj.region.ny;
							tmpdeg(i,j)  = interp1(squeeze(lon(:,i)),squeeze(lat(:,i)),A.deg(j));
						end
					end
                                end
                                tmpdeg = repmat(tmpdeg,1,1,nz); tmpdeg = permute(tmpdeg,[1 3 2]);
				
				% - Apply mask
				tmpdepth = squeeze(tmpdepth(:,:,1,1));
				for j = 1:ns
					masktmp         = squeeze(tmpmask(:,:,1,j));
					tmp             = tmpdeg(:,:,j);
					tmp(masktmp==0) = NaN;
					tmpdeg(:,:,j)   = tmp;
					for rcrd = 1:nl
						tmp                  = squeeze(tmpslice(:,:,rcrd,j));
						tmp(masktmp==0)      = NaN;
						tmpslice(:,:,rcrd,j) = tmp;
					end
				end

				% - Save results
				idx  = strmatch(A.vars{ff},{obj.info.Variables.Name},'exact');
				iidx = strmatch('units',{obj.info.Variables(idx).Attributes.Name}); 
                                obj.romsData.(A.vars{ff}).slicedata  = tmpslice;
                                obj.romsData.(A.vars{ff}).slicedepth = tmpdepth;
                                obj.romsData.(A.vars{ff}).slicedeg   = tmpdeg;
                                if dmsn == obj.region.nx;
                                        obj.romsData.(A.vars{ff}).slicecoord = 'latitude';
                                elseif dmsn == obj.region.ny;
                                        obj.romsData.(A.vars{ff}).slicecoord = 'longitude';
                                end
                                obj.romsData.(A.vars{ff}).slicesect = A.deg;
				obj.romsData.(A.vars{ff}).units     = obj.info.Variables(idx).Attributes(iidx).Value;
                                fprintf(['\n']);
			end
			delete(gcp('nocreate'));

			% Clear tmp
			tmp = [];

			% - Check if validation data needs to be sliced
			if sum(A.vali) > 0

				% - Perform slices of each variable
				for ff = 1:length(A.vars);

					% - Check if validation is requested
					if A.vali(ff) == 0
						continue
					end

					disp(' '); disp(['Slicing validation ',A.vars{ff},' data...']);disp(' ');
					% - get path and coords for current variable
					curVar  = obj.paths.diag.(A.vars{ff});
					
					% - if file is empty, fill with NaNs and move on
					if isempty(curVar.file)
						if strcmp(A.choice,'lat');
							obj.diagData.(A.vars{ff}).slicedata  = nan(obj.region.nx,length(obj.grid.z_avg_dep),obj.region.nt,length(A.deg));
							obj.diagData.(A.vars{ff}).slicedepth = tmpdepth;
							obj.diagData.(A.vars{ff}).slicedeg   = tmpdeg;
							obj.diagData.(A.vars{ff}).slicecoord = 'latitude';
							obj.diagData.(A.vars{ff}).slicesect  = A.deg;
						elseif strcmp(A.choice,'lon');				
							obj.diagData.(A.vars{ff}).slicedata  = nan(obj.region.ny,length(obj.grid.z_avg_dep),obj.region.nt,length(A.deg));
							obj.diagData.(A.vars{ff}).slicedepth = tmpdepth;
							obj.diagData.(A.vars{ff}).slicedeg   = tmpdeg;
							obj.diagData.(A.vars{ff}).slicecoord = 'longitude';
							obj.diagData.(A.vars{ff}).slicesect  = A.deg;	 
						end
						continue
					end
					
					% - if file exists, load data
					if strcmp(curVar.type,'nc');
						tmp.lon   = ncread(curVar.file,curVar.lon); tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
						tmp.lat   = ncread(curVar.file,curVar.lat);
						tmp.depth = ncread(curVar.file,curVar.zvar);
						if nanmean(tmp.depth(:)) > 0
							tmp.depth = -tmp.depth;
						end
						tmp.data  = ncread(curVar.file,curVar.var);
						% convert predicted N2O from nmol/l to mmol/m3
						if strcmp(A.vars{ff},'N2O') 
							tmp.data = tmp.data./1000; % nMol/l --> mMol/m3
						end		
					elseif strcmp(curVar.type,'mat');
						tmp.lon   	   = load(curVar.file,curVar.lon);
						tmp.lat   	   = load(curVar.file,curVar.lat);
						tmp.depth 	   = load(curVar.file,curVar.zvar);
						tmp.data           = load(curVar.file,curVar.file);
						tmp.lon    	   = tmp.lon.(curVar.lon);
						tmp.lat            = tmp.lat.(curVar.lat);
						tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
						tmp.depth          = tmp.depth.(curvar.zvar);
					end

					% - make lon/lat data gridded if it isnt setup that way
					if sum(size(tmp.lon(:,:))==1) > 0
						[tmp.latg,tmp.long,tmp.depthg] = meshgrid(tmp.lat,tmp.lon,tmp.depth);
					end
					
					% - go through each requested slice
					for j = 1:ns
						% - go through each time record
						for rcrd = 1:nl
							% - clear variables to fill
							tmp.latr   = []; 
							tmp.lonr   = [];
							tmp.depthr = [];
							tmp.datar  = [];
							% - record progress
							fprintf([num2str(rcrd),'...']);
							% - Get reduced matrix to feed into scatteredInterpolant
							if strcmp(A.choice,'lat');
								idx = find(abs([tmp.lat-A.deg(j)]) == min(abs([tmp.lat-A.deg(j)])));
								tmp.latr   = tmp.latg(:,idx,:);
								tmp.lonr   = tmp.long(:,idx,:);
								tmp.depthr = tmp.depthg(:,idx,:);
								tmp.datar  = tmp.data(:,idx,:,rcrd);
								if length(idx)==2
									tmp.latr   = squeeze(nanmean(tmp.latr,2));
									tmp.lonr   = squeeze(nanmean(tmp.lonr,2));
									tmp.depthr = squeeze(nanmean(tmp.depthr,2));
									tmp.datar  = squeeze(nanmean(tmp.datar,2));
								else
									tmp.latr   = squeeze(tmp.latr);
									tmp.lonr   = squeeze(tmp.lonr);
									tmp.depthr = squeeze(tmp.depthr);
									tmp.datar  = squeeze(tmp.datar);
								end
							elseif strcmp(A.choice,'lon');
							% - Get reduced matrix to feed into scatteredInterpolant
								idx = find(abs([tmp.lon-A.deg(j)]) == min(abs([tmp.lon-A.deg(j)])));
								tmp.latr   = tmp.latg(idx,:,:);
								tmp.lonr   = tmp.long(idx,:,:);
								tmp.depthr = tmp.depthg(idx,:,:);
								tmp.datar  = tmp.data(idx,:,:,rcrd);
								if length(idx)==2
									tmp.latr   = squeeze(nanmean(tmp.latr,1));
									tmp.lonr   = squeeze(nanmean(tmp.lonr,1));
									tmp.depthr = squeeze(nanmean(tmp.depthr,1));
									tmp.datar  = squeeze(nanmean(tmp.datar,1));
								else
									tmp.latr   = squeeze(tmp.latr);
									tmp.lonr   = squeeze(tmp.lonr);
									tmp.depthr = squeeze(tmp.depthr);
									tmp.datar  = squeeze(tmp.datar);
								end
							end
							% - build interpolant and interpolate
							if strcmp(A.choice,'lat');
								tmp.lonr      = tmp.lonr(:);
								tmp.depthr    = tmp.depthr(:);
								tmp.datar     = tmp.datar(:); 
								F             = scatteredInterpolant(tmp.lonr(~isnan(tmp.datar)),...
										tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
								tmp.out{rcrd} = F(tmpdeg(:,:,j),tmpdepth);
							elseif strcmp(A.choice,'lon');
								tmp.latr      = tmp.latr(:);
								tmp.depthr    = tmp.depthr(:);
								tmp.datar     = tmp.datar(:);
								F             = scatteredInterpolant(tmp.latr(~isnan(tmp.datar)),...
										tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
								tmp.out{rcrd} = F(tmpdeg(:,:,j),tmpdepth);
							end
						end

						% - Save results
						obj.diagData.(A.vars{ff}).slicedata(:,:,:,j)  = cat(3,tmp.out{:});
						obj.diagData.(A.vars{ff}).slicedeg(:,:,j)     = tmpdeg(:,:,j);
					end
					if dmsn == obj.region.nx;
						obj.diagData.(A.vars{ff}).slicecoord = 'latitude';
					elseif dmsn == obj.region.ny;
						obj.diagData.(A.vars{ff}).slicecoord = 'longitude';
					end
					obj.diagData.(A.vars{ff}).slicedepth = tmpdepth;
					obj.diagData.(A.vars{ff}).slicesect  = A.deg;
					fprintf(['\n']);
				end
			end

	        	% - calculate additional diags (e.g. N*)
			if sum(ismember(A.vars,{'NO3','PO4'}))>=2;
            			disp('Calculating N*')
				obj.romsData.nstar           = obj.romsData.NO3;
				obj.romsData.nstar.slicedata = obj.romsData.NO3.slicedata - 16.*obj.romsData.PO4.slicedata + 2.9;
				if sum(A.vali)>0
					obj.diagData.nstar           = obj.diagData.NO3;
					obj.diagData.nstar.slicedata = obj.diagData.NO3.slicedata - 16.*obj.diagData.PO4.slicedata + 2.9;
				end
			end
		end % end methods sliceROMS

		%--------------------------------------------------------------------------------
		function obj = getMLD(obj,varargin)
        		% ------------------
			% Computes mixed layer depth from ROMS data for each month
			% Compares against an ARGO only or combo (WOCE/NODC/ARGO) dataset
                        % 
                        % Usage:
                        % - obj = getMLD(obj,varargin);
                        % 
                        % Inputs (varargin):
                        % - outer = degrees of lon/lat to pad interpolant, default = 3
			%
			% Example
			% - obj = getMLD(obj);
			% - obj = getMLD(obj,'outer',5);
			% -------------------
			
			% - process inputs
			A.outer   = 3;
         		A         = parse_pv_pairs(A,varargin);
		
			% - check that object is initialized(init) and initDiag has run
			try; x = obj.info.var3d; obj.paths.diag.woamasks; catch; obj = init(obj); obj = initDiag(obj); end

			% - load depth for z_avg data
			depth = -(ncread(obj.paths.zavg,'depth'));
			
			% - load temp/salt data then calculate density
			temp = ncread(obj.paths.zavg,'temp',[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
			salt = ncread(obj.paths.zavg,'salt',[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
			% - blank bad data
			temp(abs(temp)>1e4) = NaN;
			salt(abs(salt)>1e4) = NaN; 
			dens = real(sw_dens0(salt,temp));

			% - initiate matrices and dimensions
                        nx      = obj.region.nx;
                        ny      = obj.region.ny;
			nz      = length(depth);
                        nl      = obj.region.nt;
			fillmat = [nx ny nz];
 
			% - Compute ROMS MLD for each month
                        disp(' ');
                        disp('Calculating ROMS MLD for time record');
                        fprintf(['\n record:']);
                        MLD     = NaN(nx,ny,nl);
			for rcrd = 1:nl
                                % - Announce progress
				fprintf([num2str(rcrd),'...']);
                                % - Squeeze dens for that time record
				tmpdens = squeeze(dens(:,:,:,rcrd));
				% - Get index of 10m
				zind = find(depth == 10);
				% - Find density at 10m
				dens10 = squeeze(tmpdens(:,:,zind));
				% - Go through each lat/lon and find where dens > (dens10 + 0.03)
				for i = 1:nx
					for j = 1:ny
						if ~isnan(tmpdens(i,j,zind))
							dind = find(tmpdens(i,j,:) >= (dens10(i,j)+0.03));
							if ~isempty(dind);
								MLD(i,j,rcrd) = depth(dind(1));
							end
						end
					end
				end
			end
			delete(gcp);
			fprintf(['\n']);

			% - Save ROMS data (make positive)
			obj.romsData.MLD.surfdata = MLD;

			% - Grab both MLD validation data sets
			for dat = 1:2
				% - Calculate MLD from WOCE, NODC, and ARGO data from 1961 - 2008
				if dat == 1
					% - Get netcdf file
					fname  = ['/data/project1/demccoy/ROMS/validation/MLD/mld_DR003_c1m_reg2.0.nc'];
					var    = ['mld'];
					% - Get coordinates
					tmp.lon = ncread(fname,'lon');
					tmp.lat = ncread(fname,'lat');
				elseif dat == 2
					% - Get netcdf file
					fname = ['/data/project1/demccoy/ROMS/validation/MLD/Argo_mixedlayers_monthlyclim_12112019.nc'];
					var   = ['mld_da_mean'];
					% - Get coordinates
					tmp.lon = ncread(fname,'lon');
					tmp.lat = ncread(fname,'lat');
					% - Fix longitude (-180:180 --> 0:360)
					indn = find(tmp.lon<0);
					indp = find(tmp.lon>0);
					tmp.lon = [tmp.lon(indp);tmp.lon(indn)];
					tmp.lon = obj.lon360(tmp.lon);
				end

				% - Get in grid form
				if size(tmp.lon,2)==1
					[tmp.lat,tmp.lon] = meshgrid(tmp.lat,tmp.lon);
				end

				% - get indeces for reduced domain interpolation
            			% - note: this breaks if ROMS boundary longitude is close to 0 or 360
  				idx = find(tmp.lon(:) > obj.region.minlon_rho-A.outer ...
            	 			 & tmp.lon(:) < obj.region.maxlon_rho+A.outer ...
              				 & tmp.lat(:) > obj.region.minlat_rho-A.outer ...
              				 & tmp.lat(:) < obj.region.maxlat_rho+A.outer);
					
				% - get climatology for each month
				disp(' ');
				if dat == 1
					disp(['Calculating MLD from ',obj.plots.diag.mld_ifremer.data.varsource]);
				elseif dat == 2
					disp(['Calculating MLD from ',obj.plots.diag.mld_argo.data.varsource]);
				end
				fprintf(['\n month:'])
				for mnth = 1:12
					fprintf([num2str(mnth),'...']);
					if dat == 1
						tmp.mld = squeeze(ncread(fname,var,[1,1,mnth],[inf,inf,1]));
						% - ignore mask and missing values
						tmp.mld(tmp.mld==10^9)  = NaN;
						tmp.mld(tmp.mld==-9999) = NaN;
					elseif dat == 2
						tmp.mld = squeeze(ncread(fname,var,[mnth,1,1],[1,inf,inf]));
						% - reorganize data due to longitude issues
						tmp.mld = [tmp.mld(indp,:);tmp.mld(indn,:)];
					end
			      	        % - build interpolant and interpolate
             				F 			       = scatteredInterpolant(double(tmp.lon(idx)),double(tmp.lat(idx)),...
									 double(tmp.mld(idx)),'linear','nearest');
             				tmpout{mnth} 		       = F(double(obj.region.lon_rho),double(obj.region.lat_rho));
				end
				fprintf(['\n']);

				if dat == 1
					% - Save validation MLD (make all positive)
                			obj.diagData.MLD_ifremer.surfdata = cat(3, tmpout{:});
         				obj.diagData.MLD_ifremer.surfdata = single(obj.diagData.MLD_ifremer.surfdata);
				elseif dat == 2
					% - Save validation MLD (make all positive)
                			obj.diagData.MLD_argo.surfdata = cat(3, tmpout{:});
         				obj.diagData.MLD_argo.surfdata = single(obj.diagData.MLD_argo.surfdata);
				end
			end
		end % end method getMLD

		%--------------------------------------------------------------------------------
		function obj = cravatteCMPR(obj)
			% -----------------------
			% Compares equatorial current structure with Cravette2017 results
			% 
			% Returns depth slices of u-velocity along 195 longitude
			% -----------------------

			% - Load Cravette(2017) mean zonal currents from 179E - 160W
			fname = '/data/project1/data/Tropical_Currents_Cravatte_2017/Mean_zonal_currents_179E-160W.cdf';
			datu  = ncread(fname,'U_SADCPD');
			depu  = ncread(fname,'DEPTH');
			latu  = ncread(fname,'LATI');

			% - Grid that shit
			[depu,latu] = meshgrid(-depu,latu);

			% - Get 3D u velocity data
			romsu  = ncread(obj.paths.zavg,'u',[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
			romsln = obj.region.lon_rho;
			romslt = obj.region.lat_rho;
			romsdp = obj.grid.woa1p0.depth;

			% - get reduce longitude matrix (u points)
			[a,b] = size(romsln);
			TMPLN = NaN(a-1,b);
			TMPLT = NaN(a-1,b);
			for i = 1:b
				tmpln(:,i) = romsln(:,i);
				tmplt(:,i) = romslt(:,i);
				TMPLN(:,i) = [tmpln(1:end-1,i) + tmpln(2:end,i)]/2;
				TMPLT(:,i) = [tmplt(1:end-1,i) + tmplt(2:end,i)]/2;
			end
			romsln = TMPLN; clear TMPLN tmpln;
			romslt = TMPLT; clear TMPLT tmplt;

			% - Take slice of data for each month
			dmsn = obj.region.ny;
			nl   = obj.region.nt;
			nz   = obj.grid.woa0p25.depth; nzz = length(nz);
			disp(' '); disp(['Slicing ROMS u data']);disp(' ');
			tmpslice = NaN(dmsn,nzz,nl);
			%parfor mnth = 1:nl
			for mnth = 1:nl
				fprintf([num2str(mnth),'...']);
				% - Only get data for that month
				tmpdata = squeeze(romsu(:,:,:,mnth));
				% - Interp to lon/lat line
				for i = 1:dmsn
					% - Fill tmpslice			
					for z = 1:nzz;
						tmpslice(i,z,mnth) = interp1(squeeze(romsln(:,i)),squeeze(tmpdata(:,i,z)),195);
					end
				end
			end

			% - Set up grid
			for i = 1:dmsn
				tmpdeg(i) = interp1(squeeze(romsln(:,i)),squeeze(romslt(:,i)),195);
			end
		 	tmpdepth = -nz;
			[tmpdepth,tmpdeg] = meshgrid(tmpdepth,tmpdeg);
			
			% Interpolate results to validation grid	
			tmpdepth = tmpdepth(:);
			tmpdeg   = tmpdeg(:);
			for mnth = 1:12
				tmpu        = tmpslice(:,:,mnth);
				idx         = find(isnan(tmpu)==0);
				F           = scatteredInterpolant(tmpdeg(idx),tmpdepth(idx),tmpu(idx));
				gridu{mnth} = F(latu,depu);
			end

			% Save results in format for plotZDiags
			obj.romsData.u.slicedata  = cat(3,gridu{:});
			obj.romsData.u.slicedeg   = latu;
			obj.romsData.u.slicedepth = depu;
			obj.romsData.u.slicecoord = 'longitude';
			obj.romsData.u.slicesect  = 195;
			obj.diagData.u.slicedata  = datu/100;
			obj.diagData.u.slicedeg   = latu;
			obj.diagData.u.slicedepth = depu;
			obj.diagData.u.slicecoord = 'longitude';
			obj.diagData.u.slicesect  = 195;
		end % end method cravatteCMPR

		%--------------------------------------------------------------------------------
      		function obj = interp2DData(obj,varargin)
        		% ------------------
			% Loads and interpolates monthly 2D validation data to the ROMS grid
			% See initDiags for available fields
			%
			% Usage:
			% - obj = interp2DData(obj,varargin);
			%
			% Inputs:
			% - vars  = validation variable to load, paths set in initDiag
			% - outer = degrees of lon/lat to pad interpolant, default = 3
			% - depth = depths to interpolate data, will round to nearest z_avg_dep
			%
			% Example:
			% - obj = interp2DData(obj,'vars',{'temp','salt','O2','NO3'},'depth',0);
        		% -------------------

			% - define variables to interpolate
	     		A.vars  = [];  % variables
         		A.outer = [3]; % interpolant padding (degrees)
			A.depth = [];  % depth surface (m), can use 0 for surface
         		A=parse_pv_pairs(A,varargin);

			% - check that object is initialized(init) and initDiag has run
			try; x = obj.info.var3d; obj.paths.diag.woamasks; catch; obj = init(obj); obj = initDiag(obj); end
			
			% - check inputs
			diagfields = fieldnames(obj.paths.diag); ind = find(strcmp(diagfields,'woamasks'))==1; diagfields(ind) = []; clear ind
			if ~isempty(A.vars)
				for i = 1:length(A.vars)
					if ~strcmp(A.vars{i},diagfields);
						disp(' ');
						disp([A.vars{i},' is not a diagnostic variable']);
						disp(' '); return
					end
				end
			end
			if ~isempty(A.depth)
				if A.depth < 0 | A.depth > max(obj.grid.z_avg_dep)
					disp(' ');
					disp('Check depth input');
					disp(' '); return
				end
				for i = 1:length(A.depth)
					diffd = abs(A.depth(i) - obj.grid.z_avg_dep);
					ind   = find(diffd == min(diffd));
					A.depth(i) = obj.grid.z_avg_dep(ind);
				end
			end

			% - process variable input
			if isempty(A.vars) % validation variables
				disp('-------------');
				for i = 1:length(diagfields)
					if i < 10
						disp([num2str(i),'  == ',diagfields{i}]);
					else
						disp([num2str(i),' == ',diagfields{i}]);
					end
				end
                        	q0 = input('---------------------------\nChoose from above:\n---------------------------\n>> ');
				if max(q0) <= length(diagfields) & min(q0) > 0
				else
					disp(' ');
					disp('Bad choice');
					disp(' '); return
				end
				A.vars = diagfields(q0);
			end	

			% - process depth input
			if isempty(A.depth) % depth 
				disp('-------------------------');
				disp('z_avg Depths:');
				for i = 1:length(obj.grid.z_avg_dep)
					if i < 10
						disp([num2str(i),'   == ',num2str(obj.grid.z_avg_dep(i)),'m']);
					elseif i < 100
						disp([num2str(i),'  == ', num2str(obj.grid.z_avg_dep(i)),'m']);
					else
						disp([num2str(i),' == ',  num2str(obj.grid.z_avg_dep(i)),'m']);
					end
				end
                        	q0 = input('---------------------------\nChoose from above:\n---------------------------\n>> ');
				if ~ismember(q0,[1:length(obj.grid.z_avg_dep)])
					disp(' ');
					disp('Bad choice');
					disp(' '); return
				end
				A.depth = obj.grid.z_avg_dep(q0);
			end

			% - get temporary output longitude and latitude
         		tmp.outlon = obj.region.lon_rho;
         		tmp.outlat = obj.region.lat_rho;

			% - process each variable
         		for ff = 1 : length(A.vars)		
       		 		fprintf(['\n Processing ', A.vars{ff}]);

				% - ignore empty fields for data like NH4
				if isempty(obj.paths.diag.(A.vars{ff}).file)
					disp(['No diagnostic data for ',A.vars{ff}]);
					disp(['Filling with NaNs']);
					if ismember(0,A.depth)
						obj.diagData.(A.vars{ff}).surfdata = nan(obj.region.ndim_xyt);
					end
					if max(A.depth)>0
						ll  = length(A.depth(A.depth>0));
						obj.diagData.(A.vars{ff}).zsurfdata  = nan(obj.region.nx,obj.region.ny,ll,obj.region.nt);
						obj.diagData.(A.vars{ff}).zsurfdepth = A.depth(A.depth>0);
					end
					continue
				end

			      	% - get path and coords for current variable
      	  			curVar  = obj.paths.diag.(A.vars{ff});
				if strcmp(curVar.type,'nc');
  		     	        	tmp.lon = romsMaster.lon360(ncread(curVar.file,curVar.lon));
      	 	     			tmp.lat = ncread(curVar.file,curVar.lat);
				elseif strcmp(curVar.type,'mat');
					tmp.lon = load(curVar.file,curVar.lon);
					tmp.lat = load(curVar.file,curVar.lat);
					tmp.lon = tmp.lon.(curVar.lon);
					tmp.lat = tmp.lat.(curVar.lat);
					tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
				end

				% - make lon/lat data gridded if it isnt setup that way
       		     		if sum(size(tmp.lon(:,:))==1) > 0
       		        		[tmp.lat,tmp.lon] = meshgrid(tmp.lat, tmp.lon);
       		     		end;
            			tmp = romsMaster.struct2double(tmp);
            
            			% - get indeces for reduced domain interpolation
            			% - note: this breaks if ROMS boundary longitude is close to 0 or 360
  			        idx = find(tmp.lon(:) > obj.region.minlon_rho-A.outer ...
            	 			 & tmp.lon(:) < obj.region.maxlon_rho+A.outer ...
              				 & tmp.lat(:) > obj.region.minlat_rho-A.outer ...
              				 & tmp.lat(:) < obj.region.maxlat_rho+A.outer);

				% - initialize interpolated output
				if ismember(0,A.depth)
            				obj.diagData.(A.vars{ff}).surfdata = nan(obj.region.ndim_xyt);
				end
				if max(A.depth)>0
					obj.diagData.(A.vars{ff}).zsurfdata = nan(obj.region.nx,obj.region.ny,length(A.depth(A.depth ~= 0)),obj.region.nt);
				end

    			        % - interpolate data on ROMS coords for each month
				lvls = length(A.depth);
        		  	fprintf('\n month:');
				%delete(gcp('nocreate'));
				%parpool(8);
        	    		%parfor l = 1 : 12
				for l = 1:12
        	       			tmpdata = []; tmpdepth = [];
           	   			fprintf([num2str(l),'...']);

		              		% - get data
					if strcmp(curVar.type,'nc')
             					if strcmp(curVar.dim,'xyt')
               	   					tmpdata  = squeeze(ncread(curVar.file,curVar.var,...
               	      					 	  [1,1,l],[inf,inf,1]));
							tmpdepth = 0;
               					elseif strcmp(curVar.dim,'xyzt')
							tmpdepth = ncread(curVar.file,curVar.zvar);
							tmpdata  = squeeze(ncread(curVar.file,curVar.var,...
								   [1,1,1,l],[inf,inf,inf,1]));
						end
					elseif strcmp(curVar.type,'mat')
						if strcmp(curVar.dim,'xyt');
							tmpdata  = load(curVar.file,curVar.var);
							tmpdata  = tmpdata.(curVar.var)(:,:,l);
							tmpdepth = 0;
						elseif strcmp(curVar.dim,'xyzt')
							tmpdepth = load(curVar.file,curVar.zvar);
							tmpdepth = tmpdepth.(curVar.zvar);
							tmpdata  = load(curVar.file,curVar.var);
							tmpdata  = tmpdata.(curVar.var)(:,:,:,l);
						end
					end
		
					% Convert nMol/l --> mMol/m3 for N2O
					if strcmp(A.vars{ff},'N2O') 
						tmpdata = tmpdata./1000; % nMol/l --> mMol/m3
					end		

					% - build interpolant and interpolate over all depths
					for z = 1:lvls
						data  = tmpdata;
						depth = tmpdepth; 
						if nanmean(depth) < 0
							depth = -depth;
						end
						zind  = find(abs([depth-A.depth(z)]) == min(abs([depth-A.depth(z)])));
						if ~strcmp(curVar.dim,'xyt')
							data  = data(:,:,zind);
						end
						if strcmp(curVar.dim,'xyt') & z > 1
							disp(['No z-data for ',A.vars{ff},', skipping']);
							continue
						end
	 					% - interpolate
						F         = scatteredInterpolant(double(tmp.lon(idx)),double(tmp.lat(idx)),...
                                                                                 double(data(idx)),'linear','nearest');
						tmpout{z,l} = F(double(tmp.outlon),double(tmp.outlat));
                                        	tmpout{z,l}(isnan(obj.region.mask_rho)) = nan;
					end
     				end
			
				% - save interpolated data
				surfind  = find(A.depth == 0);
				zsurfind = find(A.depth > 0);
				if ~isempty(surfind)
	                		obj.diagData.(A.vars{ff}).surfdata = cat(3, tmpout{surfind,:});
           				obj.diagData.(A.vars{ff}).surfdata = single(obj.diagData.(A.vars{ff}).surfdata);
				end
				if ~isempty(zsurfind)
					for ll = 1:length(zsurfind);
						obj.diagData.(A.vars{ff}).zsurfdata(:,:,ll,:)  = cat(3, tmpout{zsurfind(ll),:});
						obj.diagData.(A.vars{ff}).zsurfdata(:,:,ll,:)  = single(obj.diagData.(A.vars{ff}).zsurfdata(:,:,ll,:));
					end
					obj.diagData.(A.vars{ff}).zsurfdepth = A.depth(zsurfind);
				end
			end
	        	fprintf('\n');

	        	% - calculate additional diags (e.g. N*)
			if sum(ismember(A.vars,{'NO3','PO4'}))>=2;
            			disp('Calculating N*')
				if ~isempty(surfind)
					obj.diagData.nstar.surfdata = obj.diagData.NO3.surfdata - 16.*obj.diagData.PO4.surfdata + 2.9;
				end
         			if ~isempty(zsurfind)
					obj.diagData.nstar.zsurfdata = obj.diagData.NO3.zsurfdata - 16.*obj.diagData.PO4.zsurfdata + 2.9;
				end
			end
		end % end method interp2DData

		%--------------------------------------------------------------------------------
 	        function obj = get2DDiags(obj,var) 
        		% ------------------
			% Prepare surface ROMS variables for plotting against diagnostic data
			% For data with z-levels, averages data over surfNlevs, defined in:
			% obj.plots.diag.surfNlevs (called in initDiag).
			% If input is ws (wind stress) or wsc (wind stress curl), calculate it!
			% Called in plot2DDiags.m
			% 
			% Usage:
			% - obj = get2DDiags(obj,var)
			%
			% Inputs:
			% - var = ROMS variable to load and/or average
			% 
			% Other Inputs allowed:
			% - 'nstar', 'npp', 'ssh', 'zws', 'mws','wsc','ws'
			%
			% Example
			% - obj = get2DDiags(obj,'temp');
			% ------------------

			% - squeeze variables	
			if strcmp(var,'nstar')
 				obj.romsData.(var).surfdata = squeeze(nanmean(ncread(obj.paths.zavg,obj.plots.diag.NO3.roms.var,...
               			    			      [obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,...
							      obj.plots.diag.surfNlevs,inf])- 16.* ...
               			    			      ncread(obj.paths.zavg,obj.plots.diag.PO4.roms.var,...
               			    			      [obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,...
          						      obj.plots.diag.surfNlevs,inf]),3)) + 2.9;
			elseif strcmp(var,'npp')
            			tmp                         = ncread(obj.paths.avg,'TOT_PROD',[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],...
                                                              [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
            			tmpHz                       = obj.region.Hz; 
            			obj.romsData.(var).surfdata = squeeze(nansum(tmp.*tmpHz,3)).*3600*24*12;
			elseif ismember(var,{'ssh','zws','mws'})
				obj.romsData.(var).surfdata = ncread(obj.paths.avg,obj.plots.diag.(var).roms.var,[obj.region.lon_lim(1),obj.region.lat_lim(1),1],...
    							      [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf]);
				if strcmp(var,'ssh');
					slacorr = nanmedian(obj.diagData.ssh.surfdata(:)) - nanmedian(obj.romsData.ssh.surfdata(:));
					disp(' '); disp(['Adding correction of ',num2str(slacorr),'m to ROMS SSH']);
					obj.romsData.ssh.surfdata = obj.romsData.ssh.surfdata + slacorr;
				end
			elseif ismember(var,{'ws','wsc'})
				tmpu   = ncread(obj.paths.avg,obj.plots.diag.zws.roms.var,[obj.region.lon_lim(1),obj.region.lat_lim(1),1],...
                                         [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf]);
				tmpv   = ncread(obj.paths.avg,obj.plots.diag.mws.roms.var,[obj.region.lon_lim(1),obj.region.lat_lim(1),1],...
                                         [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf]);
				tmpang = obj.region.angle; 
				if strcmp(var,'ws')
					[obj.romsData.ws.surfdata,~]  = obj.WindStress(tmpu,tmpv,obj.region.lon_rho,obj.region.lat_rho,tmpang);
				elseif strcmp(var,'wsc')
					[~,obj.romsData.wsc.surfdata] = obj.WindStress(tmpu,tmpv,obj.region.lon_rho,obj.region.lat_rho,tmpang);
				end
			else
 		        	obj.romsData.(var).surfdata = squeeze(nanmean(ncread(obj.paths.zavg,obj.plots.diag.(var).roms.var,...
 	                	    			      [obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],[diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,...
 							      obj.plots.diag.surfNlevs,inf]),3));
			end	
      		end % end method get2DDiags

		%--------------------------------------------------------------------------------
		function obj = plotZDiags(obj,varargin);
			% -----------------------
			% - Plot depth sections and comparisons extracted from sliceROMS
			% ------------------------
			
			A.vars  = {'temp','salt','O2','NO3','PO4','nstar'};
			A.xlims = [];
			A.zlims = [];
			A.plt   = 1; % final pdf
			A       = parse_pv_pairs(A,varargin);

			% - Gather savepath
			savepath = obj.paths.plots.diag.zsecfigs;

			% - Add nstar to plots if NO3 and PO4 are selected
			if ismember('NO3',A.vars) & ismember('PO4',A.vars) & ~ismember('nstar',A.vars)
				A.vars = [A.vars,{'nstar'}];
			end
			
			% - Loop through variables and plot!
			for ff = 1 : length(A.vars)
				disp(['Plotting ',A.vars{ff}]);
				% - Loop through slices
				for d = 1:length(obj.romsData.(A.vars{ff}).slicesect);
					disp(['slice #',num2str(d)]);
					% - Grab data, metadata, and options
					tmp.roms = obj.romsData.(A.vars{ff}).slicedata(:,:,:,d);
					tmp.data = obj.diagData.(A.vars{ff}).slicedata(:,:,:,d);
					gridopt  = obj.romsData.(A.vars{ff});
					gridopt  = rmfield(gridopt,'slicedata');
					gridopt.slicedepth = gridopt.slicedepth;
					gridopt.slicedeg   = gridopt.slicedeg(:,:,d);
					gridopt.slicesect  = gridopt.slicesect(d);
					pltopt   = obj.plots.diag.(A.vars{ff});
					% annual
					tag='annual';
					tmproms = nanmean(tmp.roms,3);
					tmpdata = nanmean(tmp.data,3);
					romsMaster.std_Zplots(tmproms,tmpdata,gridopt,pltopt,savepath,tag,A.xlims,A.zlims,A.plt);
					if strcmp(A.vars{ff},'u');
						continue
					end
					% winter
					tag='DJF';
					tmproms = nanmean(tmp.roms(:,:,[12,1,2]),3);
					tmpdata = nanmean(tmp.data(:,:,[12,1,2]),3);
					romsMaster.std_Zplots(tmproms,tmpdata,gridopt,pltopt,savepath,tag,A.xlims,A.zlims,A.plt);
					% spring
					tag='MAM';
					tmproms = nanmean(tmp.roms(:,:,[3,4,5]),3);
					tmpdata = nanmean(tmp.data(:,:,[3,4,5]),3);
					romsMaster.std_Zplots(tmproms,tmpdata,gridopt,pltopt,savepath,tag,A.xlims,A.zlims,A.plt);
					% summer
					tag='JJA';
					tmproms = nanmean(tmp.roms(:,:,[6,7,8]),3);
					tmpdata = nanmean(tmp.data(:,:,[6,7,8]),3);
					romsMaster.std_Zplots(tmproms,tmpdata,gridopt,pltopt,savepath,tag,A.xlims,A.zlims,A.plt);
					% fall
					tag='SON';
					tmproms = nanmean(tmp.roms(:,:,[9,10,11]),3);
					tmpdata = nanmean(tmp.data(:,:,[9,10,11]),3);
					romsMaster.std_Zplots(tmproms,tmpdata,gridopt,pltopt,savepath,tag,A.xlims,A.zlims,A.plt);
				end					
			end
		end % end method plotZDiags

		%--------------------------------------------------------------------------------
      		function obj = plot2DDiags(obj,varargin) 
        		% ------------------
			% - Plots annual and seasonal maps of observed, simulated NO3, PO4, CHL, NPP
			% - as well as their difference. The scale for CHL is logarithmic 
			% ------------------

			A.vars   = {'NO3','nstar','PO4','chl','nppcbpm','nppvgpm',... % BGC VARS
				    'temp','salt','ssh','zws','mws','ws','wsc'};      % PHY VARS
         		A.depth  = 0; % surface
			A.plt    = 1; % final pdf
			A        = parse_pv_pairs(A,varargin);
         	
			% - Get index of 'fill' data
			mask = obj.region.mask_rho;

			% - Add nstar to plots if NO3 and PO4 are selected
			if ismember('NO3',A.vars) & ismember('PO4',A.vars) & ~ismember('nstar',A.vars)
				A.vars = [A.vars,{'nstar'}];
			end
	
			% - Loop through variables and plot
			for ff = 1 : length(A.vars)
            			fprintf(['\n Plotting ', A.vars{ff}]);
            			if sum(strcmp(A.vars{ff},'nppcbpm')) > 0 | sum(strcmp(A.vars{ff},'nppvgpm')) > 0 | sum(strcmp(A.vars{ff},'nppcafe')) > 0
               				tmpvar = 'npp';
            			else
               				tmpvar = A.vars{ff};
            			end
				% - Loop through depths
				ddcnt = 0;
				for d = 1:length(A.depth)
					if A.depth(d) == 0
						if ~isfield(obj.romsData,tmpvar)
							try
								obj = get2DDiags(obj,tmpvar);
							catch
								disp(' '); disp([A.vars{ff},' missing from ROMS file'])
								continue
							end
						end	
						tmp.roms = obj.romsData.(tmpvar).surfdata;
						tmp.data = obj.diagData.(A.vars{ff}).surfdata;
						savepath = obj.paths.plots.diag.surfacefigs;
					elseif ~isnan(A.depth(d))
						if strcmp(tmpvar,'nstar')==1
							tmpdepth 		     = ncread(obj.paths.diag.NO3.file,obj.paths.diag.NO3.zvar);
							zind     		     = find(tmpdepth == A.depth(d));
							tmpno3   		     = squeeze(ncread(obj.paths.zavg,'NO3',[obj.region.lon_lim(1) obj.region.lat_lim(1) zind 1],...
										       [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,zind,1]));
							tmppo4   		     = squeeze(ncread(obj.paths.zavg,'PO4',[obj.region.lon_lim(1) obj.region.lat_lim(1) zind 1],...
										       [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,zind,1]));
							obj.romsData.nstar.zsurfdata = tmpno3 - 16.*tmppo4 + 2.9;
						else
							tmpdepth                        = obj.grid.z_avg_dep;
							zind                            = find(tmpdepth == A.depth(d));
							obj.romsData.(tmpvar).zsurfdata = squeeze(ncread(obj.paths.zavg,obj.plots.diag.(tmpvar).roms.var,...
											  [obj.region.lon_lim(1) obj.region.lat_lim(1) zind 1],...
											  [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 1 inf]));
						end
						ddcnt    = ddcnt + 1;
						tmp.roms = obj.romsData.(tmpvar).zsurfdata;
						tmp.data = obj.diagData.(A.vars{ff}).zsurfdata;
						savepath = obj.paths.plots.diag.zsurfacefigs;
					elseif isnan(A.depth(d)) & strcmp(A.vars{ff},'mld_ifremer')
						tmp.roms = obj.romsData.MLD.surfdata;
						tmp.data = obj.diagData.MLD_ifremer.surfdata;
						savepath = obj.paths.plots.diag.mldfigs;
					elseif isnan(A.depth(d)) & strcmp(A.vars{ff},'mld_argo');
						tmp.roms = obj.romsData.MLD.surfdata;
						tmp.data = obj.diagData.MLD_argo.surfdata;
						savepath = obj.paths.plots.diag.mldfigs;
					end
					options      = obj.plots.diag.(A.vars{ff});
					options.mask = mask;
					
					% - remove mask data (due to z_avg file)
					tr = []; td = [];
					for mnth = 1:12
						tr{mnth}          = tmp.roms(:,:,mnth);
						tr{mnth}(mask==0) = NaN;
						if size(tmp.data,4) == 12
							td{mnth}          = tmp.data(:,:,ddcnt,mnth);
							td{mnth}(mask==0) = NaN;
						else
							td{mnth}          = tmp.data(:,:,mnth);
							td{mnth}(mask==0) = NaN;
						end
					end
					tmp.roms = cat(3,tr{:});
					tmp.data = cat(3,td{:});
					% annual
					tag='annual';
					tmproms = nanmean(tmp.roms,3);
					tmpdata = nanmean(tmp.data,3);
					romsMaster.std_2dplots(tmproms,tmpdata,obj.region.lon_rho,obj.region.lat_rho,options,tag,savepath,A.depth(d),A.plt);
					% winter
					tag='DJF';
					tmproms = nanmean(tmp.roms(:,:,[12,1,2]),3);
					tmpdata = nanmean(tmp.data(:,:,[12,1,2]),3);
					romsMaster.std_2dplots(tmproms,tmpdata,obj.region.lon_rho,obj.region.lat_rho,options,tag,savepath,A.depth(d),A.plt);
					% spring
					tag='MAM';
					tmproms = nanmean(tmp.roms(:,:,[3,4,5]),3);
					tmpdata = nanmean(tmp.data(:,:,[3,4,5]),3);
					romsMaster.std_2dplots(tmproms,tmpdata,obj.region.lon_rho,obj.region.lat_rho,options,tag,savepath,A.depth(d),A.plt);
					% summer
					tag='JJA';
					tmproms = nanmean(tmp.roms(:,:,[6,7,8]),3);
					tmpdata = nanmean(tmp.data(:,:,[6,7,8]),3);
					romsMaster.std_2dplots(tmproms,tmpdata,obj.region.lon_rho,obj.region.lat_rho,options,tag,savepath,A.depth(d),A.plt);
					% fall
					tag='SON';
					tmproms = nanmean(tmp.roms(:,:,[9,10,11]),3);
					tmpdata = nanmean(tmp.data(:,:,[9,10,11]),3);
					romsMaster.std_2dplots(tmproms,tmpdata,obj.region.lon_rho,obj.region.lat_rho,options,tag,savepath,A.depth(d),A.plt);
				end	
			end
      		end % end method plot2DDiags

		%--------------------------------------------------------------------------------
      		function obj = make_zavg(obj,varargin)
        		% ------------------
			% - Vertically interpolates the sigma coordinates ROMS file to 
			% - WOA constant depth levels
			% ------------------	

			A.avg      = [];
			A.gridFile = [];
			A          = parse_pv_pairs(A,varargin);
			
			% - get initial directory
			od = pwd;

			% - create z_avg_dir if it doesnt exist
			z_avg_dir = [obj.paths.simPath,'z_avg'];
			if exist(z_avg_dir) ~= 7
				mkdir(z_avg_dir);
			end
			
			% - run z_slice_onwoa if avg contains 'h'
			try
				kill
				h = ncread(obj.paths.avg,'h');
				cd(z_avg_dir);
				cmd = ['z_slice_onwoa ', obj.paths.avg];
         			system(cmd);
			catch
				path1 = [obj.paths.simPath,'avg'];
				path2 = [obj.paths.simPath,'z_avg'];
				cmd = ['cp ',obj.paths.grid,' ',path1,'/.'];
				system(cmd);
				cd(path1);
				dps = obj.grid.woa0p25.depth;
				ss  = [];
				for i = 1:length(dps)
					if i < length(dps)
						str    = [num2str(dps(i)),' '];
					else
						str    = [num2str(dps(i))];
					end
					str = {str};
					ss     = strcat(ss,str);
				end
				[pathstr, gname, gext] = fileparts(obj.paths.grid);
				cmd = ['zslice ',ss{:},' ',[gname,gext],' ',A.avg]; 
				system(cmd);
				cmd = ['mv z_',A.avg,' ',path2];
				system(cmd)
			end

			% - return to original directory
			cd(od);

      		end % end method make_zavg
	end % end methods declarations
	%----------------------------------------------------------------------------------------

	%----------------------------------------------------------------------------------------
	% Utility functions (static methods)
	methods (Static)
		%--------------------------------------------------------------------------------
		function std_Zplots(tmproms,tmpdata,gridopt,options,savepath,tag,xlimit,zlimit,plttype);
			% -----------------------
			% - produces depth section plots along a prescriped parallel
			% -----------------------

			% - load colormap and set up axes positions
         		cmap = load('/data/project1/yangsi/MATLAB/functions/colormaps/cmap.mat');
			pos1 = [0.1100 0.7133 0.7527 0.2267];
			pos2 = [0.1100 0.3867 0.7527 0.2267];
			pos3 = [0.1100 0.0800 0.7527 0.2267];
			cbp1 = [0.8792 0.7132 0.0317 0.2272];
			cbp2 = [0.8792 0.3866 0.0317 0.2272];
			cbp3 = [0.8792 0.0799 0.0317 0.2272];

			% - get some other plot info
			indx = find(~nansum(tmproms')==0);
			indz = find(~nansum(tmproms)==0);
			if isempty(xlimit)
				xl   = gridopt.slicedeg(:,1); xl = xl(indx); xl = [min(xl) max(xl)];
			else
				xl   = xlimit;
			end
			if isempty(zlimit)
				zl   = gridopt.slicedepth(1,:); zl = zl(indz); zl = [min(zl)-100 0];
			else
				zl   = zlimit;
			end

			% - define anomalies
			tmpdiff = (tmproms - tmpdata);

			% - Automatically set caxis limits
			if isempty(xlimit) & isempty(zlimit)
				rclim    = prctile(tmproms(:),[1 99]);
				vclim    = prctile(tmpdata(:),[1 99]);
				rvclim   = [min([rclim(1) vclim(1)]) max([rclim(2) vclim(2)])];
				dclim    = prctile(tmpdiff(:),[1 99]);
				dclim      = [-max(abs(dclim)) max(abs(dclim))];
				[absCaxis ,absLevs ] = romsMaster.oom_levs(rvclim);
                                if sum(~isnan(dclim)) == 0
                                        diffCaxis = [-1 1];
                                        diffLevs  = [-1:0.5:1];
                                else
                                        [diffCaxis,diffLevs] = romsMaster.oom_levs(dclim,'diff',1);
                                end
			elseif isempty(xlimit) & ~isempty(zlimit)
				depth  = gridopt.slicedepth(1,:);
				zind   = find(zlimit(1) <= depth & depth <= zlimit(2));
				rdat   = tmproms(:,zind);
				vdat   = tmpdata(:,zind);
				ddat   = tmpdiff(:,zind);
				rclim  = prctile(rdat(:),[1 99]);
				vclim  = prctile(vdat(:),[1 99]);
				rvclim = [min([rclim(1) vclim(1)]) max([rclim(2) vclim(2)])];
				dclim  = prctile(ddat(:),[1 99]);
				dclim  = [-max(abs(dclim)) max(abs(dclim))];
				[absCaxis ,absLevs ] = romsMaster.oom_levs(rvclim);
				if sum(~isnan(dclim)) == 0
					diffCaxis = [-1 1];
					diffLevs  = [-1:0.5:1];
				else
					[diffCaxis,diffLevs] = romsMaster.oom_levs(dclim,'diff',1);
				end
			elseif ~isempty(xlimit) & ~isempty(zlimit);
				depth  = gridopt.slicedepth(1,:);
				deg    = gridopt.slicedeg(:,1);
				zind   = find(zlimit(1) <= depth & depth <= zlimit(2));
				xind   = find(xlimit(1) <= deg   & deg   <= xlimit(2));
				rdat   = tmproms(xind,zind);
				vdat   = tmpdata(xind,zind);
				ddat   = tmpdiff(xind,zind);
				rclim  = prctile(rdat(:),[1 99]);
				vclim  = prctile(vdat(:),[1 99]);
				rvclim = [min([rclim(1) vclim(1)]) max([rclim(2) vclim(2)])];
				dclim  = prctile(ddat(:),[1 99]);
				dclim  = [-max(abs(dclim)) max(abs(dclim))];
				[absCaxis ,absLevs ] = romsMaster.oom_levs(rvclim);
	                        if sum(~isnan(dclim)) == 0
                                        diffCaxis = [-1 1];
                                        diffLevs  = [-1:0.5:1];
                                else
                                        [diffCaxis,diffLevs] = romsMaster.oom_levs(dclim,'diff',1);
                                end
			end
			
			fig=figure;
         		set(fig, 'Units', 'centimeters', 'Position', [0, 0, 17.8,22], 'PaperUnits',...
            			 'centimeters', 'PaperSize', [17.8, 22]);
         		set(gcf, 'color','w','Visible','off');
         		ax1=axes('Position',pos1,'Units','normalized');
			    pcolor(gridopt.slicedeg,gridopt.slicedepth,tmproms);
			    shading flat
			    c1 = colorbar('Position',cbp1);
            		    tiVar = options.roms.var; tiSrc = 'ROMS'; tiTag = tag;
            		    title([tiSrc,':',' ',tiTag,' ',tiVar,' along ',num2str(gridopt.slicesect),'{^o} ',gridopt.slicecoord])
			    ylabel('Depth (m)');
			    if strcmp(gridopt.slicecoord,'latitude')
				xlabel('Longitude (^o)');
			    elseif strcmp(gridopt.slicecoord,'longitude');
				xlabel('Latitude (^o)');
			    end
			    ylabel(c1,options.units);
			    xlim(xl); ylim(zl);
			    caxis(absCaxis);
			    c1.XTick = absLevs;
			    c1.XTickLabel = absLevs;
         		ax2=axes('Position',pos2,'Units','normalized');
			    pcolor(gridopt.slicedeg,gridopt.slicedepth,tmpdata);
			    shading flat
			    c2 = colorbar('Position',cbp2);
			    tiVar = options.roms.var; tiSrc = options.data.varsource; tiTag = tag;
			    title([tiSrc,':',' ',tiTag,' ',tiVar,' along ',num2str(gridopt.slicesect),'{^o}' ,gridopt.slicecoord])
			    ylabel('Depth (m)');
			    if strcmp(gridopt.slicecoord,'latitude')
				xlabel('Longitude (^o)');
			    elseif strcmp(gridopt.slicecoord,'longitude');
				xlabel('Latitude (^o)');
			    end
			    ylabel(c2,options.units);
         		    caxis(absCaxis);
			    xlim(xl); ylim(zl);
			    c2.XTick = absLevs;
			    c2.XTickLabel = absLevs;
			ax3=axes('Position',pos3,'Units','normalized');
			    pcolor(gridopt.slicedeg,gridopt.slicedepth,(tmproms - tmpdata));
			    shading flat
			    c3 = colorbar('Position',cbp3);
			    title(['Difference']);
			    ylabel('Depth (m)');
			    if strcmp(gridopt.slicecoord,'latitude')
				xlabel('Longitude (^o)');
			    elseif strcmp(gridopt.slicecoord,'longitude');
				xlabel('Latitude (^o)');
			    end
			    ylabel(c3,options.units);
			    xlim(xl); ylim(zl);
			    colormap(ax3,cmap.cmap_ferretdiff);
			    caxis(diffCaxis);
			    c3.XTick = diffLevs;
			    c3.XTickLabel = diffLevs;
			if strcmp(gridopt.slicecoord,'latitude');
         			ssave = [savepath,tiTag,'-',tiVar,'_vs_obs-',num2str(gridopt.slicesect),'lat'];
			elseif strcmp(gridopt.slicecoord,'longitude');
				ssave = [savepath,tiTag,'-',tiVar,'_vs_obs-',num2str(gridopt.slicesect),'lon'];
			end
			if plttype == 1
				export_fig(fig,ssave,'-png'); close(fig);
			elseif plttype == 2
				pltshow; close(fig);
			end
		end % end static method std_Zplots

		%--------------------------------------------------------------------------------
		function std_2dplots(tmproms,tmpdata,tlon,tlat,options,tag,savepath,depth,plttype)
        		% ------------------
			% - produces horizontal 2D plots, given input options 
			% ------------------	
	
			% - get surface level
			if depth == 0
				tiDepth = ' @ SFC';
				fdepth  = '-sfc';
			elseif ~isnan(depth)
				tiDepth = [' @ ',num2str(depth),'m'];
				fdepth  = ['-',num2str(depth)];
			elseif isnan(depth)
				tiDepth = [];
				fdepth  = [];
			end

			% - define map limits and ticks
         		latbounds = [floor(nanmin(tlat(:))), floor(nanmax(tlat(:)))+1];
         		latticks  = (latbounds(1):floor(range(latbounds)/10):latbounds(2));
         		lonbounds = [floor(nanmin(tlon(:))), floor(nanmax(tlon(:)))+1];
         		lonticks  = (lonbounds(1):floor(range(lonbounds)/10):lonbounds(2));

			% - define anomalies, get yticks/caxis limits
         		tmpdiff = tmproms - tmpdata;
			
			% - get caxis limits automatically
         		if options.surf.logplt
            			tmproms 		= real(log10(tmproms));
            			tmproms(isnan(tmproms)) = 0;
            			tmpdata 		= real(log10(tmpdata));
            			tmpdata(isnan(tmpdata)) = 0;
            			tmpdiff 		= romsMaster.dfloglevs(tmpdiff,0.01);
            			tmpdiff(isnan(tmpdiff)) = 0;
            			absLevs			= log10(options.surf.absLevs);
            			diffLevs 		= romsMaster.dfloglevs(options.surf.diffLevs,0.01);
            			absCaxis 		= real(log10(options.surf.absCaxis));
            			diffCaxis 		= [diffLevs(1),diffLevs(end) ];
				% - reapply mask
				tmproms(options.mask==0) = NaN;
				tmpdata(options.mask==0) = NaN;
				tmpdiff(options.mask==0) = NaN;	
			%elseif depth == 0
			%	absCaxis  	     = options.surf.absCaxis;
			%	absLevs              = options.surf.absLevs;
			%	[absCaxis,absLevs]   = romsMaster.oom_levs(absCaxis);
            		%	diffCaxis 	     = options.surf.diffCaxis;
			%	diffLevs             = options.surf.diffLevs;
				%[diffCaxis,diffLevs] = romsMaster.oom_levs(diffCaxis,'diff',1);
			else
				rclim                = prctile(tmproms(:),[1 99]);
				vclim                = prctile(tmpdata(:),[1 99]);
				absCaxis             = [min([rclim(1) vclim(1)]) max([rclim(2) vclim(2)])];
				[absCaxis,absLevs]   = romsMaster.oom_levs(absCaxis);
				dclim                = prctile(tmpdiff(:),[1 99]);
				diffCaxis            = [-max(abs(dclim)) max(abs(dclim))];
				if ismember(1,isnan(diffCaxis)) % empty
					diffCaxis = [-1 1];
					diffLevs  = [-1 1];
				else
					[diffCaxis,diffLevs] = romsMaster.oom_levs(diffCaxis,'diff',1);
				end
			end
			
			% - load colormap and set up axes positions
         		cmap = load('/data/project1/yangsi/MATLAB/functions/colormaps/cmap.mat');
         		lmy  = 0.02; tmy = 0.04; dsy = (1 - 3*tmy-lmy)./3;
         		mx   = 0.02; dsx = 1 - 2*mx;
         		pos1 = [mx,lmy+2*tmy+2*dsy,dsx,dsy];
         		pos2 = [mx,lmy+tmy+dsy   ,dsx,dsy];
         		pos3 = [mx,lmy         ,dsx,dsy];
         		
			fig=figure;
         		set(fig, 'Units', 'centimeters', 'Position', [0, 0, 17.8,22], 'PaperUnits',...
            		         'centimeters', 'PaperSize', [17.8, 22]);
         		set(gcf,'color','w','Visible','off');
         		ax1=axes('Position',pos1,'Units','normalized');
            	    	    m_proj('mercator','lat',latbounds, 'lon', lonbounds);
            	    	    mintmp = nanmin(tmproms(:)); tmplevs = absLevs;
			    tmproms(tmproms>1e16) = NaN;
			    sg = m_pcolor(tlon,tlat,tmproms); shading flat
			    m_coast('patch',[0,0,0]+0.5,'edgecolor','k');
            	    	    m_grid('box','fancy','linestyle', 'none','xtick',lonticks,...
               	    	           'ytick',latticks,'backcolor',[.9 .99 1],'yticklabels',[],'xticklabels',[]);
            	    	    tiVar = options.roms.var; tiSrc = 'ROMS'; tiTag = tag;
            	    	    tiUnits = options.units;
            	    	    title([tiSrc,':',' ',tiTag,' ',tiVar,tiDepth],'Interpreter','none');
            	    	    set(ax1,'FontSize',10);
            	    	    c1=colorbar(ax1,'XTick', absLevs, 'XTickLabel', absLevs);
			    if options.surf.logplt
				c1.XTickLabel = options.surf.absLevs;
			    end
            	    	    caxis(ax1,absCaxis);
			    ylabel(c1,tiUnits);
         		ax2=axes('Position',pos2,'Units','normalized');
            	    	    m_proj('mercator','lat',latbounds, 'lon', lonbounds);
            	    	    mintmp = nanmin(tmpdata(:)); tmplevs = absLevs;
            	    	    sg2 = m_pcolor(tlon,tlat,tmpdata); shading flat
			    m_coast('patch',[0,0,0]+0.5,'edgecolor','k');
            	    	    m_grid('box','fancy','linestyle', 'none','xtick',lonticks,...
               	    	           'ytick',latticks,'backcolor',[.9 .99 1],'yticklabels',[],'xticklabels',[]);
            	    	    tiVar = options.roms.var; tiSrc = options.data.varsource; tiTag = tag;
            	    	    tiUnits = options.units;
            	    	    title([tiSrc,':',' ',tiTag,' ',tiVar,tiDepth],'Interpreter','none');
            	    	    set(ax2,'FontSize',10);
            	    	    c2=colorbar(ax2,'XTick', absLevs, 'XTickLabel', absLevs);
			    if options.surf.logplt
				c2.XTickLabel = options.surf.absLevs;
			    end
            	    	    caxis(ax2,absCaxis);
			    ylabel(c2,tiUnits);
         		ax3=axes('Position',pos3,'Units','normalized');
            	    	    m_proj('mercator','lat',latbounds, 'lon', lonbounds);
            	    	    mintmp = nanmin(tmpdiff(:)); tmplevs = diffLevs;
            	    	    sg3 = m_pcolor(tlon,tlat,tmpdiff); shading flat
			    m_coast('patch',[0,0,0]+0.5,'edgecolor','k');
            	    	    m_grid('box','fancy','linestyle', 'none','xtick',lonticks,...
               	    	           'ytick',latticks,'backcolor',[.9 .99 1],'yticklabels',[],'xticklabels',[]);
            	    	    tiUnits = options.units;
			    title('Difference'); 
			    set(ax1,'FontSize',10);
            	    	    c3=colorbar(ax3,'XTick', diffLevs,'XTickLabel', diffLevs);
			    if options.surf.logplt
				c3.XTickLabel = options.surf.diffLevs;
			    end
            	    	    caxis(ax3,diffCaxis);
			    ylabel(c3,tiUnits);
            	    	    colormap(ax3,cmap.cmap_ferretdiff);
			if ~isnan(depth)
         			ssave = [savepath,tiTag,'-',tiVar,'_vs_obs-',options.data.var,fdepth];
			else
				ssave = [savepath,tiTag,'-',tiVar,'_vs_obs-',options.data.var];
			end
			if plttype == 1
				export_fig(fig,ssave,'-png'); close(fig);
			elseif plttype == 2
				pltshow; close(fig);
			end
		end % end static method std_2dplots
      
		%--------------------------------------------------------------------------------
		function [smean sstd scount sprob]  = scatter_density(xdata,ydata,xbounds,ybounds)
        		% ------------------
			% - grid data for scatter density plots  
			% ------------------
         		
			xbounds = xbounds(:); 
			ybounds = ybounds(:);
         		scount  = nan(length(xbounds)-1,length(ybounds)-1);
         		sprob   = scount;
         		smean   = nan(length(xbounds)-1,1);
         		sstd    = smean;
         
         		for i = 1 : length(xdata(:))
            			if mod(i,10000) == 0
               				i
            			end
            			indx = find(xdata(i)>=xbounds(1:end-1) & xdata(i)<xbounds(2:end));
            			indy = find(ydata(i)>=ybounds(1:end-1) & ydata(i)<ybounds(2:end));
            			if isempty(indx)
               				if xdata(i) <= xbounds(1)
                  				indx=1;
               				elseif xdata(i)>=xbounds(end)
                  				indx=length(xbounds)-1;
               				end
            			end
         
            			if isempty(indy)
               				if ydata(i) <= ybounds(1)
                  				indy=1;
               				elseif ydata(i)>=ybounds(end)
                  				indy=length(ybounds)-1;
               				end
            			end
            			if isnan(scount(indx,indy))
                			scount(indx,indy)=1;
            			else
                			scount(indx,indy)=scount(indx,indy)+1;
            			end
            			if isnan(smean(indx))
                			smean(indx) = ydata(i);
            			else
                			smean(indx) = smean(indx) + ydata(i);
            			end
         		end
         
         		smean = smean ./nansum(scount,2);
         		sprob = scount./repmat(nansum(scount,2),1,length(ybounds)-1)*100;
		end % end static method scatter_density

		%--------------------------------------------------------------------------------
		function x = lon360(x)
        		% ------------------
			% - corrects negative longitudes such that all lon are between 0:360
			% ------------------
			
			x(x<0) = x(x<0) + 360;
      		end % end static method lon360

		%--------------------------------------------------------------------------------
      		function x = lon180(x)
        		% ------------------
			% - corrects longitudes such that all lon are between -180:180  
			% ------------------
         		
			x(x>180) = x(x>180) - 360;
      		end % end static method lon180

		%--------------------------------------------------------------------------------
      		function o2_corr = woao2corr(o2)
        		% ------------------
			% -  applies Bianchi 2012 correction to O2
			% ------------------
         		
			o2_corr=1.009*o2-2.523;
         		o2_corr(o2_corr<0)=0;
      		end % end static method woao2corr

		%--------------------------------------------------------------------------------
      		function x = struct2double(x)
        		% ------------------
			% - converts structure fields to double  
			% ------------------
         		
			ffields = fields(x);
         		for ff = 1 : length(ffields)
            			if isa(x.(ffields{ff}),'double')
               				x.(ffields{ff}) = double(x.(ffields{ff}));
            			end
         		end
      		end % end static method struct2double

		%--------------------------------------------------------------------------------
      		function loglevs=dfloglevs(levs,logminparam)
        		% ------------------
			% - converts positive levels into symmetrical centered around 0 levels
			% - in log space. For use in difference plots with log scales  
			% ------------------
         		
			levs(abs(levs)<logminparam)=sign(levs(abs(levs)<logminparam)).*logminparam;
         		loglevs=levs;
         		loglevs(levs>0)=+abs(log10(logminparam) - (log10(levs(levs>0))));
         		loglevs(levs<0)=-abs(log10(logminparam) - log10(abs(levs(levs<0))));
      		end % end static method dfloglevs

		%--------------------------------------------------------------------------------
      		function x = struct2single(x)
        		% ------------------
			% - converts structure fields to single  
			% ------------------
         		
			ffields = fields(x);
         		for ff = 1 : length(ffields)
            			if isa(x.(ffields{ff}),'double')
               				x.(ffields{ff}) = single(x.(ffields{ff}));
            			end
         		end
		end % end static method struct2single

		%--------------------------------------------------------------------------------
		function [mnmx,levs] = oom_levs(mnmx,varargin);
			% ----------------
			% - Gets order of magnitude estimate and sets difference levels for axis labels
			% - Essentially, used to automatically set axes such that there are 10 levels
			% - and the levels are nice rounded values rather than 16.666666 etc
			% ----------------
		
			A.diff = 0; % default to not centered around 0
			A = parse_pv_pairs(A,varargin);
	
			% Adjust min and max levels
			di     = diff(mnmx);
			di_oom = floor(log10(di));
			if di_oom >= 0
				mnmx = round(mnmx,-di_oom+1);
			else
				mnmx = round(mnmx,-di_oom);
			end 
			if A.diff == 1		
				mnmx = [-max(abs(mnmx)) max(abs(mnmx))];
			end

			% Get levels
			dabs = [diff(mnmx)/9];
			oom  = floor(log10(dabs));
			doom = [dabs/10^(oom)];
			doom = round(doom);
			dabs = doom*10^(oom);
			levs = [mnmx(1):dabs:mnmx(2)];
			if mnmx(2) > levs(end)
				levs = [levs levs(end)+dabs];
			end
			if A.diff == 0
				mnmx = [levs(1) levs(end)];
			elseif A.diff == 1
				tmplevs = [0:dabs:mnmx(2)];
				if mnmx(2) > tmplevs(end)
					tmplevs = [tmplevs tmplevs(end)+dabs];
				end
				levs = unique([-tmplevs tmplevs]);
				mnmx = [levs(1) levs(end)];
			end	
		end % end static method oom_levs

		%--------------------------------------------------------------------------------
		function [ws,wsc] = WindStress(u,v,lon,lat,ang)
        		% ------------------
			% - calculates wind stress and wind stress curl from zonal/meridional wind stress 
			% ------------------
			
			% - calculate wind stress
			ws = sqrt(u.*u + v.*v);

			% - check that lon/lat are gridded
			[a,b] = size(lon);
			if a == 1 | b == 1
				[lon,lat] = meshgrid(lon,lat);
			end
		
			% - check for yearly or monthly file
			[a,b,c] = size(u);

			% - track progress
			d        = a*b*c;
			progress = [d/10:d/10:d];
			cnt      = 0;
			pcnt     = 10;
			% - calculate wind stress curl
			E   = cos(ang).*u - sin(ang).*v;
			N   = sin(ang).*u + cos(ang).*v;
			wsc = ws.*NaN;
			for i = 2:size(lon,1)-1;
				for j = 2:size(lon,2)-1;
					dx = sw_dist([lat(i+1,j) lat(i-1,j)],...
						     [lon(i+1,j) lon(i-1,j)],'km')*1000;
					dy = sw_dist([lat(i,j+1) lat(i,j-1)],...
						     [lon(i,j+1) lon(i,j-1)],'km')*1000;
					for t = 1:c
						wsc(i,j,t) = (N(i+1,j,t)-N(i-1,j,t))/dx - ...
							     (E(i,j+1,t)-E(i,j-1,t))/dy ;
						cnt = cnt + 1;
						if ismember(cnt,progress)
							disp([num2str(pcnt),'% complete']);
							pcnt = pcnt + 10;
						end
					end
				end
			end
      		end % end static method WindStress

		%--------------------------------------------------------------------------------
		function [lon,lat] = glodap_coord
			% -------------------
			% - Function to grab allowable coordinates for GLODAPv2 comparisons
			% -------------------

			% Load GLODAPv2.2020 data
			load('/data/project1/demccoy/ROMS/validation/GLODAPv2/glodap_v2_2020_format.mat');
			tmplon  = [glodap.lon]; tmplon(tmplon<0) = tmplon(tmplon<0)+360;
			tmplat  = [glodap.lat];
			tmptime = [glodap.time];

			% Plot map of locations
			figure
			set(gcf,'Position',[451          75        1376         946],'color','w');
			plot(tmplon,tmplat,'.k'); hold on; plot_coast; xlim([0 360]); ylim([-90 90]);
			set(gca,'XTick',[0:10:360],'YTick',[-90:10:90]);
			grid on

			% Get user input of lon/lat
			disp('Choose location for 5x5 grid box');
			for i = 1:10
				if i > 1
					disp('Try again');
					delete(r);
				end

				% Query lon/lat box
				x = []; y = []; idx = [];
                       		[x] = input('---------------------------\nLon/lat estimate ([lon,lat]):\n---------------------------\n>> ');
				[y] = x(2); 
				[x] = x(1);

				% Find number of profiles within 5x5 box
				idx   = find(x-2.5 <= tmplon & tmplon <= x+2.5 & y-2.5 <= tmplat & tmplat <= y+2.5);
				r     = rectangle('Position',[x-2.5 y-2.5 5 5],'FaceColor','b');
				
				% Query keep or reject
				disp([num2str(length(idx)),' profiles found']);
                       		q = input('---------------------------\nOK? (1 == YES):\n---------------------------\n>> ');
				if q == 1 & ~isempty(idx)
					break
				else
					continue
				end
			end

			% Save results
			lon = x;
			lat = y;

		end % end static method glodap_coord
	end % end static methods declarations
	%----------------------------------------------------------------------------------------
end % end class
%------------------------------------------------------------------------------------------------
