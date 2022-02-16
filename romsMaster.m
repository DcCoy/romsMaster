%------------------------------------------------------------------------------------------------
classdef romsMaster
%-------------------
% Class used to load ROMS data and make various plots,
% comparison figures against validation data, and many
% other actions.
%
% To begin, simply initialize romsMaster using 'sim' from current_sims.mat
% - obj = init(romsMaster,sim)
%
% To view available routines
% - methods(obj)
%
% For help
% - help obj.(method)
%------------------
	%----------------------------------------------------------------------------------------
	properties
		info     % contains variable information
		paths    % paths to data, directories, diagnostic data, etc
		grid     % contains grid data and extra fields
		region   % contains grid data in defined region
		slice    % contains slice coordinates
		profile  % contains profile coordinates
		romsData % contains ROMS data for comparisons
		diagData % contains validation data for comparions
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
			%            'manual' will not auto-run defineRegion (do it after calling init(obj,...))
			%            'default' will use default region in current_sims.mat
			% - year   = override the default year (set in current_sims.m)
			%
			% Example:
			% - obj = init(obj,'peru')         <-- if obj defined
			% - obj = init(romsMaster,'peru')  <-- if obj undefined
			% -------------------
			
			%  Begin
			disp(' ');
			disp('---------------------------------');
			disp('---------------------------------');
			disp('            romsMaster           ');
			disp('---------------------------------');
			disp('---------------------------------');
			disp(' ');               	

			%  Input test
			if nargin < 2
				% Load default settings
				sim = 'default';	
				A.region = 'full';
				A.year   = [];
			elseif nargin < 3
				% Load full domain
				A.region = 'full';
				A.year   = [];
			else
				% Optional inputs...
				A.region = ['default'];
				A.year   = [];
				A        = parse_pv_pairs(A,varargin);
			end
 
			% Suppress warning messages and addpath to Danny's scripts
			warning off
			addpath /data/project1/demccoy/matlab_scripts/
			addpath /data/project1/demccoy/ROMS/tools/Roms_tools_MF/Preprocessing_tools/
			
			% Run current_sims
			current_sims;
			
			% Load simulation
			load('current_sims.mat',sim);	
			eval(['unpactStruct(',sim,');']);

			% override year?
			if ~isempty(A.year)
				avgfile = ['avg_',num2str(A.year),'.nc'];
				ext     = [num2str(A.year)];
				year    = A.year;
			end
			
			% grab paths according to inputs
			obj.paths.simPath = simPath; % root path

			% grab file paths
			obj.paths.avg    = [obj.paths.simPath,'avg/',avgfile];
			obj.paths.his    = [obj.paths.simPath,'his/','his_',ext,'.nc'];
			obj.paths.flux   = [obj.paths.simPath,'phys_flux/','phys_flux_avg_',ext,'.nc'];
			obj.paths.zavg 	 = [obj.paths.simPath,'z_avg/','z_',avgfile];
			obj.paths.grid   = gfile;

			% initiate directories if they dont exist
			[pathstr, name, ext] = fileparts(obj.paths.avg);
			mkdir([obj.paths.simPath,'Figures']);
			mkdir([obj.paths.simPath,'Figures/',name,'/']);
			mkdir([obj.paths.simPath,'Figures/',name,'/Diagnostic']);
			mkdir([obj.paths.simPath,'Figures/',name,'/Comparison']);
		
			% grab plot paths
			obj.paths.plots.diag	= [obj.paths.simPath,'Figures/',name,'/Diagnostic/'];
			obj.paths.plots.comp	= [obj.paths.simPath,'Figures/',name,'/Comparison/'];
			obj.paths.plots.tmpfigs = ['/data/project1/demccoy/tmpfigs/'];

			% Get info for ROMS variables
			obj = romsInfo(obj,r_tag,year);
			
			% Load grid
			obj = loadGrid(obj);

			% make z_avg file if it doesn't exist
			if exist([simPath,'z_avg/z_',avgfile]) ~= 2
				disp('---------------------------------');
				disp('Create z_avg Climatology File');
				disp('---------------------------------');
				obj = make_zavg(obj,avgfile);
			end

			% Save subregion indices
			obj.region.lat_lim   = rlati;
			obj.region.lon_lim   = rloni;
			obj.region.dep_lim   = rdepi;
			obj.region.coast_lim = rcoast;

			% Grab subregion
			if strcmp(A.region,'full')
				lnl = [1 obj.grid.nx-1]; %-1 for x-advection
				ltl = [1 obj.grid.ny-1]; %-1 for y-advection
				cst = [-inf];
				dpl = [-inf inf];
				obj = defineRegion(obj,'lon_lim',[lnl],'lat_lim',[ltl],'dep_lim',[dpl],'coast_lim',[cst]);
			elseif strcmp(A.region,'manual');
				% Do nothing, user will define region
			elseif strcmp(A.region,'default');
				obj = defineRegion(obj);
			end

			% Get paths for diagnostic products
			obj = initDiag(obj);
		end % end methods init

		%--------------------------------------------------------------------------------
		function [obj] = loadGrid(obj)
			% ----------------------
			% Loads grid information into obj.grid
			%
			% NOTE:
			% Also loads z_r, z_w, Hz, so obj.paths.avg needs to be called (see init)
			% ----------------------

			rmpath('/data/project1/demccoy/ROMS/ROMS_tools/nc_tools/');
			addpath('/usr/local/MATLAB/R2019b/toolbox/matlab/imagesci/');

			% get grid coordinates
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

			% get more coord diags 
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

			% Load vertical coordinates
			tmpinfo = ncinfo(obj.paths.avg);
			ind = strmatch('theta_s',{tmpinfo.Attributes.Name});
			theta_s = tmpinfo.Attributes(ind).Value;
			ind = strmatch('theta_b',{tmpinfo.Attributes.Name});
			theta_b = tmpinfo.Attributes(ind).Value;
			ind = strmatch('hc',{tmpinfo.Attributes.Name});
			hc = tmpinfo.Attributes(ind).Value;
			try
				obj.grid.Hz       = double(ncread(obj.paths.avg,'Hz'));
				obj.grid.z_r      = double(ncread(obj.paths.avg,'z_r'));
				obj.grid.z_w      = double(ncread(obj.paths.avg,'z_w'));
			catch
				tmpzeta = ncread(obj.paths.avg,'zeta');
				for i = 1:obj.grid.nt
					tmpzr = zlevs4(obj.grid.h,tmpzeta(:,:,i),theta_s,theta_b,hc,obj.grid.nz,'r','new2012');
					tmpzw = zlevs4(obj.grid.h,tmpzeta(:,:,i),theta_s,theta_b,hc,obj.grid.nz,'w','new2012');
					obj.grid.z_r(:,:,:,i) = permute(tmpzr,[2 3 1]);
					obj.grid.z_w(:,:,:,i) = permute(tmpzw,[2 3 1]);
				end
				obj.grid.Hz = diff(obj.grid.z_w,1,3);
			end
			% Fix masks
			obj.grid.mask_rho(obj.grid.mask_rho==0) = NaN;

			% Get grid polygon
			obj.grid.polygon(:,1) = [obj.grid.lon_rho(1,:) obj.grid.lon_rho(:,end)'...
									 fliplr(obj.grid.lon_rho(end,:)) flipud(obj.grid.lon_rho(:,1))'];
			obj.grid.polygon(:,2) = [obj.grid.lat_rho(1,:) obj.grid.lat_rho(:,end)'...
									 fliplr(obj.grid.lat_rho(end,:)) flipud(obj.grid.lat_rho(:,1))'];

			% Make coord data all single
			obj.grid = romsMaster.struct2double(obj.grid);
		end % end methods loadGrid

		%--------------------------------------------------------------------------------
		function [obj] = clearROMS(obj) 
			% ----------------------
			% Clears loaded data and coordinates
			%
			% Usage:
			% - [obj] = clearROMS(obj) 
			% ----------------------
			
			% Clear fields
			obj.slice = [];
			obj.profile = [];
			obj.romsData = [];
			obj.diagData = [];

		end % end methods clearROMS
		%--------------------------------------------------------------------------------
		function [obj] = romsInfo(obj,r_tag,year)
			% ----------------------
			% Obtains info from roms file
			%
			% Usage:
			% - [obj] = romsInfo(obj,r_tag,year) 
			%
			% Inputs:
			% - r_tag = name of simulation (i.e. peru_chile_0p1)
			% - year  = simulation year (i.e. 2000)
			% ----------------------

			rmpath /data/project1/demccoy/ROMS/ROMS_tools/nc_tools/

			% list variables
			obj.info = ncinfo(obj.paths.avg);
			cnt2d = 1;
			cnt3d = 1;
			nt = obj.info.Dimensions(find(strcmp('time',{obj.info.Dimensions.Name})==1)).Length;
			for i = 1:length(obj.info.Variables)
				if length(obj.info.Variables(i).Size)==3 & obj.info.Variables(i).Size(3) == nt;
					obj.info.var2d{cnt2d}  = obj.info.Variables(i).Name;
					for j = 1:length(obj.info.Variables(i).Attributes);
						if strcmp(obj.info.Variables(i).Attributes(j).Name,'long_name');
							nameidx = j;
						elseif strcmp(obj.info.Variables(i).Attributes(j).Name,'units');
							unitidx = j;
						end
					end
					obj.info.name2d{cnt2d} = obj.info.Variables(i).Attributes(nameidx).Value;
					obj.info.unit2d{cnt2d} = obj.info.Variables(i).Attributes(unitidx).Value;
					obj.info.idx2d(cnt2d)  = i;
					cnt2d                  = cnt2d + 1;
				elseif length(obj.info.Variables(i).Size)==4
					obj.info.var3d{cnt3d} = obj.info.Variables(i).Name;
					for j = 1:length(obj.info.Variables(i).Attributes);
						if strcmp(obj.info.Variables(i).Attributes(j).Name,'long_name');
							nameidx = j;
						elseif strcmp(obj.info.Variables(i).Attributes(j).Name,'units');
							unitidx = j;
						end
					end
					obj.info.name3d{cnt3d} = obj.info.Variables(i).Attributes(nameidx).Value;
					obj.info.unit3d{cnt3d} = obj.info.Variables(i).Attributes(unitidx).Value;
					obj.info.idx3d(cnt3d)  = i;
					cnt3d                  = cnt3d + 1;
				end
			end

			% Add some other info
			obj.info.r_tag = r_tag;

			% Grab time info from average file
			fieldnames = {obj.info.Attributes.Name};
			fieldvalue = {obj.info.Attributes.Value};
			dt         = fieldvalue(strcmp(fieldnames,'dt'));
			dt         = double(dt{1});
			navg       = fieldvalue(strcmp(fieldnames,'navg'));
			navg       = double(navg{1});
			ntimes     = fieldvalue(strcmp(fieldnames,'ntimes'));
			ntimes     = double(ntimes{1});

			% Compute dt according to output frequency
			obj.info.Year   = year;
			obj.info.Freq   = dt*navg/86400;
			obj.info.Ntimes = ((dt*ntimes)/86400)/(obj.info.Freq);
			if obj.info.Freq > 27 & obj.info.Freq < 32
				obj.info.time_string     = ['Monthly'];
				obj.info.time_string_idv = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
			elseif obj.info.Freq == 1
				obj.info.time_string = ['Daily'];
				obj.info.time_string_idv = num2cell(1:31);
			elseif obj.info.Freq < 1
				obj.info.time_string = ['Hourly'];
				obj.info.time_string_idv = num2cell(1:24);
			end

			% Replace Unit Strings with Latex version
			strings_to_replace = {'meter second-1','Celsius','PSU','mMol P','mMol O2','mMol N','mMol',...
								  'mMol C','mg Chl-a','mMol CaCO3','meter2 second-1','W m-2','mmol C/m3/s',...
								  'mmol/m3/s','mg Chl/m3','mmol N/m3/s','mmol N2O/m3/s','meter',...
								  'ppm','mmol/m2/s','mmol/m3'};

			replacement_string = {'$m$ $s^{-1}$','$^{o}C$','$PSU$','$mmol$ $P$','$mmol$ $O_2$','$mmol$ $N$','$mmol$',...
								  '$mmol$ $C$','$mg$ $Chl-A$','$mmol$ $CaCO_3$','$m^{2}$ $s^{-1}$','$W$ $m^{-2}$','$mmol$ $C$ $m^{-3}$ $s^{-1}$',...
								  '$mmol$ $m^{-3}$ $s^{-1}$','$mg$ $C$ $m^{-3}$','$mmol$ $N$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$m$',...
								  '$ppm$','$mmol$ $m^{-2}$ $s^{-1}$','$mmol$ $m^{-3}$'};
			for i = 1:length(strings_to_replace)
				obj.info.unit3d = strrep(obj.info.unit3d,strings_to_replace{i},replacement_string{i});
				obj.info.unit2d = strrep(obj.info.unit2d,strings_to_replace{i},replacement_string{i});
			end
		end % end methods romsInfo

		%--------------------------------------------------------------------------------
		function [fig,ax] = gridView(obj,varargin)
			% ----------------------
			% Plots the lat/lon indices of a grid file
			%
			% Usage:
			% - [fig,ax] = gridView(obj,varargin)
			%
			% Inputs (varargin):
			% - dx    = plot lon lines separated by dx (default = 20)
			% - dy    = plot lat lines separated by dy (default = 20)
			% - ticks = 0 (no lon/lat labels), 1 (yes), 2 (fancy box)
			% - font  = tick font size (default = 12) 
			% - save  = 1 (print), 0 (no print)
			%
			% Example:
			% - [fig,ax] = gridView(obj,'dx',20,'dy',20)
			% ----------------------

			% Grab inputs (varargin)
			A.dx    = [20];
			A.dy    = [20];
			A.ticks = [1];
			A.font  = [12];
			A.save  = [0];
			A       = parse_pv_pairs(A,varargin);

			% Plot lon/lat lines
			fig(1) = piofigs('lfig',1.5);
			ax(1)  = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks,'font',A.font);	
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
			fname = [obj.info.r_tag,'_grid'];
			if A.save == 1
				export_fig('-jpg',[obj.paths.simPath,fname]);
			end
			if nargout < 1
				pltjpg(1);
				close(fig(1));
			end
		end % end method gridView
		
		%--------------------------------------------------------------------------------
		function [fig,ax] = regionView(obj,varargin)
			% ----------------------
			% Plots a ROMS grid map along with a regional grid
			%
			% Usage:
			% - [fig,ax] = regionView(obj,varargin)
			%
			% Inputs (varargin):
			% - ticks = 0 (no lon/lat labels), 1 (yes), 2 (fancy box)
			% - font  = tick font size (default = 12) 
			% - save  = 1 (print), 0 (no print)
			%
			% Example:
			% - [fig,ax] = regionView(obj);
			% ----------------------

		% Grab inputs (varargin)
		A.ticks = [1];
		A.font  = [12];
		A.save  = [0];
		A       = parse_pv_pairs(A,varargin);
	
			% Optional plot of region
			if ~isempty(obj.region.lon_rho)
				% Generate whole map
				fig(2)  = piofigs('lfig',1.5);
				set(0,'CurrentFigure',fig(2));
				[ax(2)] = map_plot(fig(2),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks,'font',A.font);
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
				fname = [obj.info.r_tag,'_region'];
				if A.save == 1
					export_fig('-jpg',[obj.paths.simPath,fname]);
				else
					pltjpg(1);
				end
				if nargout < 1
					pltjpg(1);
					close(fig(1));
				end
			end
		end % end method regionView

		%--------------------------------------------------------------------------------
		function obj = Dist2Coast(obj)
			% --------------------
			% Call this function to calculate each grid cell's distance-to-coast
			% --------------------

			% Get grid info
			lon  = obj.grid.lon_rho;
			lat  = obj.grid.lat_rho;
			mask = obj.grid.mask_rho;
			[Mp, Lp] = size(obj.grid.lon_rho);

			% Conversion to radians
			d2r = pi/180;
			cdist = 0*lon + 1e10;
			nch(1) = 1;
			nch(2) = 3;
			nch(3) = 4;
			disp(' Getting distance from coast');
			disp(' Starting the looper')
			for it = 1:1
				ncx = nch(it);
				ncy = nch(it);
				for ic = 1:ncx
					for jc = 1:ncy
						[ic jc];
						i0 = 1 + (ic-1)*ceil(Lp/ncx);
						i1 = i0+ceil(Lp/ncx) + 20;
						j0 = 1 + (jc-1)*ceil(Mp/ncy);
						j1 = j0+ceil(Mp/ncy) + 20;
						i1 = min(i1,Lp); j1 = min(j1,Mp);
						lons = lon(j0:j1,i0:i1);
						lats = lat(j0:j1,i0:i1);
						masks= mask(j0:j1,i0:i1);
						lab = 0*masks + 1; 
						lab(2:end-1,2:end-1) = masks(1:end-2,2:end-1)+...
											   masks(3:end,2:end-1)+...
											   masks(2:end-1,1:end-2)+...
											   masks(2:end-1,3:end);
						mlon = lons;mlon(masks>0|lab<1) = [];
						mlat = lats;mlat(masks>0|lab<1) = [];
						if mlon
							for j = j0:j1
								[j j1];
								for i = i0:i1
									if mask(j,i) < 1
										cdist(j,i) = 0;
									else
										dist = romsMaster.gc_dist(lon(j,i)*d2r,lat(j,i)*d2r,mlon*d2r,mlat*d2r);
										mdist = min(min(dist));
										cdist(j,i) = min(mdist,cdist(j,i));
									end
								end
							end
						else
						end
					end
				end
			end

			% Save output
			obj.grid.coastdist = cdist;
		end % end method Dist2Coast

		%--------------------------------------------------------------------------------
		function obj = changeInputs(obj,year,varargin)
			% ----------------------
			% Easier way to change the input files
			% 
			% Usage:
			% - obj = changeInputs(obj,year)
			%
			% Inputs:
			% - year = year, as an integer
			%
			% Options:
			% - ext = optional extension, as a string, for custom runs
			%
			% Example:
			% - obj = changeInputs(obj,year,'ext','_PCref');
			% ----------------------
			
			% Process inputs		
			A.ext = [''];
			A = parse_pv_pairs(A,varargin);

			% Change filenames
			avgfile  = ['avg/avg_',num2str(year),A.ext,'.nc'];
			hisfile  = ['his/his_',num2str(year),A.ext,'.nc'];
			fluxfile = ['phys_flux/phys_flux_avg_',num2str(year),A.ext,'.nc'];
			zavgfile = ['z_avg/z_avg_',num2str(year),A.ext,'.nc'];
			% Change filepaths
			obj.paths.avg  = [obj.paths.simPath,avgfile];
			obj.paths.his  = [obj.paths.simPath,hisfile];
			obj.paths.flux = [obj.paths.simPath,fluxfile];
			obj.paths.zavg = [obj.paths.simPath,zavgfile];
			% Update info
			obj = romsInfo(obj,obj.info.r_tag,year);
			% Redefine grid (mostly for z_r, z_w, Hz)
			obj = loadGrid(obj);
			% Redefine region
			obj.region.z_r        = obj.grid.z_r(obj.region.lon_lim(1):obj.region.lon_lim(2),...
									obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
			obj.region.z_w        = obj.grid.z_w(obj.region.lon_lim(1):obj.region.lon_lim(2),...
									obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
			obj.region.Hz         = obj.grid.Hz(obj.region.lon_lim(1):obj.region.lon_lim(2),...
									obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);	
			% Initiate directories if they dont exist
			[pathstr, name, ext] = fileparts(obj.paths.avg);
			mkdir([obj.paths.simPath,'Figures']);
			mkdir([obj.paths.simPath,'Figures/',name,'/']);
			mkdir([obj.paths.simPath,'Figures/',name,'/Diagnostic']);

			% grab plot paths
			obj.paths.plots.diag    = [obj.paths.simPath,'Figures/',name,'/Diagnostic/'];
			obj.paths.plots.tmpfigs = ['/data/project1/demccoy/tmpfigs/'];
		end % end methods changeInputs

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
			% - lon_lim   = x-grid indices (i.e. [200 400]) 
			% - lat_lim   = y-grid indices (i.e. [200 400])
			% - dep_lim   = z-grid depth limits (i.e. [-600 0]) for 0 - 600m
			% - coast_lim = minimum distance-to-coast, in km (see Dist2Coast)
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
			A.lon_lim   = [];
			A.lat_lim   = [];
			A.dep_lim   = [];
			A.coast_lim = [];
			A           = parse_pv_pairs(A,varargin);

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
			if ~isempty(A.coast_lim)
				obj.region.coast_lim = [A.coast_lim];
			elseif isempty(obj.region.coast_lim)
				obj.region.coast_lim = [-inf];
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
			catch; disp('Something went wrong with defineRegion...'); return; end
			obj.region.dx		  = 1./obj.region.pm;
			obj.region.dy		  = 1./obj.region.pn;
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
			
			% (OPTIONAL)
			if obj.region.coast_lim>(-inf) % coast_lim is defined in current_sims.m
				% Calculate distance-from-coast (m)
				disp('Refining mask via Dist2Coast');
				obj = Dist2Coast(obj);
				% Reduce to region
				obj.region.coastdist = obj.grid.coastdist(obj.region.lon_lim(1):obj.region.lon_lim(2),...
														  obj.region.lat_lim(1):obj.region.lat_lim(2));
				% Update mask
				obj.region.mask_rho(obj.region.coastdist < obj.region.coast_lim*1000) = NaN;
			end

			% Make 3D grid_area, grid_mask with depth limits applied
			for z = 1:obj.grid.nz
				% Grab temporary mask/area
				tmp_mask = obj.region.mask_rho;
				tmp_area = obj.region.grid_area;
				for t = 1:obj.grid.nt
					% Apply dep limits?
					if isempty(obj.region.dep_lim)
						obj.region.mask_rho3d(:,:,z,t) = obj.region.mask_rho;
						obj.region.area3d(:,:,z,t)     = obj.region.grid_area;
					else
						% Depth limits toggled	
						% Go through each time-entry, if any grid is above/below dep_lim, NaN it
						% This will build a mask that will be applied to all times equally 
						tmp_mask(obj.region.z_r(:,:,z,t)<obj.region.dep_lim(1)) = NaN;
						tmp_mask(obj.region.z_r(:,:,z,t)>obj.region.dep_lim(2)) = NaN;
						tmp_area(obj.region.z_r(:,:,z,t)<obj.region.dep_lim(1)) = NaN;
						tmp_area(obj.region.z_r(:,:,z,t)>obj.region.dep_lim(2)) = NaN;
					end
				end
				for t = 1:obj.grid.nt
					% Apply mask now	
					obj.region.mask_rho3d(:,:,z,t) = tmp_mask;
					obj.region.area3d(:,:,z,t)     = tmp_area;
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

			% Get region polygon
			obj.region.polygon(:,1) = [obj.region.lon_rho(1,:) obj.region.lon_rho(:,end)' ...
				fliplr(obj.region.lon_rho(end,:)) flipud(obj.region.lon_rho(:,1))'];
			obj.region.polygon(:,2) = [obj.region.lat_rho(1,:) obj.region.lat_rho(:,end)' ...
				fliplr(obj.region.lat_rho(end,:)) flipud(obj.region.lat_rho(:,1))'];
		end % end method defineRegion

		%--------------------------------------------------------------------------------
		function obj = getBudg(obj,varname,varargin)
			% --------------------
			% Main method to perform budget analysis on a ROMS tracer (varname).
			%
			% Usage:
			% - obj = getBudg(obj,varname,varargin)
			%
			% Required Inputs:
			% - varname = 'NO2','NO3','NH4','N2O','N2O_decomp','N2' for different budgets
			%
			% Optional Inputs:
			% - year = option to choose a different file
			%
			% Example:
			% - obj = getBudg(obj,'N2O')
			% This will instrcut budget methods (i.e. intRates, computeDzDt) on which
			% tracer to perform budget on. Here, tracking N2O.
			%
			% Available tracer budgets:
			% NO3
			% NO2
			% N2O
			% N2O_decomp
			% N2
			% --------------------

			disp('---------------------------------');
			disp('Get parameters');

			% Toggles
			A.year     = [];
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

			% Change year/file?
			if ~isempty(A.year)
				obj = changeInputs(obj,A.year);
			end

			% Get fields for budget calc
			if strcmp(varname,'NO3')
				vars   = {'NO3'}; 
				tits   = {'$NO_3$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'NITROX','DENITRIF1','SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE'};
				smseq  = [(1) (-1) (-1) (-1) (-1)];
				fluxes = {'SED_DENITRIF'};
				lvls   = {'sed'};
			elseif strcmp(varname,'NO2')
				vars   = {'NO2'};
				tits   = {'$NO_2$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'AMMOX','N2OAMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2',...
						  'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
				smseq  = [(1) (-2) (-1) (-1) (1) (-1) (-1) (-1) (-1)]; 
				fluxes = {[]};
				lvls   = {[]};
			elseif strcmp(varname,'NH4')
				vars   = {'NH4'};
				rates  = {[]}; 
				fluxes = {[]};
			elseif strcmp(varname,'N2')
				vars   = {'N2','N2_SED'};
				tits   = {'$N_2$','$N_2_{sed}$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'DENITRIF3','ANAMMOX'};
				smseq  = [(1) (1)];
				fluxes = {'FG_N2','SED_DENITRIF'};
				lvls   = {  'sfc',         'sed'};
			elseif strcmp(varname,'N2O')
				vars   = {'N2O'};
				tits   = {'$N_2O$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
				smseq  = [(0.5) (1) (-1)];
				fluxes = {'FG_N2O'};
				lvls   = {   'sfc'};
			elseif strcmp(varname,'N2O_SODEN')
				vars   = {'N2O_SODEN'};
				tits   = {'$N_2O_{den}$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'DENITRIF2','N2OSODEN_CONS'};
				smseq  = [(0.5) (-1)];
				fluxes = {'FG_N2O_SODEN'};
				lvls   = {   'sfc'};
			elseif strcmp(varname,'N2O_AO1')
				vars   = {'N2O_AO1'};
				tits   = {'$N_2O_{nit}$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'N2OAMMOX','N2OAO1_CONS'};
				smseq  = [(1) (-1)];
				fluxes = {'FG_N2O_AO1'};
				lvls   = {   'sfc'};		
			elseif strcmp(varname,'N2O_ATM')
				vars   = {'N2O_ATM'};
				tits   = {'$N_2O_{atm}$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'N2OATM_CONS'};
				smseq  = [(-1)];
				fluxes = {'FG_N2O_ATM'};
				lvls   = {   'sfc'};
			elseif strcmp(varname,'N2O_SIDEN')
				vars   = {'N2O_SIDEN'};
				tits   = {'$N_2O_{bou}$'};
				units  = '$mmol$ $N$ $m^{-2}$ $s^{-1}$'; 
				rates  = {'N2OSIDEN_CONS'};
				smseq  = [(-1)];
				fluxes = {'FG_N2O_SIDEN'};
				lvls   = {   'sfc'};
			end

			% Compute dzdt
			obj = computeDzDt(obj);

			% Load all 3D output terms
			terms = [vars,rates,fluxes];
			terms = [terms(find(~cellfun(@isempty,terms)))];
			obj = loadData(obj,terms,'type','raw');

			% Integrate 3D variables, rates vertically
			terms = [vars,rates];
			terms = [terms(find(~cellfun(@isempty,terms)))];
			obj = intVar(obj,terms);

			% Load and process 2D fluxes
			obj = getFluxes(obj,varname,fluxes,lvls);

			% Get dCdt, advection, sms, net terms
			obj = computeDcDt(obj,vars);
			obj = computeXYZflux(obj,vars);
			obj = computeSMS(obj,varname,rates,smseq);
			obj = computeNet(obj,varname);

			% Integrate vertically (and horizontally)
			obj = intBudg(obj,varname); 
		end % end method getBudg

		%--------------------------------------------------------------------------------
		function obj = getFluxes(obj,varname,fluxes,lvls);
			% -------------------
			% Grab 2D fluxes for budget, convert to 3D 
			% Called in getBudg
			% Fluxes in mmol/m2/s
			%
			% Usage:
			% - obj = getFluxes(obj,varname,fluxes,lvls)
			%
			% Inputs:
			% - varname = BGC budget that you are closing (i.e. 'NO2')
			% - fluxes = 2D fluxes (air-sea, sediment, etc)
			% - lvls   = levels to apply 2D flux ('sfc','sed', as a cell array)
			%
			% Example:
			% - obj = getFluxes(obj,'N2O',{'FG_N2O'},{'sfc'})
			% -------------------
			
			disp('---------------------------------');
			disp('Get 2D fluxes, apply to 3D grid')

			% Convert 2D flux to 3D based on lvls
			for i = 1:length(fluxes)
			
				% Initialize matrices-to-fill
				obj.romsData.(varname).budget.fg = zeros(size(obj.region.mask_rho3d)).*obj.region.mask_rho3d;
				obj.romsData.(varname).budget.sed = zeros(size(obj.region.mask_rho3d)).*obj.region.mask_rho3d;
	
				% Apply 2D mask to 2D data
				obj.romsData.(fluxes{i}).data = obj.romsData.(fluxes{i}).data .* obj.region.mask_rho;
				
				% Apply 2D flux to correct z-level to make 3D
				if strcmp(lvls{i},'sfc')
					tmpfg = zeros(size(obj.region.mask_rho3d));
					% Apply value into 3D grid
					tmpfg(:,:,obj.region.nz,:) = obj.romsData.(fluxes{i}).data .* obj.region.mask_rho;
					% Divide by z, save as 3D rate
					obj.romsData.(varname).budget.fg = tmpfg ./ obj.region.dz;
				elseif strcmp(lvls{i},'sed')
					tmpsed = zeros(size(obj.region.mask_rho3d));
					% Apply value into 3D grid
					tmpsed(:,:,1,:) = obj.romsData.(fluxes{i}).data .* obj.region.mask_rho;
					% Divide by z, save as 3D rate
					obj.romsData.(varname).budget.sed = tmpsed ./ obj.region.dz;
				end
			end
		end % end method getFluxes

		%--------------------------------------------------------------------------------
		function obj = computeDzDt(obj)
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
			
			% Grab time info from average file
			stats      = ncinfo(obj.paths.avg);
			fieldnames = {stats.Attributes.Name};
			fieldvalue = {stats.Attributes.Value};
			dt         = fieldvalue(strcmp(fieldnames,'dt'));
			dt         = double(dt{1});
			navg       = fieldvalue(strcmp(fieldnames,'navg'));
			navg       = double(navg{1});

			% Compute dt according to output frequency
			dt         = dt*navg;

			% Calculate dzdt
			avg_z = obj.region.z_w;
			
			% Save levels
			obj.region.dz = diff(avg_z,1,3);
			obj.region.dz = obj.region.dz .* obj.region.mask_rho3d; % apply 3D mask
			obj.region.dt = dt;
		end % end methods computeDzDt

		%--------------------------------------------------------------------------------
		function obj = computeDcDt(obj,vars)
			% -------------------
			% Compute change in concentration with time
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeDcDt(obj,vars)
			%
			% Inputs:
			% - vars = variable to calculated dC/dt for
			%
			% Example:
			% - obj = computeDcDt(obj,'NO2');
			% -------------------
	
			disp('---------------------------------');
			disp('Computing dC/dt')
			
			% Load data on either side of 'history' (first,last snapshot)
			% Scroll through each variable and get dCdt
			for v = 1:length(vars)
				
				% Get start/finish from history
				his1 = double(ncread(obj.paths.his,vars{v},...
						[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1],...
					   [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf obj.region.nt]));
				
				his2 = double(ncread(obj.paths.his,vars{v},...
					   [obj.region.lon_lim(1) obj.region.lat_lim(1) 1 2],...
					   [diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf obj.region.nt]));

				% Divide by dt, apply 3d mask
				obj.romsData.(vars{v}).budget.dcdt = ((his2 - his1)/obj.region.dt) .* obj.region.mask_rho3d;
			end 
		end % end method computeDcDt

		%--------------------------------------------------------------------------------
		function obj = computeXYZflux(obj,vars)
			% -------------------
			% Compute HorXAdvFlux, HorYAdvFlux, and top/bottom advection
			% Also get diffusion, if it is available
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeXYZflux(obj,vars)
			%
			% Inputs:
			% - vars = variable to get flux terms for 
			%
			% Example:
			% - obj = computeXYZflux(obj,'NO2');
			% -------------------
			
			disp('---------------------------------');
			disp('Get advective/diffusion terms')
			
			% Cycle through and load XYfluxes
			for f = 1:length(vars)
				% Load XYZ advection
				temp.X  = double(ncread(obj.paths.flux,['HorXAdvFlux_',vars{f}],...
								[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1  ],...
								[diff(obj.region.lon_lim)+2 diff(obj.region.lat_lim)+1 inf inf]));
				temp.Y  = double(ncread(obj.paths.flux,['HorYAdvFlux_',vars{f}],...
								[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1  ],...
								[diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+2 inf inf]));
				temp.Z  = double(ncread(obj.paths.flux,['VertAdvFlux_',vars{f}],...
								[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1  ],...
								[diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf inf]));
				temp.Zd = double(ncread(obj.paths.flux,['VertDiffFlux_',vars{f}],...
								[obj.region.lon_lim(1) obj.region.lat_lim(1) 1 1  ],...
								[diff(obj.region.lon_lim)+1 diff(obj.region.lat_lim)+1 inf inf]));

				% Initiate matrices-to-fill, compute adv X and Y
				% Converts from mmol/s to mmol/m3/s by dividing by grid area and height of cell (volume)
				% Get dimensions
				nx = obj.region.nx;
				ny = obj.region.ny;
				nz = obj.region.nz;
				nt = obj.region.nt;
				
				% X advection
				adx(1:nx,1:ny,1:nz,1:nt) = NaN;
				adx(1:nx,:,:,1:nt)       = (temp.X(1:nx,:,:,:) - temp.X(2:nx+1,:,:,:)) ./ ...
											obj.region.area3d(1:nx,:,:,:) ./ obj.region.dz(1:nx,:,:,:);
				
				% Y advection
				ady(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				ady(:,1:ny,:,:)          = (temp.Y(:,1:ny,:,:) - temp.Y(:,2:ny+1,:,:)) ./ ...
											obj.region.area3d(:,1:ny,:,:) ./ obj.region.dz(:,1:ny,:,:);
				
				% Z advection
				adz(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				adz(:,:,:,:)             = (temp.Z(:,:,1:nz,:) - temp.Z(:,:,2:nz+1,:)) ./ obj.region.dz(:,:,:,:);
				
				% Z diffusion
				dfz(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				dfz(:,:,:,:)             = (temp.Zd(:,:,1:nz,:) - temp.Zd(:,:,2:nz+1,:)) ./ obj.region.dz(:,:,:,:);
				
				% Apply 3D mask
				obj.romsData.(vars{f}).budget.adx = adx .* obj.region.mask_rho3d(:,:,:,:);
				obj.romsData.(vars{f}).budget.ady = ady .* obj.region.mask_rho3d(:,:,:,:);
				obj.romsData.(vars{f}).budget.adz = adz .* obj.region.mask_rho3d(:,:,:,:);
				obj.romsData.(vars{f}).budget.dfz = dfz .* obj.region.mask_rho3d(:,:,:,:);
			end
		end % end methods computeXYZFlux

		%--------------------------------------------------------------------------------
		function obj = computeSMS(obj,varname,rates,smseq)
			% ---------------------
			% Gathers sources and sinks
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeSMS(obj,rates,smseq)
			%
			% Inputs:
			% - varname = budget variable (i.e. 'NO2')
			% - rates = BGC rates, set in getBudg
			% - smseq = S-M-S equation, set in getBudg
			%
			% Example:
			% - varname = 'N2O_AO1';
			% - rates = {'N2OAMMOX','N2OAO1_CONS');
			% - smseq = [(1) (-1)];
			% - obj = computeSMS(obj,varname,rates,smseq);
			% ---------------------

			disp('---------------------------------');
			disp('Computing sources-minus-sinks (SMS)');

			% Write equation
			eq = [];
			smseq_orig = smseq;
			for i = 1:length(smseq)
				if smseq(i) > 0
					sym = ' + ';
				elseif smseq(i) < 0 
					sym = ' - ';
					smseq(i) = -smseq(i);
				end
				if i == 1
					eq = ['(',num2str(smseq(i)),')*',rates{i}];
				else
					eq = [eq,sym,'(',num2str(smseq(i)),')*',rates{i}];
				end
			end
			smseq = smseq_orig;
			obj.romsData.(varname).budget.info.sms_eq = eq;
			obj.romsData.(varname).budget.info.sms_factors = smseq;

			% Get SMS 
			dims = ndims(obj.romsData.(rates{1}).data);
			eq = [];
			for i = 1:length(rates)
				eq{i} = [(smseq(i).*obj.romsData.(rates{i}).data)];
			end
			tmpsms = sum(cat(dims+1,eq{:}),dims+1);
			obj.romsData.(varname).budget.sms = tmpsms;

			% Get production 
			ind = find(smseq>0);
			if ~isempty(ind);
				tmprates = rates(ind);
				obj.romsData.(varname).budget.info.prod_sources = tmprates;
				obj.romsData.(varname).budget.info.prod_factors = smseq(ind);
				eq = [];
				for i = 1:length(tmprates)
					eq{i} = [(smseq(ind(i)).*obj.romsData.(tmprates{i}).data)];
				end
				tmpprod = sum(cat(dims+1,eq{:}),dims+1);
				obj.romsData.(varname).budget.prod = tmpprod;
			else
				obj.romsData.(varname).budget.info.prod_sources = [];
				obj.romsData.(varname).budget.info.prod_factors = []; 
				obj.romsData.(varname).budget.prod = [];
			end

			% Get consumption 
			ind = find(smseq<0);
			if ~isempty(ind)
				tmprates = rates(ind);
				obj.romsData.(varname).budget.info.cons_sources = tmprates;
				obj.romsData.(varname).budget.info.cons_factors = smseq(ind);
				eq = [];
				for i = 1:length(tmprates)
					eq{i} = [(smseq(ind(i)).*obj.romsData.(tmprates{i}).data)];
				end
				tmpcons = sum(cat(dims+1,eq{:}),dims+1);
				obj.romsData.(varname).budget.cons = tmpcons;
			else
				obj.romsData.(varname).budget.info.cons_sources = [];
				obj.romsData.(varname).budget.info.cons_factors = []; 
				obj.romsData.(varname).budget.cons = [];
			end
		end % end method getSMS

		%--------------------------------------------------------------------------------
		function obj = computeNet(obj,varname)
			% ---------------------
			% Computes remainder (net) from budget equation
			% dC/dt = adv + dfz + sms + fg
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeNet(obj,varname)
			%
			% Inputs:
			% - varname = budget to close (i.e 'NO2');
			%
			% Example:
			% - obj = computeNet(obj,'NO2');
			% ---------------------

			disp('---------------------------------');
			disp('Computing budget remainder (net)');

			% Calculate remainder (net)
			obj.romsData.(varname).budget.net =   obj.romsData.(varname).budget.dcdt - (obj.romsData.(varname).budget.adx + ...
												  obj.romsData.(varname).budget.ady  +  obj.romsData.(varname).budget.adz + ...
												  obj.romsData.(varname).budget.dfz  +  obj.romsData.(varname).budget.sms + ...
												  obj.romsData.(varname).budget.fg   +  obj.romsData.(varname).budget.sed);
		end % end method computeNet
	
		%--------------------------------------------------------------------------------
		function obj = intVar(obj,vars)
			% ------------------
			% Vertically integrate 3D variable(s) 
			%
			% Usage:
			% - obj = intVar(obj,vars)
			%
			% Inputs:
			% - vars = 3D variables to integrate, as a cell array
			%
			% Example:
			% - obj = loadData(obj,{'NO3'},'type','raw');
			% - obj = computeDzDt(obj);
			% - obj = intVar(obj,{'NO3'});
			% ------------------

			disp('---------------------------------');
			disp('Integrating 3D variables');
			
			% Go through each 3D rate, integrate vertically (and totally)
			for i = 1:length(vars)		
				tmpdata           	       = obj.romsData.(vars{i}).data .* obj.region.dz .* obj.region.mask_rho3d; % mmol/m3/s --> mmol/m2/s
				tmpdata           	       = squeeze(nansum(tmpdata,3));
				tmpdata                    = tmpdata .* obj.region.mask_rho;
				obj.romsData.(vars{i}).int = tmpdata; 
				for t = 1:obj.region.nt
					obj.romsData.(vars{i}).tot(t) = nansum(obj.romsData.(vars{i}).int(:,:,t) .*obj.region.grid_area,'all');
				end
			end
		end % end method intVar

		%--------------------------------------------------------------------------------
		function obj = intBudg(obj,varname,terms)
			% ------------------
			% Vertically integrate budget terms
			% Called in getBudg
			% Terms are all in mmol/m3/s, get it in mmol/m2/s by using dz
			%
			% Usage:
			% - obj = intBudg(obj)
			%
			% Inputs:
			% - varname = budget variable (i.e. 'NO2');
			% - terms = (optional) cell array of budget terms to integrate
			%		   ...use if you want to exclude other terms
			%		   ...to integrate all terms, exclude from command
			%
			% Example:
			% - obj = intBudg(obj,'NO2',{'sms'}); <-- only sms 
			% - obj = intBudg(obj,'NO2');         <-- all terms
			% ------------------

			if nargin >1 & nargin <3
				terms = {'dcdt','adx','ady','adz','dfz','sms','prod','cons','fg','sed','net'};;
			elseif nargin <1
				disp('Must supply varname');
				return
			end

			disp('---------------------------------');
			disp('Computing vertical integration...');
			disp('...and getting grid area fluxes');
	
			% Grab terms, including remainder
			for i = 1:length(terms)
				eval([terms{i},' = obj.romsData.(varname).budget.',terms{i},' .* obj.region.mask_rho3d;']);
			end

			% Integrate vertically (fg term should match real flux)
			% ...mmol/m3/s to mmol/m2/s
			for i = 1:length(terms)
				eval(['obj.romsData.(varname).budget.int',terms{i},' = squeeze(nansum(',terms{i},'.*obj.region.dz,3));']);
			end
			
			% ...mmol/m2/s to mmol/s
			for i = 1:length(terms);
				for t = 1:obj.region.nt
					eval(['obj.romsData.(varname).budget.tot',terms{i},'(t)  = nansum(obj.romsData.(varname).budget.int',terms{i},...
						  '(:,:,t) .*obj.region.grid_area,''all'');']); 
				end
			end

			% Check for total advection
			if sum(ismember({'adx','ady','adz'},terms)) == 3
				obj.romsData.(varname).budget.intadv  =  obj.romsData.(varname).budget.intadx + ...
														 obj.romsData.(varname).budget.intady + ...
														 obj.romsData.(varname).budget.intadz;
				obj.romsData.(varname).budget.totadv  =	 obj.romsData.(varname).budget.totadx + ...
														 obj.romsData.(varname).budget.totady + ...
														 obj.romsData.(varname).budget.totadz;
			end
		end

		%--------------------------------------------------------------------------------
		function plotBudg(obj,varargin);
			% ------------------
			% Plot the intergrated budget terms
			%
			% Usage:
			% - plotBudg(obj,varargin)
			%
			% Inputs (varargin):
			% - time = time record to plot (default--> all)
			% - prc  = percentile to limit colorbar
			%	   ...if prc == 5, then caxis = [5th %, 95th %]
			%
			% Example:
			% - plotBudg(obj,'time',10,'prc',5)
			% ------------------

			% defaults for optional  arguments
			A.time      = [];
			A.prc       = [2];
			A.lonbounds = [floor(min(obj.region.lon_rho(:))) ceil(max(obj.region.lon_rho(:)))];
			A.latbounds = [floor(min(obj.region.lat_rho(:))) ceil(max(obj.region.lat_rho(:)))];
			A.tfont     = 18;
			A.cbfont    = 14;
			A           = parse_pv_pairs(A,varargin); % parse method arguments to A
			
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
					dcdt = obj.budget.intdcdt.(vars{i})(:,:,tt);
					adv  = obj.budget.intadx.(vars{i})(:,:,tt) + ...
						   obj.budget.intady.(vars{i})(:,:,tt) + ...
						   obj.budget.intadz.(vars{i})(:,:,tt);
					dfz  = obj.budget.intdfz.(vars{i})(:,:,tt);
					sms  = obj.budget.intsms.(vars{i})(:,:,tt);
					fg   = obj.budget.intfg.(vars{i})(:,:,tt);
					sed  = obj.budget.intsed.(vars{i})(:,:,tt);
					net  = obj.budget.intnet.(vars{i})(:,:,tt);
					
					% Blank land
					dcdt = dcdt .* obj.region.mask_rho;
					adv  = adv  .* obj.region.mask_rho;
					dfz  = dfz  .* obj.region.mask_rho;
					sms  = sms  .* obj.region.mask_rho;
					fg   = fg   .* obj.region.mask_rho;
					sed  = sed  .* obj.region.mask_rho;
					net  = net  .* obj.region.mask_rho;

					% Get colobars
					cbar_dcdt = romsMaster.prclims(dcdt,'prc',A.prc);
					cbar_adv  = romsMaster.prclims(adv,'prc',A.prc);
					cbar_dfz  = romsMaster.prclims(dfz,'prc',A.prc);
					cbar_sms  = romsMaster.prclims(sms,'prc',A.prc);
					cbar_fg   = romsMaster.prclims(fg,'prc',A.prc);
					cbar_sed  = romsMaster.prclims(sed,'prc',A.prc);
					cbar_net  = romsMaster.prclims(net,'prc',A.prc);
					if t == 1
						cbar.dcdt = cbar_dcdt;
						cbar.adv  = cbar_adv;
						cbar.dfz  = cbar_dfz;
						cbar.sms  = cbar_sms;
						cbar.fg   = cbar_fg;
						cbar.sed  = cbar_sed;
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
					if max(cbar_sed) > max(cbar.sed);
						cbar.sed  = cbar_sed;
					end
					if max(cbar_net) > max(cbar.net)
						cbar.net = cbar_net;
					end
				end

				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	

					% Gather terms
					dcdt = obj.budget.intdcdt.(vars{i})(:,:,tt);
					adv  = obj.budget.intadx.(vars{i})(:,:,tt) + ...
						   obj.budget.intady.(vars{i})(:,:,tt) + ...
						   obj.budget.intadz.(vars{i})(:,:,tt);
					dfz  = obj.budget.intdfz.(vars{i})(:,:,tt);
					sms  = obj.budget.intsms.(vars{i})(:,:,tt);
					fg   = obj.budget.intfg.(vars{i})(:,:,tt);
					sed  = obj.budget.intsed.(vars{i})(:,:,tt);
					net  = obj.budget.intnet.(vars{i})(:,:,tt);

					% Blank land
					dcdt = dcdt .* obj.region.mask_rho;
					adv  = adv  .* obj.region.mask_rho;
					dfz  = dfz  .* obj.region.mask_rho;
					sms  = sms  .* obj.region.mask_rho;
					fg   = fg   .* obj.region.mask_rho;
					sed  = sed  .* obj.region.mask_rho;
					net  = net  .* obj.region.mask_rho;
					
					% Start plots
					for j = 1:7
				
						% Initiate map figure
						fig = piofigs('lfig',1);
						if j == 1
							dat    = dcdt;
							titstr = ['Integrated ',tits{i},': change over time'];
							fstr   = ['dcdt'];
							clevs  = cbar.dcdt;
						elseif j == 2
							dat    = adv;
							titstr = ['Integrated ',tits{i},': advection'];
							fstr   = ['adv'];
							clevs  = cbar.adv;
						elseif j == 3
							dat    = dfz;
							titstr = ['Integrated ',tits{i},': diffusion'];
							fstr   = ['diff'];
							clevs  = cbar.dfz;
						elseif j == 4
							dat    = sms;
							titstr = ['Integrated ',tits{i},': sources-minus-sinks'];
							fstr   = ['sms'];
							clevs  = cbar.sms;
						elseif j == 5
							dat    = fg;
							titstr = ['Integrated ',tits{i},': air-sea flux'];
							fstr   = ['fg'];
							clevs  = cbar.fg;
						elseif j == 6
							dat    = sed;
							titstr = ['Integrated ',tits{i},': sediment flux'];
							fstr   = ['sed'];
							clevs  = cbar.sed;
						elseif j == 7
							dat    = net;
							titstr = ['Integrated ',tits{i},': net remainder'];
							fstr   = ['net'];
							clevs  = cbar.net;
						end
				
						% Skip empty datasets
						if isempty(clevs) | diff(clevs)==0	
							continue;
						end
						
						% Plot results
						[ax] = map_plot(fig(1),obj.region.lon_rho,obj.region.lat_rho,...
									'lonbounds',A.lonbounds,'latbounds',A.latbounds);
						clevs        	    = clevs(1):(diff(clevs)/31):clevs(2);
						dat(dat<clevs(1))   = clevs(1);
						dat(dat>clevs(end)) = clevs(end);
						[tmp hc]            = m_contourf(obj.region.lon_rho,obj.region.lat_rho,dat,clevs,'LineStyle','none');
						% Plot region box
						m_plot(obj.region.lon_rho(1,:),obj.region.lat_rho(1,:),'--k','linewidth',2);
						m_plot(obj.region.lon_rho(:,1),obj.region.lat_rho(:,1),'--k','linewidth',2);
						m_plot(obj.region.lon_rho(end,:),obj.region.lat_rho(end,:),'--k','linewidth',2);
						m_plot(obj.region.lon_rho(:,end),obj.region.lat_rho(:,end),'--k','linewidth',2);
						cb = colorbar;
						title({titstr,[obj.info.time_string_idv{tt},' ',num2str(obj.info.Year)]},...
							  'Interpreter', 'Latex','FontSize',A.tfont);
						ylabel(cb,units,'Interpreter','Latex');
						caxis([clevs(1) clevs(end)]);
						colormap(gca,cmocean('balance'));
						ax.FontSize = A.cbfont;
						cb.FontSize = A.cbfont;
						% Print	
						if strcmp(obj.info.time_string,'Monthly')
							fname = [vars{i},'_',fstr,'_M',num2str(tt)];
						elseif strcmp(obj.info.time_string,'Daily');
							fname = [vars{i},'_',fstr,'_D',num2str(tt)];
						end
						export_fig('-jpg',[obj.paths.plots.budget,fname]);
						vecrast(gcf,[obj.paths.plots.budget,fname], 300, 'bottom', 'pdf');
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
			% Inputs (varargin):
			% - vars - variable(s) to plot
			% - time - time record to plot (default--> all)
			% - lon  - longitude point (can be a vector)
			% - lat  - latitude point (can be a vector)
			% - zlim - limit in depth space (i.e. [1200] for 0 - 1200)
			%
			% - Example:
			%   plot1D(obj,'vars',{'NO2','N2O'},'time',1:12,'lon',200,'lat',0)
			% ------------------

			% Defaults for optional  arguments
			A.vars  = [];
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
			
			% Load data
			obj = loadData(obj,A.vars,'type','raw');	

			% Get index of nearest lat/lon grid cell
			lon_idx = [];
			lat_idx = [];
			for i = 1:length(A.lon)
				all_dist   = distance(A.lat(i),A.lon(i),obj.region.lat_rho,obj.region.lon_rho);
				[lon_idx(i),lat_idx(i)]  = find(all_dist == min(all_dist(:)));
			end

			% Get colormix
			clrs = colormap(cmocean('phase'));
			idx  = floor(linspace(1,length(clrs),length(A.time)+1));
			clrs = clrs(idx,:); clrs(end,:) = [];

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

			% Plot results
			for i = 1:length(A.lon)
				for j = 1:length(A.vars)
					for k = 1:length(A.time)
				
						% Gather data
						tmpdata = squeeze(obj.romsData.(A.vars{j}).data(lon_idx(i),lat_idx(i),:,A.time(k)));

						% Gather depth
						depth = -squeeze(obj.region.z_r(lon_idx(i),lat_idx(i),:,A.time(k)));

						% Make plot(s)
						if k == 1
							fig = piofigs('lfig',1.5);
							set(gca,'YDir','Reverse');
							xlabel(obj.romsData.(A.vars{j}).units);
							ylabel('Depth (m)');
							hold on; grid on
						end
						if length(A.time) == 1
							title({[A.vars{j},' vs. Depth: ',tstr{A.time(k)},' Average'],...
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
							fname = [A.vars{j},'_vs_z_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_',tstr{A.time(k)}];
							print('-djpeg',[obj.paths.plots.roms.profile,fname])
							close all
						end
					end
					
					% Add averages if many months plotted
					if length(A.time) == obj.region.nt
						% Gather data
						tmpdata = squeeze(obj.romsData.(A.vars{j}).data(lon_idx(i),lat_idx(i),:,:));
						tmpdata = nanmean(tmpdata,2);
						plot(tmpdata,depth,'--k','linewidth',4);
						lstr = [tstr];
						lstr(13) = {'Avg'};
						legend(lstr,'Location','Southeast');
						title({[A.vars{j},' vs. Depth'],...
							   ['Lon=',num2str(A.lon(i)),', Lat=',num2str(A.lat(i))]});
						fname = [A.vars{j},'_vs_z_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_annual'];
						print('-djpeg',[obj.paths.plots.roms.profile,fname])
						close all
					end
				end

				% Plot locations
				fig = piofigs('lfig',1);
				plot(A.lon,A.lat,'sk','markersize',10);
				hold on
				xlim([obj.grid.minlon_rho obj.grid.maxlon_rho]);
				ylim([obj.grid.minlat_rho obj.grid.maxlat_rho]);
				plot_coast
				pltjpg(i);
				disp('Map saved in: /data/project1/demccoy/tmpfigs');
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
			
			% check that object is initialized
			try; x = obj.info.var3d; catch; obj = init(obj); end

			% temperature
			obj.paths.diag.temp.file  = {'/data/project3/data/woa18/temperature/0p25/temp_woa18_clim.nc'};
			obj.paths.diag.temp.type  = {'nc'};
			obj.paths.diag.temp.var   = {'temp'};
			obj.paths.diag.temp.zvar  = {'depth'};
			obj.paths.diag.temp.dim   = {'xyzt'};
			obj.paths.diag.temp.lon   = {'lon'};
			obj.paths.diag.temp.lat   = {'lat'};
			obj.paths.diag.temp.name  = {'WOA-18 Temperature'};
			obj.paths.diag.temp.units = {'$^oC$'};

			% salinity
			obj.paths.diag.salt.file  = {'/data/project3/data/woa18/salinity/0p25/salt_woa18_clim.nc'};
			obj.paths.diag.salt.type  = {'nc'};
			obj.paths.diag.salt.var   = {'salt'};
			obj.paths.diag.salt.zvar  = {'depth'};
			obj.paths.diag.salt.dim   = {'xyzt'};
			obj.paths.diag.salt.lon   = {'lon'};
			obj.paths.diag.salt.lat   = {'lat'};
			obj.paths.diag.salt.name  = {'WOA-18 Salinity'};
			obj.paths.diag.salt.units = {'$PSU$'};

			% density
			obj.paths.diag.sigma.file  = {'/data/project3/data/woa18/density/0p25/sigma_woa18_clim.nc'};
			obj.paths.diag.sigma.type  = {'nc'};
			obj.paths.diag.sigma.var   = {'sigma'};
			obj.paths.diag.sigma.zvar  = {'depth'};
			obj.paths.diag.sigma.dim   = {'xyzt'};
			obj.paths.diag.sigma.lon   = {'lon'};
			obj.paths.diag.sigma.lat   = {'lat'};
			obj.paths.diag.sigma.name  = {'WOA-18 Density'};
			obj.paths.diag.sigma.units = {'$kg$ $m^{-3}$'};

			% sea surface height
			obj.paths.diag.SSH.file  = {'/data/project1/demccoy/ROMS/validation/AVISO/monthly_AVISO.mat'};
			obj.paths.diag.SSH.type  = {'mat'};
			obj.paths.diag.SSH.var   = {'adt_month_av'};
			obj.paths.diag.SSH.zvar  = {[]};
			obj.paths.diag.SSH.dim   = {'xyt'};
			obj.paths.diag.SSH.lon   = {'lon_aviso'};
			obj.paths.diag.SSH.lat   = {'lat_aviso'};
			obj.paths.diag.SSH.name  = {'AVISO Absolute Dynamic Topography'};
			obj.paths.diag.SSH.units = {'$m$'};

			% zonal wind stress
			obj.paths.diag.zws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
			obj.paths.diag.zws.type  = {'mat'};
			obj.paths.diag.zws.var   = {'u'};
			obj.paths.diag.zws.zvar  = {[]};
			obj.paths.diag.zws.dim   = {'xyt'};
			obj.paths.diag.zws.lon   = {'lon'};
			obj.paths.diag.zws.lat   = {'lat'};
			obj.paths.diag.zws.name  = {'SCOW-2010 Zonal Wind Stress'};
			obj.paths.diag.zws.units = {'$N$ $m^{-2}$'};

			% meridional wind stress
			obj.paths.diag.mws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
			obj.paths.diag.mws.type  = {'mat'};
			obj.paths.diag.mws.var   = {'v'};
			obj.paths.diag.mws.zvar  = {[]};
			obj.paths.diag.mws.dim   = {'xyt'};
			obj.paths.diag.mws.lon   = {'lon'};
			obj.paths.diag.mws.lat   = {'lat'};
			obj.paths.diag.mws.name  = {'SCOW-2010 Meridional Wind Stress'};
			obj.paths.diag.mws.units = {'$N$ $m^{-2}$'};
			
			% wind stress
			obj.paths.diag.ws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
			obj.paths.diag.ws.type  = {'mat'};
			obj.paths.diag.ws.var   = {'ws'};
			obj.paths.diag.ws.zvar  = {[]};
			obj.paths.diag.ws.dim   = {'xyt'};
			obj.paths.diag.ws.lon   = {'lon'};
			obj.paths.diag.ws.lat   = {'lat'};
			obj.paths.diag.ws.name  = {'SCOW-2010 Wind Stress'};
			obj.paths.diag.ws.units = {'$N$ $m^{-2}$'};
		
			% wind stress curl
			obj.paths.diag.wsc.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
			obj.paths.diag.wsc.type  = {'mat'};
			obj.paths.diag.wsc.var   = {'wsc'};
			obj.paths.diag.wsc.zvar  = {[]};
			obj.paths.diag.wsc.dim   = {'xyt'};
			obj.paths.diag.wsc.lon   = {'lon'};
			obj.paths.diag.wsc.lat   = {'lat'};
			obj.paths.diag.wsc.name  = {'SCOW-2010 Wind Stress Curl'};
			obj.paths.diag.wsc.units = {'$N$ $m^{-2}$'};

			% u velocity
			obj.paths.diag.u.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
			obj.paths.diag.u.type  = {'mat'};
			obj.paths.diag.u.var   = {'U'};
			obj.paths.diag.u.zvar  = {'dep'};
			obj.paths.diag.u.dim   = {'xyzt'};
			obj.paths.diag.u.lon   = {'lon'};
			obj.paths.diag.u.lat   = {'lat'};
			obj.paths.diag.u.name  = {'GODAS U-Velocity'};
			obj.paths.diag.u.units = {'$m$ $s^{-1}'}; 
	
			% v velocity
			obj.paths.diag.v.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
			obj.paths.diag.v.type  = {'mat'};
			obj.paths.diag.v.var   = {'V'};
			obj.paths.diag.v.zvar  = {'dep'};
			obj.paths.diag.v.dim   = {'xyzt'};
			obj.paths.diag.v.lon   = {'lon'};
			obj.paths.diag.v.lat   = {'lat'};
			obj.paths.diag.v.name  = {'GODAS V-Velocity'};
			obj.paths.diag.v.units = {'$m$ $s^{-1}'}; 

			% oxygen (only 1.00 available)
			obj.paths.diag.O2.file  = {'/data/project3/data/woa18/oxygen/1p0/o2_woa18_clim.nc'};
			obj.paths.diag.O2.type  = {'nc'};
			obj.paths.diag.O2.var   = {'o2'};
			obj.paths.diag.O2.zvar  = {'depth'};
			obj.paths.diag.O2.dim   = {'xyzt'};
			obj.paths.diag.O2.lon   = {'lon'};
			obj.paths.diag.O2.lat   = {'lat'};
			obj.paths.diag.O2.name  = {'WOA-18 Oxygen'};
			obj.paths.diag.O2.units = {'$mmol$ $m^{-3}$'};

			% nitrate (only 1.00 available)
			obj.paths.diag.NO3.file  = {'/data/project3/data/woa18/nitrate/1p0/no3_woa18_clim.nc'};
			obj.paths.diag.NO3.type  = {'nc'};
			obj.paths.diag.NO3.var   = {'no3'};
			obj.paths.diag.NO3.zvar  = {'depth'};
			obj.paths.diag.NO3.dim   = {'xyzt'};
			obj.paths.diag.NO3.lon   = {'lon'};
			obj.paths.diag.NO3.lat   = {'lat'};
			obj.paths.diag.NO3.name  = {'WOA-18 Nitrate'};
			obj.paths.diag.NO3.units = {'$mmol$ $m^{-3}$'};

			% phosphate (only 1.00 available)
			obj.paths.diag.PO4.file  = {'/data/project3/data/woa18/phosphate/1p0/po4_woa18_clim.nc'};
			obj.paths.diag.PO4.type  = {'nc'};
			obj.paths.diag.PO4.var   = {'po4'};
			obj.paths.diag.PO4.zvar  = {'depth'};
			obj.paths.diag.PO4.dim   = {'xyzt'};
			obj.paths.diag.PO4.lon   = {'lon'};
			obj.paths.diag.PO4.lat   = {'lat'};
			obj.paths.diag.PO4.name  = {'WOA-18 Phosphate'};
			obj.paths.diag.PO4.units = {'$mmol$ $m^{-3}$'};
 
			% N* (only 1.00 available)
			obj.paths.diag.nstar.file  = {'/data/project3/data/woa18/nstar/1p0/nstar_woa18_clim.nc'};
			obj.paths.diag.nstar.type  = {'nc'};
			obj.paths.diag.nstar.var   = {'nstar'};
			obj.paths.diag.nstar.zvar  = {'depth'};
			obj.paths.diag.nstar.dim   = {'xyzt'};
			obj.paths.diag.nstar.lon   = {'lon'};
			obj.paths.diag.nstar.lat   = {'lat'};
			obj.paths.diag.nstar.name  = {'WOA-18 N*'};
			obj.paths.diag.nstar.units = {'$mmol$ $m^{-3}$'};
						 
			% silicate (only 1.00 available)
			obj.paths.diag.SiO3.file  = {'/data/project3/data/woa18/silicate/1p0/si_woa18_clim.nc'};
			obj.paths.diag.SiO3.type  = {'nc'};
			obj.paths.diag.SiO3.var   = {'si'};
			obj.paths.diag.SiO3.zvar  = {'depth'};
			obj.paths.diag.SiO3.dim   = {'xyzt'};
			obj.paths.diag.SiO3.lon   = {'lon'};
			obj.paths.diag.SiO3.lat   = {'lat'};
			obj.paths.diag.SiO3.name  = {'WOA-18 Silicate'};
			obj.paths.diag.SiO3.units = {'$mmol$ $m^{-3}$'};

			% N2O (only 1.00 available, reconstructed)
			%obj.paths.diag.N2O.file   = {'/data/project1/dclements/N_Interior/n2opredRF_02-21-2021.nc'};
			%obj.paths.diag.N2O.file   = {'/data/project1/demccoy/ROMS/validation/n2o/n2o_NN_format.mat'};
			%obj.paths.diag.N2O.type   = {'mat'};
			obj.paths.diag.N2O.file   = {'/data/project2/yangsi/analysis/n2oInterior/processed/n2opredRF_05-27-2020.nc'};
			obj.paths.diag.N2O.type   = {'nc'};
			obj.paths.diag.N2O.var    = {'n2o'};
			obj.paths.diag.N2O.zvar   = {'depth'};
			obj.paths.diag.N2O.dim    = {'xyzt'};
			obj.paths.diag.N2O.lon    = {'lon'};
			obj.paths.diag.N2O.lat    = {'lat'};
			obj.paths.diag.N2O.name   = {'Clements et al. Nitrous Oxide'};
			obj.paths.diag.N2O.units  = {'$mmol$ $m^{-3}$'};
			obj.paths.diag.N2O.factor = {0.001}; %nmol to mmol

			% NO2 (only 1.00 available, reconstructed)
			%obj.paths.diag.NO2.file  = {'/data/project1/demccoy/ROMS/validation/no2/no2_NN_format.mat'};
			%obj.paths.diag.NO2.type  = {'mat'};
			obj.paths.diag.NO2.file  = {'/data/project2/yangsi/analysis/no2Interior/processed/no2predRF_05-27-2020.nc'};
			obj.paths.diag.NO2.type  = {'nc'};
			obj.paths.diag.NO2.var   = {'no2'};
			obj.paths.diag.NO2.zvar  = {'depth'};
			obj.paths.diag.NO2.dim   = {'xyzt'};
			obj.paths.diag.NO2.lon   = {'lon'};
			obj.paths.diag.NO2.lat   = {'lat'};
			obj.paths.diag.NO2.name  = {'Clements et al. Nitrite'};
			obj.paths.diag.NO2.units = {'$mmol$ $NO_2 m^{-3}$'};

			% Chl 
			obj.paths.diag.SFC_CHL.file  = {'/data/project1/data/MODIS-Aqua/backup/res0p25/chl_clim_0p25.nc'};
			obj.paths.diag.SFC_CHL.type  = {'nc'};
			obj.paths.diag.SFC_CHL.dim   = {'xyt'};
			obj.paths.diag.SFC_CHL.var   = {'chl'};
			obj.paths.diag.SFC_CHL.zvar  = {[]};
			obj.paths.diag.SFC_CHL.lon   = {'lon'};
			obj.paths.diag.SFC_CHL.lat   = {'lat'};
			obj.paths.diag.SFC_CHL.name  = {'MODIS-Aqua Chlorophyll'};
			obj.paths.diag.SFC_CHL.units = {'$mg$ $C$ $m^{-3}$'};
				
			% NPP products (VGPM, CBPM, CAFE)	
			obj.paths.diag.NPP.file  = {'/data/project2/yangsi/analysis/NPPcode/std-VGPM/std-VGPMnpp_MODIS_clim2002-2018.nc',...
										'/data/project2/yangsi/analysis/NPPcode/CbPM2/CbPM2npp_MODIS_clim2002-2018.nc',...
										'/data/project1/data/MODIS-Aqua/CAFE_NPP/climatology/nppclim_CAFE_MODIS_9km.nc'}; 
			obj.paths.diag.NPP.type  = {'nc','nc','nc'};
			obj.paths.diag.NPP.var   = {'npp','npp','npp'};
			obj.paths.diag.NPP.zvar  = {[],[],[]};
			obj.paths.diag.NPP.dim   = {'xyt','xyt','xyt'};
			obj.paths.diag.NPP.lon   = {'lon','lon','lon'};
			obj.paths.diag.NPP.lat   = {'lat','lat','lat'};
			obj.paths.diag.NPP.name  = {'NPP-VGPM','NPP-CBPM','NPP-CAFE'};
			obj.paths.diag.NPP.units = {'$mg$ $C$ $m^{-2}$ $d^{-1}$','$mg$ $C$ $m^{-2}$ $d^{-1}$','$mg$ $C$ $m^{-2}$ $d^{-1}$'};

			% MLD products (WOCE/NODC/ARGO or ARGO only);
			obj.paths.diag.MLD.file  = {'/data/project1/demccoy/ROMS/validation/MLD/Argo_mixedlayers_monthlyclim_12112019.nc',...
										'/data/project1/demccoy/ROMS/validation/MLD/mld_DR003_c1m_reg2.0.nc'};
			obj.paths.diag.MLD.type  = {'nc','nc'};
			obj.paths.diag.MLD.var   = {'mld_da_mean','mld'};
			obj.paths.diag.MLD.zvar  = {[],[]};
			obj.paths.diag.MLD.dim   = {'txy','xyt'};
			obj.paths.diag.MLD.lon   = {'lon','lon'};
			obj.paths.diag.MLD.lat   = {'lat','lat'};
			obj.paths.diag.MLD.name  = {'Argo Mixed Layer Depth','IFREMER Mixed Layer Depth'}; 
			obj.paths.diag.MLD.units = {'$m$','$m$'};

		end % end methods initPlots

		%--------------------------------------------------------------------------------
		function obj = loadData(obj,vars,varargin)
			% -------------------
			% Method to Load ROMS data. Load either raw data,
			% or data interpolated to standard WOCE depths (z_avg)
			%
			% Usage:
			% - obj = loadData(obj,vars,varargin)
			%
			% Inputs:
			% - vars = variable(s) to load, as a cell array
			%
			% Inputs (varargin):
			% - type    = 'z_avg' or 'raw'
			% - file    = option to load a different file than what was initialized
			% - depth   = can be used with type 'z_avg' to extract specific depth(s) only
			% 
			% Examples:
			% - obj = loadData(obj)
			% - obj = loadData(obj,{'temp','salt'},'type','z_avg');
			% - obj = loadData(obj,{'temp','salt'},'type','z_avg','file','avg_2011');
			% -------------------
			
			% defaults for optional  arguments
			A.type  = ['raw']; % input variable type
			A.file  = [];
			A.depth = []; 
			A      = parse_pv_pairs(A,varargin); % parse method arguments to A
		
			% check that region has been defined
			if isempty(obj.region)
				obj = defineRegion(obj);
			end
	
			% check inputs
			if ~isempty(vars) % vars
				for i = 1:length(vars)
					if ~iscell(vars(i));
						disp(' ');
						disp('vars must be a cell array');
						disp(' '); 
						return
					end
				end
			end			
			if ~isempty(A.file)
				tmpavg         = obj.paths.avg; % file
				tmpzavg        = obj.paths.zavg; % file
				obj.paths.avg  = [obj.paths.simPath,'avg/',A.file,'.nc'];
				obj.paths.zavg = [obj.paths.simPath,'z_avg/z_',A.file,'.nc'];
				avgfile        = [A.file,'.nc'];
			else
				avgfile = ['avg_',num2str(obj.info.Year),'.nc'];
			end
		
			% check depth input
			if ~isempty(A.depth)
				A.type = 'z_avg';
			end

			% if type is defined, select correct choices
			% if not, then define type
			if strcmp(A.type,'raw');
				q0 = 1;
			elseif strcmp(A.type,'z_avg');
				q0 = 2;
				% make z_avg file if it doesn't exist
				if exist([obj.paths.simPath,'z_avg/z_',avgfile]) ~= 2
					disp('---------------------------------');
					disp('Create z_avg Climatology File');
					disp('---------------------------------');
					obj = make_zavg(obj,avgfile);
				end

			elseif isempty(A.type)
				q0 = input('---------------------------\nraw(1) or z_avg(2)\n---------------------------\n>> ');
			end

			% load variables 
			if ~isempty(vars) % user provided inputs, load variables
				for i = 1:length(vars)
					% Check if variable is in ROMS output
					idx = find(ismember(obj.info.var2d,vars{i})==1);
					idx = obj.info.idx2d(idx);
					typ = 2;
					if isempty(idx)
						try
							idx = find(ismember(obj.info.var3d,vars{i})==1);
							idx = obj.info.idx3d(idx);
						catch
							idx = [];
						end
						typ = 3;
					end
					% Check for derived fields			
					if isempty(idx)						
						obj = computeVar(obj,vars(i),'type',A.type);
						% Force end of loop
						continue
					end
					disp('-------------------------');
					disp(['Loading ',vars{i},' data...']);
					if q0 == 1
						tmpdata = ncread(obj.paths.avg,vars{i});
					elseif q0 == 2
						tmpdata = ncread(obj.paths.zavg,vars{i});
					end
					% Reduce to region
					if typ == 2
						tmpdata 		            = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
															  obj.region.lat_lim(1):obj.region.lat_lim(2),:);
						tmpdata                     = tmpdata .* obj.region.mask_rho;
						obj.romsData.(vars{i}).data = tmpdata;
					elseif typ == 3
						tmpdata = tmpdata(obj.region.lon_lim(1):obj.region.lon_lim(2),...
										  obj.region.lat_lim(1):obj.region.lat_lim(2),:,:);
						if q0 == 1 
							tmpdata = tmpdata .* obj.region.mask_rho3d;
						elseif q0 == 2
							tmpdata = tmpdata .* obj.region.mask_rhoz3d;
						end
						obj.romsData.(vars{i}).data = tmpdata;
					end
					disp(['...done!']);
					disp('-------------------------');
					if typ == 2
						ind = find(strcmp(vars{i},obj.info.var2d)==1);
						obj.romsData.(vars{i}).name = obj.info.name2d{ind};
						obj.romsData.(vars{i}).units = obj.info.unit2d{ind};
					elseif typ == 3
						ind = find(strcmp(vars{i},obj.info.var3d)==1);
						obj.romsData.(vars{i}).name = obj.info.name3d{ind};
						obj.romsData.(vars{i}).units = obj.info.unit3d{ind};
					end
				end
			else % no input, show available fields
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
							obj.romsData.(obj.info.var2d{q2(j)}).name  = obj.info.name2d{q2(j)};
							obj.romsData.(obj.info.var2d{q2(j)}).units = obj.info.unit2d{q2(j)};					
						end
					end
				end	
			end

			% Reduce to specific depths if requested	
			if strcmp(A.type,'z_avg') & ~isempty(A.depth);
				zind = find(ismember(obj.grid.z_avg_dep,A.depth)==1);
				for i = 1:length(vars)
					if ndims(obj.romsData.(vars{i}).data)==4
						obj.romsData.(vars{i}).data  = squeeze(obj.romsData.(vars{i}).data(:,:,zind,:));
						obj.romsData.(vars{i}).depth = [obj.grid.z_avg_dep(zind)'];
					end
				end
			end
	
			% Return original avg file
			if ~isempty(A.file)
				obj.paths.avg  = tmpavg;
				obj.paths.zavg = tmpzavg;
			end
		end % end methods loadData

		%--------------------------------------------------------------------------------
		function [obj] = getAvgData(obj,vars,yr_range,varargin)
			% -------------------
			% Get an average field across multiple years
			% 
			% Usage:
			% - obj = getAvgData(obj,vars,yr_range)
			%
			% Inputs:
			% - vars     = cell array of variables to load and average
			% - yr_range = years to average across (i.e. 2000:2005) 
			%
			% Optional (varargin):
			% - lvl      = level(s) to grab (see obj.grid.z_avg_dep)
			%
			% Example:
			% - obj = getAvgData(obj,{'O2','NO3'},2000:2009);
			% - obj = getAvgData(obj,{'O2','NO3'},2000:2009,'lvl',[300 500]);
			% -------------------

			% defaults for optional  arguments
			A.lvl  = [];
			A      = parse_pv_pairs(A,varargin); % parse method arguments to A

			% Loop through years and load variables
			for i = 1:length(vars)
				if ~ismember(vars{i},obj.info.var2d) % 3D variable
					if isempty(A.lvl)
						tmp.(vars{i}) = nan([obj.region.nx obj.region.ny length(obj.grid.z_avg_dep) obj.region.nt length(yr_range)]);
					else
						tmp.(vars{i}) = nan([obj.region.nx obj.region.ny length(A.lvl) obj.region.nt length(yr_range)]);
					end
				else % 2D variable
					tmp.(vars{i}) = nan([obj.region.nx obj.region.ny obj.region.nt]);
				end
				for j = 1:length(yr_range)
					% Get filename
					current_file = ['avg_',num2str(yr_range(j))];
					disp(['...',num2str(yr_range(j))]);
					if ~ismember(vars{i},obj.info.var2d)
						obj = loadData(obj,vars(i),'type','z_avg','file',current_file);
						if isempty(A.lvl)
							tmp.(vars{i})(:,:,:,:,j) = obj.romsData.(vars{i}).data; 
						else
							ind = [];
							for k = 1:length(A.lvl)
								ind(k) = find(obj.grid.z_avg_dep == A.lvl(k));
							end
							tmp.(vars{i})(:,:,:,:,j) = obj.romsData.(vars{i}).data(:,:,ind,:);   
						end
						obj.romsData.(vars{i}).data = [];
					else
						obj                         = loadData(obj,vars(i),'type','raw','file',current_file);
						tmp.(vars{i})(:,:,:,j)      = obj.romsData.(vars{i}).data; 
						obj.romsData.(vars{i}).data = [];
					end
				end
				if ~ismember(vars{i},obj.info.var2d)
					obj.romsData.(vars{i}).data = nanmean(tmp.(vars{i}),5);
				else
					obj.romsData.(vars{i}).data = nanmean(tmp.(vars{i}),4);
				end
			end
		end % end methods getAvgData	

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
			%		      i.e. {'DJF','MAM','JJA','SON'} for seasonal plots
			% - conv    = array of length(vars) to apply conversion factors
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

			% defaults for optional  arguments
			A.vars     = []; % input variable
			A.plttype  = []; % plot type, tmpfig by default
			A.depths   = []; % depths to plot
			A.levels   = []; % levels to plot
			A.time     = []; % time averaging or month to plot
			A.conv     = []; % optional conversion factors
			A.caxis    = []; % optional caxis input (not recommended)
			A          = parse_pv_pairs(A,varargin); % parse method arguments to A

			% month string
			mstr = {'Jan','Feb','Mar','Apr','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

			% depth and level indexes
			lidx = []; didx = [];

			% check that object is initialized
			try; x = obj.info.var3d; catch; obj = init(obj); end

			% check inputs
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
			if isempty(A.conv)
				A.conv = ones(length(A.vars));
			else
			end

			% process depth input
			if ~isempty(A.depths)
				for i = 1:length(A.depths)
					diffd       = abs(A.depths(i) - obj.grid.z_avg_dep);
					ind         = find(diffd == min(diffd));
					didx(i)     = ind;
					dstr{i}     = [num2str(obj.grid.z_avg_dep(ind)),'m'];
				end
			end

			% process levels input
			if ~isempty(A.levels)
				lidx = A.levels;
				for i = 1:length(lidx)
					lstr{i} = ['Level ',num2str(lidx(i))];
				end
			end
			
			% if depth/levels empty, choose
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

			% process variable input
			if ~isempty(A.vars) % user provided input
				for i = 1:length(A.vars)
					if isempty(didx) & isempty(lidx)
						obj = loadData(obj,A.vars(i),'type','raw');
					elseif isempty(didx)
						obj = loadData(obj,A.vars(i),'type','raw');
					elseif isempty(lidx)
						obj = loadData(obj,A.vars(i),'type','z_avg');
					end
				end
			else % no user input
				if isempty(didx)
					obj    = loadData(obj,'type','raw');
				else
					obj    = loadData(obj,'type','z_avg');
				end
				A.vars = fieldnames(obj.romsData);
			end

			% process time input
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

			% get lonbounds/latbounds		
			lonbounds = [floor(obj.region.minlon_rho) ceil(obj.region.maxlon_rho)];
			latbounds = [floor(obj.region.minlat_rho) ceil(obj.region.maxlat_rho)];
			
			% get levels or depths
			if isempty(lidx) & ~isempty(didx)
				zidx = didx; zstr = dstr;
			elseif isempty(didx) & ~isempty(lidx)
				zidx = lidx; zstr = lstr;
			else
				zidx = [0]; zstr = {'SFC'};
			end
			
			% make plot for each variable, depth, and time combo
			pcnt = [0];
			for i = 1:length(A.vars)
				for j = 1:length(zidx)
					for k = 1:length(A.time);
						try % 3D
							data = squeeze(obj.romsData.(A.vars{i}).data(:,:,zidx(j),tidx{k}));
						catch % 2D
							data = squeeze(obj.romsData.(A.vars{i}).data(:,:,tidx{k}));
						end
						
						% apply averaging
						data = nanmean(data,3);
	
						% apply conversion 
						data = data .* A.conv(i);

						% mask coast
						data = data .* obj.region.mask_rho;

						% colorbar limits
						if ~isempty(A.caxis)
							datarng = [A.caxis(1) A.caxis(2)];
						else
							datarng = prctile(data(:),[1 99]);
						end
						[axisLims axisTicks]   = romsMaster.oom_levs(datarng);	
						lvls                   = [axisLims(1):(diff(axisLims)/255):axisLims(2)];	
						data(data>axisLims(2)) = axisLims(2);
						data(data<axisLims(1)) = axisLims(1);
						
						% initiate figure
						fig(i,j,k)            = piofigs('lfig',1);
						set(0,'CurrentFigure',fig(i,j,k));
						[ax(i,j,k)] = map_plot(fig(i,j,k),obj.region.lon_rho,obj.region.lat_rho);
						m_contourf(obj.region.lon_rho,obj.region.lat_rho,data,lvls,'LineStyle','none');
						cb(i,j,k)   = colorbar;
						title([tstr{k},' ',obj.romsData.(A.vars{i}).name,' @ ',zstr{j}]);
						ylabel(cb(i,j,k),obj.romsData.(A.vars{i}).units);
	
						% axis limits
						cb(i,j,k).XLim		 = axisLims;
						cb(i,j,k).XTick      = axisTicks;
						cb(i,j,k).XTickLabel = axisTicks;
						caxis([axisLims]);
					
						% save
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
						
						% rinse and repeat
					end
				end
			end
		end % end methods plotData
	
		%--------------------------------------------------------------------------------
		function obj = sliceROMS(obj,vars,choice,deg,varargin);
			% -------------------
			% Takes 2D depth slice of ROMS along a given latitude or longitude
			% Options are set by user
			% 
			% Usage:
			% - obj = sliceROMS(obj,vars,choice,deg,varargin);
			% 
			% Inputs:
			% - vars   = ROMS variable(s) to slice, as a cell array
			% - choice = 'lon','lat','y', or 'x' 
			%			  (lat/lon slices along a given lat/lon degree)
			%			  (x/y slices along a given x or y index (lon_rho or lat_rho))
			% - deg    = lon/lat degree or x/y index
			%
			% Inputs (varargin):
			% - type     = 'raw' (ROMS levels) or 'z_avg' (WOA-18 depths)
			% - yr_range = years to average across (i.e. 2045:2049) 
			% 
			% Example
			% - obj = sliceROMS(obj,{'temp','salt'},'lon',0);
			% 
			% This will slice temp and salt data at 0 degrees longitude, and will also return
			% sliced temperature (but not salinity) data from the validation dataset
			% -------------------
	
			% Check inputs
			if nargin<4
				disp('Incorrect number of inputs');
				help sliceROMS
				return
			end	

			% Grab user inputs
			A.type     = ['raw'];
			A.yr_range = [];
			A          = parse_pv_pairs(A,varargin);

			% Force z_avg if asking for lon/lat
			if strcmp(choice,'lon') | strcmp(choice,'lat')
				A.type = 'z_avg';
			end

			% Grab longitude/latitude data
			lon  = obj.region.lon_rho;
			lat  = obj.region.lat_rho;

			% Load variables
			if isempty(A.yr_range)
				obj = loadData(obj,vars,'type',A.type);
			else
				obj = getAvgData(obj,vars,A.yr_range);
			end

			% Choose latitude, longitude, x-idx or y-idx
			if strcmp(choice,'lat')
				lonlat = lat;
				dmsn   = obj.region.nx; 
				if strcmp(A.type,'raw')
					nz = obj.region.nz;
				elseif strcmp(A.type,'z_avg');
					nz = length(obj.grid.z_avg_dep); 
				end
				nl     = obj.region.nt; 
				ns     = length(deg); 
				dstr   = ['latitude'];
			elseif strcmp(choice,'lon')
				lonlat = lon; 
				dmsn   = obj.region.ny; 
				if strcmp(A.type,'raw')
					nz = obj.region.nz;
				elseif strcmp(A.type,'z_avg');
					nz = length(obj.grid.z_avg_dep); 
				end
				nl     = obj.region.nt; 
				ns     = length(deg); 
				dstr   = ['longitude'];
			elseif strcmp(choice,'x');
				for v = 1:length(vars)
					obj.romsData.(vars{v}).slice  = squeeze(obj.romsData.(vars{v}).data(deg,:,:,:));
					if strcmp(A.type,'raw')
						obj.slice.depth = squeeze(obj.region.z_r(deg,:,:,:));
					elseif strcmp(A.type,'z_avg');
						obj.slice.depth = permute(repmat(obj.grid.z_avg_dep,1,obj.region.ny,obj.region.nt),[2 3 1]);
					end
					slicelon        = obj.region.lon_rho(deg,:);
					slicelon        = slicelon';
					obj.slice.lon   = repmat(slicelon,1,obj.region.nz);
					slicelat        = obj.region.lat_rho(deg,:);
					slicelat        = slicelat';
					obj.slice.lat   = repmat(slicelat,1,obj.region.nz);
					obj.slice.idx   = deg;
					obj.slice.coord = 'x (lon_rho)';
				end
				% Killing routine here
				return
			elseif strcmp(choice,'y');
				for v = 1:length(vars)
					obj.romsData.(vars{v}).slice  = squeeze(obj.romsData.(vars{v}).data(:,deg,:,:));
					if strcmp(A.type,'raw')
						obj.slice.depth = squeeze(obj.region.z_r(:,deg,:,:));
					elseif strcmp(A.type,'z_avg');
						obj.slice.depth = permute(repmat(obj.grid.z_avg_dep,1,obj.region.nx,obj.region.nt),[2 3 1]);
					end
					slicelon        = obj.region.lon_rho(:,deg);
					obj.slice.lon   = repmat(slicelon,1,obj.region.nz);
					slicelat        = obj.region.lat_rho(:,deg);
					obj.slice.lat   = repmat(slicelat,1,obj.region.nz);
					obj.slice.idx   = deg;
					obj.slice.coord = 'y (lat_rho)'; 
				end
				% Killing routine here
				return
			end

			% Dimensions to fill
			fillmat = [dmsn nz ns];	

			% Copy mask
			mask = obj.region.mask_rho;
			mask(isnan(mask)) = 0;
			mask = repmat(mask,1,1,nz);
		
			% Organize 3D ROMS data and slice
			for ff = 1:length(vars)
				tmpdata{ff} = obj.romsData.(vars{ff}).data;			
			end
			data = cat(5,tmpdata{:});
			clear tmpdata;

			% Take slice of data for each month
			disp(' '); disp(['Slicing ROMS data @ ',num2str(deg),'deg ',dstr,'...']);disp(' ');
			dims = [dmsn nz nl ns obj.region.nx obj.region.ny length(vars)];
			[tmpslice,tmpmask,tmpdepth] = romsMaster.lonlat_slice(dims,lonlat,data,mask,deg,obj.grid.z_avg_dep,fillmat);
				 
			% Get lon/lat data (outside parfor)
			tmpdeg  = NaN(dmsn,ns);
			for i = 1:dmsn
				for j = 1:ns	
					if dmsn == obj.region.nx;
						tmpdeg(i,j)  = interp1(squeeze(lat(i,:)),squeeze(lon(i,:)),deg(j));
					elseif dmsn == obj.region.ny;
						tmpdeg(i,j)  = interp1(squeeze(lon(:,i)),squeeze(lat(:,i)),deg(j));
					end
				end
			end
			tmpdeg = repmat(tmpdeg,1,1,nz); tmpdeg = permute(tmpdeg,[1 3 2]);
			
			% Apply mask
			tmpdepth = squeeze(tmpdepth(:,:,1,1));
			for j = 1:ns
				masktmp         = squeeze(tmpmask(:,:,1,j));
				tmp             = tmpdeg(:,:,j);
				tmp(masktmp==0) = NaN;
				tmpdeg(:,:,j)   = tmp;
				for ff = 1:length(vars)
					for rcrd = 1:nl
						tmp                     = squeeze(tmpslice(:,:,rcrd,j,ff));
						tmp(masktmp==0)         = NaN;
						tmpslice(:,:,rcrd,j,ff) = tmp;
					end
				end
			end

			% Save results
			for ff = 1:length(vars);
				obj.romsData.(vars{ff}).slice  = squeeze(tmpslice(:,:,:,:,ff));
				obj.slice.depth = tmpdepth;
				obj.slice.deg   = tmpdeg;
				if dmsn == obj.region.nx;
						obj.slice.coord = 'latitude';
				elseif dmsn == obj.region.ny;
						obj.slice.coord = 'longitude';
				end
				obj.slice.sect = deg;
				fprintf(['\n']);
			end
		end % end methods sliceROMS

		%--------------------------------------------------------------------------------
		function obj = sliceDiag(obj,vars,choice,deg,varargin);
			% -------------------
			% Takes 2D depth slice of diagnostic data along a given latitude or longitude
			% Options are set by user
			% 
			% Usage:
			% - obj = sliceDiag(obj,vars,choice,deg,varargin);
			% 
			% Inputs:
			% - vars   = diagnostic variable(s) to slice, as a cell array
			%			 (see obj.paths.diag)
			% - choice = 'lon','lat','y', or 'x' 
			%			  (lat/lon slices along a given lat/lon degree)
			%			  (x/y slices along a given x or y index (lon_rho or lat_rho))
			% - deg    = lon/lat degree or x/y index
			%
			% Inputs (varargin):
			% - type = 'raw' (ROMS levels) or 'z_avg' (WOA-18 depths)
			% 
			% Example
			% - obj = sliceDiag(obj,{'temp','salt'},'lon',0,'vali',[1 0]);
			% 
			% This will slice temp and salt data at 0 degrees longitude, and will also return
			% sliced temperature (but not salinity) data from the validation dataset
			% -------------------

			% Grab longitude/latitude data
			lon  = obj.region.lon_rho;
			lat  = obj.region.lat_rho;

			% Choose latitude, longitude
			if strcmp(choice,'lat')
				lonlat = lat;
				dmsn   = obj.region.nx; 
				nz     = length(obj.grid.z_avg_dep); 
				nl     = obj.region.nt; 
				ns     = length(deg); 
				dstr   = ['latitude'];
			elseif strcmp(choice,'lon')
				lonlat = lon; 
				dmsn   = obj.region.ny; 
				nz     = length(obj.grid.z_avg_dep); 
				nl     = obj.region.nt; 
				ns     = length(deg); 
				dstr   = ['longitude'];
			end

			% Get lon/lat data (outside parfor)
			tmpdeg  = NaN(dmsn,ns);
			for i = 1:dmsn
				for j = 1:ns	
					if dmsn == obj.region.nx;
						tmpdeg(i,j)  = interp1(squeeze(lat(i,:)),squeeze(lon(i,:)),deg(j));
					elseif dmsn == obj.region.ny;
						tmpdeg(i,j)  = interp1(squeeze(lon(:,i)),squeeze(lat(:,i)),deg(j));
					end
				end
			end
			tmpdeg = repmat(tmpdeg,1,1,nz); tmpdeg = permute(tmpdeg,[1 3 2]);

			% Dimensions to fill
			fillmat  = [dmsn nz ns];
			tmpdepth = obj.grid.z_avg_dep'.*ones(fillmat);

			% Perform slices of each variable
			for ff = 1:length(vars);
				disp(' '); disp(['Slicing validation ',vars{ff},' data...']);disp(' ');
				% get path and coords for current variable
				curVar  = obj.paths.diag.(vars{ff});
				for i = 1:length(curVar.file); 
					if strcmp(curVar.type{i},'nc');
						tmp.lon   = ncread(curVar.file{i},curVar.lon{i}); tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
						tmp.lat   = ncread(curVar.file{i},curVar.lat{i});
						tmp.depth = ncread(curVar.file{i},curVar.zvar{i});
						if nanmean(tmp.depth(:)) < 0
							tmp.depth = -tmp.depth;
						end
						tmp.data  = ncread(curVar.file{i},curVar.var{i});
						% Apply any factors (eg nMol to uMol)	
						if isfield(curVar,'factor');
							tmp.data = tmp.data*curVar.factor{i};
						end
					elseif strcmp(curVar.type{i},'mat');
						tmp.lon   	       = load(curVar.file{i},curVar.lon{i});
						tmp.lat   	       = load(curVar.file{i},curVar.lat{i});
						tmp.depth 	       = load(curVar.file{i},curVar.zvar{i});
						tmp.data           = load(curVar.file{i},curVar.var{i});
						tmp.lon			   = tmp.lon.(curVar.lon{i});
						tmp.lat            = tmp.lat.(curVar.lat{i});
						tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
						tmp.depth          = tmp.depth.(curVar.zvar{i});
						tmp.data           = tmp.data.(curVar.var{i});
						% Apply any factors (eg nMol to uMol)	
						if isfield(curVar,'factor');
							tmp.data = tmp.data*curVar.factor{i};
						end
					end

					% make lon/lat data gridded if it isnt setup that way
					if sum(size(tmp.lon(:,:))==1) > 0
						[tmp.latg,tmp.long,tmp.depthg] = meshgrid(tmp.lat,tmp.lon,tmp.depth);
					end
					
					% go through each requested slice
					for j = 1:ns
						% go through each time record
						for rcrd = 1:nl
							% clear variables to fill
							tmp.latr   = []; 
							tmp.lonr   = [];
							tmp.depthr = [];
							tmp.datar  = [];
							% record progress
							fprintf([num2str(rcrd),'...']);
							% Get reduced matrix to feed into scatteredInterpolant
							if strcmp(choice,'lat');
								idx = find(abs([tmp.lat-deg(j)]) == min(abs([tmp.lat-deg(j)])));
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
							elseif strcmp(choice,'lon');
							% Get reduced matrix to feed into scatteredInterpolant
								idx = find(abs([tmp.lon-deg(j)]) == min(abs([tmp.lon-deg(j)])));
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
							% build interpolant and interpolate
							if strcmp(choice,'lat');
								tmp.lonr      = double(tmp.lonr(:));
								tmp.depthr    = double(tmp.depthr(:));
								tmp.datar     = double(tmp.datar(:)); 
								F             = scatteredInterpolant(tmp.lonr(~isnan(tmp.datar)),...
												tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
								tmp.out{rcrd} = F(tmpdeg(:,:,j),tmpdepth(:,:,j));
							elseif strcmp(choice,'lon');
								tmp.latr      = double(tmp.latr(:));
								tmp.depthr    = double(tmp.depthr(:));
								tmp.datar     = double(tmp.datar(:));
								F             = scatteredInterpolant(tmp.latr(~isnan(tmp.datar)),...
												tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
								tmp.out{rcrd} = F(tmpdeg(:,:,j),tmpdepth(:,:,j));
							end
						end % end rcrd-loop
						fprintf(['\n']);
						% Save results
						obj.diagData.(vars{ff})(i).slice(:,:,:,j)  = cat(ndims(tmp.out{1})+1,tmp.out{:});
						obj.diagData.(vars{ff})(i).units = curVar.units{i};
						obj.diagData.(vars{ff})(i).name  = curVar.name{i};
						obj.slice.deg(:,:,j) = tmpdeg(:,:,j);
					end % end j-loop
					if dmsn == obj.region.nx;
						obj.slice.coord = 'latitude';
					elseif dmsn == obj.region.ny;
						obj.slice.coord = 'longitude';
					end
					obj.slice.depth = tmpdepth;
					obj.slice.sect  = deg;
				end % end i-loop
			end % end var-loop
		end % end method sliceDiag

		%--------------------------------------------------------------------------------
		function obj = computeVar(obj,vars,varargin)
			% ------------------
			% Computes additional fields like spiciness, potential density etc.
			% 
			% List of available variables:
			% - pres     (pressure)
			% - sigma    (density - 1000 kg/m3)
			% - bvf      (Brunt Vaisala frequency)
			% - sigma0   (potential density)
			% - pv       (potential vorticity)
			% - SSH      (sea-surface height, corrected via zeta) 
			% - ws       (wind-stress)
			% - wsc      (wind-stress curl)
			% - NPP      (net primary production)
			% - nstar    (N* via NO3 - 16*PO4 + 2.9);
			% - MLD      (mixed layer depth)
			% - SFC_CHL  (surface integrated chla, to compare against satellite)
			%
			% Usage:
			% - obj = computeVar(obj,vars,varargin)
			%
			% Inputs:
			% - var = variable(s) to compute, as a cell array	
			%
			% Inputs (varargin)
			% - type = 'raw' (default) or 'z_avg'
			%
			% Example:
			% - obj = computeVar(obj,{'sigma0'},'type','z_avg');
			% ------------------

			% process inputs
			A.type    = 'raw';
			A         = parse_pv_pairs(A,varargin);

			% dims for vertical computations
			nx = obj.region.nx;
			ny = obj.region.ny;
			nt = obj.region.nt;
			if strcmp(A.type,'raw');
				nz = obj.region.nz;
			elseif strcmp(A.type,'z_avg');
				nz = length(obj.grid.z_avg_dep);
			end

			% process requests
			for i = 1:length(vars)
				% pressure calc
				if strcmp(vars{i},'pres') % & ~isfield(obj.romsData,'pres');
					disp('Calculating sea pressure');
					if strcmp(A.type,'raw');
						obj.romsData.pres.data = sw_pres(-obj.region.z_r,repmat(obj.region.lat_rho,1,1,obj.region.nz,obj.region.nt));
						obj.romsData.pres.data = obj.romsData.pres.data .* obj.region.mask_rho3d;
					elseif strcmp(A.type,'z_avg');
						tmpdep = repmat(obj.grid.z_avg_dep,1,obj.region.nx,obj.region.ny,obj.region.nt);
						tmpdep = permute(tmpdep,[2 3 1 4]);
						obj.romsData.pres.data = sw_pres(tmpdep,repmat(obj.region.lat_rho,1,1,length(obj.grid.z_avg_dep),obj.region.nt));
						obj.romsData.pres.data = obj.romsData.pres.data .* obj.region.mask_rhoz3d;
					end
					obj.romsData.pres.name  = 'averaged Pressure';
					obj.romsData.pres.units = '$dbar$';
				% rho calc
				elseif strcmp(vars{i},'sigma')
					disp('Calculating sigma');
					obj = computeVar(obj,{'pres'},'type',A.type);
					if isfield(obj.romsData,'temp');
						tt = 1;
						origtemp = obj.romsData.temp;
					else
						tt = 0;
					end
					if isfield(obj.romsData,'salt');
						ss = 1;
						origsalt = obj.romsData.salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},'type',A.type);
					tmptemp = sw_temp(obj.romsData.salt.data(:),obj.romsData.temp.data(:),obj.romsData.pres.data(:),0); 
					tmpsigma = sw_dens(obj.romsData.salt.data(:),tmptemp,obj.romsData.pres.data(:));
					obj.romsData.sigma.data = reshape(tmpsigma,size(obj.romsData.temp.data))-1000;
					if strcmp(A.type,'raw');
						obj.romsData.sigma.data = real(obj.romsData.sigma.data) .* obj.region.mask_rho3d;
					elseif strcmp(A.type,'z_avg');
						obj.romsData.sigma.data = real(obj.romsData.sigma.data) .* obj.region.mask_rhoz3d;
					end
					obj.romsData.sigma.name  = 'averaged Density';
					obj.romsData.sigma.units = '$kg$ $m^{-3}$';
					if tt == 1
						obj.romsData.temp = origtemp;
					else
						obj.romsData.temp = [];
					end
					if ss == 1
						obj.romsData.salt = origsalt;
					else
						obj.romsData.salt = [];
					end
				% bvf calc
				elseif strcmp(vars{i},'bvf') | strcmp(vars{i},'pv');
					disp('Calculating bvf (Brunt Vaisala)');
					obj = computeVar(obj,{'pres'},'type',A.type);
					if isfield(obj.romsData,'temp');
						tt = 1;
						origtemp = obj.romsData.temp;
					else
						tt = 0;
					end
					if isfield(obj.romsData,'salt');
						ss = 1;
						origsalt = obj.romsData.salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},'type',A.type);
					tmptemp = reshape(sw_temp(obj.romsData.salt.data(:),obj.romsData.temp.data(:),...
											  obj.romsData.pres.data(:),0),size(obj.romsData.salt.data)); 
					if strcmp(A.type,'raw')
						tmplat  = repmat(obj.region.lat_rho,1,1,obj.region.nz,obj.region.nt);
					elseif strcmp(A.type,'z_avg');
						tmplat  = repmat(obj.region.lat_rho,1,1,length(obj.grid.z_avg_dep),obj.region.nt);
					end
					% Permute to get depth in front
					tmptemp = permute(obj.romsData.temp.data,[3 1 2 4]);
					tmpsalt = permute(obj.romsData.salt.data,[3 1 2 4]);
					tmppres = permute(obj.romsData.pres.data,[3 1 2 4]);
					tmplat  = permute(tmplat,[3 1 2 4]);
					% Reshape to MxN
					dims = size(tmptemp);
					tmptemp = reshape(tmptemp,dims(1),dims(2)*dims(3)*dims(4)); 
					tmpsalt = reshape(tmpsalt,dims(1),dims(2)*dims(3)*dims(4));	
					tmppres = reshape(tmppres,dims(1),dims(2)*dims(3)*dims(4));
					tmplat  = reshape(tmplat,dims(1),dims(2)*dims(3)*dims(4));
					[tmpbvf,tmpq,tmppav] = sw_bfrq(tmpsalt,tmptemp,tmppres,tmplat);
					% Interpolate
					bvf_rho = 0.5*(tmpbvf(1:end-1,:) + tmpbvf(2:end,:));
					pv_rho  = 0.5*(tmpq(1:end-1,:) + tmpq(2:end,:));
					% Add NaN data above and below
					bvf_rho = [nan(1,size(bvf_rho,2));bvf_rho;nan(1,size(bvf_rho,2))];
					pv_rho  = [nan(1,size(pv_rho,2));pv_rho;nan(1,size(pv_rho,2))];
					% Remake
					bvf_rho = reshape(bvf_rho,dims(1),dims(2),dims(3),dims(4));
					bvf_rho = permute(bvf_rho,[2 3 1 4]);
					pv_rho  = reshape(pv_rho,dims(1),dims(2),dims(3),dims(4));
					pv_rho  = permute(pv_rho,[2 3 1 4]);
					obj.romsData.bvf.data  = bvf_rho;
					obj.romsData.bvf.name  = 'averaged Brunt Vaisala frequency';
					obj.romsData.bvf.units = '$s^{-2}$';
					obj.romsData.pv.data   = pv_rho;
					obj.romsData.pv.name   = 'averaged Potential Vorticity';
					obj.romsData.pv.units  = '$m^{-1}$ $s^{-1}$';
					if tt == 1
						obj.romsData.temp = origtemp;
					else
						obj.romsData.temp = [];
					end
					if ss == 1
						obj.romsData.salt = origsalt;
					else
						obj.romsData.salt = [];
					end
					if strcmp(A.type,'raw');
						obj.romsData.bvf.data = obj.romsData.bvf.data .* obj.region.mask_rho3d;
						obj.romsData.pv.data  = obj.romsData.pv.data .* obj.region.mask_rho3d;
					elseif strcmp(A.type,'z_avg');
						obj.romsData.bvf.data = obj.romsData.bvf.data .* obj.region.mask_rhoz3d;
						obj.romsData.pv.data  = obj.romsData.pv.data .* obj.region.mask_rhoz3d;
					end
				% sigma0
				elseif strcmp(vars{i},'sigma0')
					disp('Calculating sigma0');
					obj = computeVar(obj,{'pres'},'type',A.type);
					if isfield(obj.romsData,'temp');
						tt = 1;
						origtemp = obj.romsData.temp;
					else
						tt = 0;
					end
					if isfield(obj.romsData,'salt');
						ss = 1;
						origsalt = obj.romsData.salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},'type',A.type);
					tmp_lon      = repmat(obj.region.lon_rho,1,1,nz,nt);
					tmp_lat      = repmat(obj.region.lat_rho,1,1,nz,nt);
					tmp_salt     = obj.romsData.salt.data;
					tmp_pres     = obj.romsData.pres.data;
					tmp_salt_abs = romsMaster.getSA(tmp_salt,tmp_pres,tmp_lon,tmp_lat);
					tmp_theta    = gsw_CT_from_pt(tmp_salt_abs(:),obj.romsData.temp.data(:));
					tmp_sigma0   = gsw_sigma0(tmp_salt_abs(:),tmp_theta(:));
					if strcmp(A.type,'z_avg')
						tmp_sigma0   = reshape(tmp_sigma0,obj.region.nx,obj.region.ny,length(obj.grid.z_avg_dep),obj.region.nt);
					elseif strcmp(A.type,'raw');
						tmp_sigma0   = reshape(tmp_sigma0,obj.region.nx,obj.region.ny,obj.region.nz,obj.region.nt);
					end
					obj.romsData.sigma0.data  = tmp_sigma0;
					if strcmp(A.type,'raw');
						obj.romsData.sigma0.data = obj.romsData.sigma0.data .* obj.region.mask_rho3d;
					elseif strcmp(A.type,'z_avg');
						obj.romsData.sigma0.data = obj.romsData.sigma0.data .* obj.region.mask_rhoz3d;
					end
					obj.romsData.sigma0.name  = 'averaged Potential Density';
					obj.romsData.sigma0.units = '$kg$ $m^{-3}$';
					if tt == 1
						obj.romsData.temp = origtemp;
					else
						obj.romsData.temp = [];
					end
					if ss == 1
						obj.romsData.salt = origsalt;
					else
						obj.romsData.salt = [];
					end
				% nstar
				elseif strcmp(vars{i},'nstar');
					obj = loadData(obj,{'NO3','PO4'},'type',A.type);
					tmpnstar = obj.romsData.NO3.data - 16.*obj.romsData.PO4.data + 2.9;
					obj.romsData.nstar.data = tmpnstar;
					if strcmp(A.type,'raw');
						obj.romsData.nstar.data = obj.romsData.nstar.data .* obj.region.mask_rho3d;
					elseif strcmp(A.type,'z_avg');
						obj.romsData.nstar.data = obj.romsData.nstar.data .* obj.region.mask_rhoz3d;
					end
					obj.romsData.nstar.name = 'averaged N*';
					obj.romsData.nstar.units = '$mmol$ $N$ $m^{-3}$';
				% NPP
				elseif strcmp(vars{i},'NPP') | strcmp(vars{i},'npp');
					obj     = loadData(obj,{'TOT_PROD'},'type','raw');
					tmpprod = obj.romsData.TOT_PROD.data;
					tmpHz   = obj.region.Hz; 
					tmpNPP  = squeeze(nansum(tmpprod.*tmpHz,3)).*3600*24*12;
					obj.romsData.NPP.data  = tmpNPP;
					obj.romsData.NPP.data  = obj.romsData.NPP.data .* obj.region.mask_rho;
					obj.romsData.NPP.name  = 'averaged Net Primary Production (NPP)';
					obj.romsData.NPP.units = '$mg$ $C$ $m^{-2}$ $d^{-1}$'; 
				% SSH	
				elseif strcmp(vars{i},'SSH') | strcmp(vars{i},'ssh');
					obj = loadData(obj,{'zeta'},'type','raw');
					obj = loadDiag(obj,{'SSH'},0);
					slacorr = nanmedian(obj.diagData.SSH.data(:)) - nanmedian(obj.romsData.zeta.data(:));
					disp(' '); disp(['Adding correction of ',num2str(slacorr),'m to ROMS SSH']);
					obj.romsData.SSH.data = obj.romsData.zeta.data + slacorr;
					obj.romsData.SSH.name = 'averaged sea-surface height';
					obj.romsData.SSH.units = '$m$'; 
				% Wind Stress or Wind Stress Curl
				elseif strcmp(vars{i},'WS') | strcmp(vars{i},'ws') | strcmp(vars{i},'WSC') | strcmp(vars{i},'wsc');
					obj    = loadData(obj,{'sustr','svstr'},'type','raw');
					tmpu   = obj.romsData.sustr.data;
					tmpv   = obj.romsData.svstr.data;
					tmpang = obj.region.angle;
					[tmpws,tmpwsc] = romsMaster.WindStress(tmpu,tmpv,obj.region.lon_rho,obj.region.lat_rho,tmpang);
					obj.romsData.ws.data  = tmpws;
					obj.romsData.wsc.data = tmpwsc;
					obj.romsData.ws.data  = obj.romsData.ws.data .* obj.region.mask_rho;
					obj.romsData.wsc.data = obj.romsData.wsc.data .* obj.region.mask_rho;
					obj.romsData.ws.name  = 'averaged wind-stress';
					obj.romsData.wsc.name = 'averaged wind-stress curl';
					obj.romsData.ws.units = '$N$ $m^{-2}$';
					obj.romsData.wsc.units = '$N$ $m^{-2}$';
				% Mixed-layer depth (MLD)
				elseif strcmp(vars{i},'MLD') | strcmp(vars{i},'mld');
					obj  = computeVar(obj,{'sigma0'},'type','z_avg');
					zind = find(obj.grid.z_avg_dep==10); 
					d10  = squeeze(obj.romsData.sigma0.data(:,:,zind,:));
					tmpmld = nan(obj.region.ndim_xyt);
					tmpdepth = permute(repmat(obj.grid.z_avg_dep,[1 obj.region.ndim_xy]),[2 3 1]);
					for t = 1:obj.region.nt;
						tmpdens = squeeze(obj.romsData.sigma0.data(:,:,:,t));
						ind1    = tmpdens>([d10(:,:,t)+0.03]);
						ind2    = isnan(tmpdens);
						glevs   = squeeze(sum(ind1,3));
						nlevs   = squeeze(sum(ind2,3));
						levs    = length(obj.grid.z_avg_dep) - (glevs + nlevs) + 1;
						for x = 1:obj.region.nx
							for y = 1:obj.region.ny
								if levs(x,y) == 1
									tmpmld(x,y,t) = NaN;
								else
									tmpmld(x,y,t) = tmpdepth(x,y,levs(x,y));
								end
								if isnan(tmpmld(x,y,t))==1 & levs(x,y)>1
									disp('bad');
								end
							end
						end
					end	

					% Save ROMS data (make positive)
					obj.romsData.MLD.data = tmpmld;	
					obj.romsData.MLD.data = obj.romsData.MLD.data .* obj.region.mask_rho;
					obj.romsData.MLD.name = 'averaged mixed-layer depth';
					obj.romsData.MLD.units = '$m$';
				% ChlA
				elseif strcmp(vars{i},'SFC_CHL') | strcmp(vars{i},'sfc_chl');
					obj     = loadData(obj,{'TOT_CHL'},'type','z_avg');
					ind     = find(obj.grid.z_avg_dep<50);
					tmpchla = obj.romsData.TOT_CHL.data;
					tmpchla = squeeze(nanmean(tmpchla(:,:,ind,:),3));
					obj.romsData.SFC_CHL.data  = tmpchla .* obj.region.mask_rho;
					obj.romsData.SFC_CHL.name  = 'averaged surface chlA';
					obj.romsData.SFC_CHL.units = '$mg$ $chlA$ $m^{-3}$';
				else
					disp([vars{i},' is already, or cant be, calculated']);
				end
			end
		end % end methods computeVar
	
		%--------------------------------------------------------------------------------
		function obj = ipslice(obj,vars,ip)
			% ------------------
			% Slices ROMS along a given isopycnal
			% 
			% Usage:
			% - obj = ipslice(obj,vars,ip);
			% 
			% Inputs:
			% - vars = cell array of ROMS variables to slice
			% - ip   = isopycnal(s) to slice along
			%
			% Example:
			% - obj = ipslice(obj,{'O2'},26.5);
			% -------------------	

			% Compute sigma0
			try; obj.romsData.sigma0.data;
			catch; obj = computeVar(obj,{'sigma0'},'type','raw');
			end

			% Load variables
			for i = 1:length(vars)
				if ismember(vars{i},obj.info.var2d) | ismember(vars{i},obj.info.var3d);
					obj = loadData(obj,vars(i),'type','raw');
				else
					obj = computeVar(obj,vars(i),'type','raw');
				end
			end

			% Interpolate variable to ip level
			for i = 1:length(vars)
				vnew = nan(obj.region.nx,obj.region.ny,length(ip),obj.region.nt);
				for z = 1:length(ip);
					for t = 1:obj.region.nt
						% Get variable for time t
						tmpvar = squeeze(obj.romsData.(vars{i}).data(:,:,:,t));
						tmpip  = squeeze(obj.romsData.sigma0.data(:,:,:,t));

						% Get sizes
						[lx,ly,lz]=size(tmpip);

						% Find level where sigma0 > ip
						a=tmpip>ip(z);
						levs=squeeze(sum(a,3));
						levs(levs==lz)=NaN;
						mask=levs./levs;
						mask(isnan(mask)) = 0;
						
						% Index
						ip1 = [];
						ip2 = [];
						for x = 1:lx
							for y = 1:ly
								if levs(x,y)>0
									ip1(x,y) = tmpip(x,y,levs(x,y)); % tmpip greater than ip
									ip2(x,y) = tmpip(x,y,levs(x,y)+1); % tmpip less than ip
									v1(x,y)  = tmpvar(x,y,levs(x,y)); % var at tmpip greater than ip
									v2(x,y)  = tmpvar(x,y,levs(x,y)+1); % var at tmpip less than ip
								else
									ip1(x,y) = NaN;
									ip2(x,y) = NaN;
									v1(x,y)  = NaN;
									v2(x,y)  = NaN;
								end
							end
						end

						% Mask data
						ip1(mask==0) = NaN;
						ip2(mask==0) = NaN;
						v1(mask==0)  = NaN;
						v2(mask==0)  = NaN;

						% Perform linear interpolation
						vnew(:,:,z,t) = [v2.*(ip1-ip(z)) + v1.*(ip(z)-ip2)]./[(ip1-ip2)];
					end
				end
				% Replace data
				obj.romsData.(vars{i}).data    = [];
				obj.romsData.(vars{i}).data	   = vnew;
				obj.romsData.(vars{i}).ipslice = ip;
			end

		end % end methods ipslice

		%--------------------------------------------------------------------------------
		function obj = equatorUcmp(obj,sect)
			% -----------------------
			% Compares equatorial current structure with Cravette2017 results or Johnson2002
			% Returns depth slices of u-velocity along specific longitude bands
			%
			% Usage:
			% - obj = equatorUcmp(obj,sect)
			%
			% Inputs:
			% - sect: Section to compare
			% choose from the following (as a string) for Cravatte 2017
			%		  '120W_90W', '135E_160E', '160E_180', '160W_120W', '179W_160W'
			% choose from the following (as an integer) for Johnson 2002
			%		  143, 156, 165, 180, 190, 205, 220, 235, 250, 265
			%
			% -----------------------
			
			% Check input
			if isnumeric(sect)
				% Load Johnson(2002) mean zonal currents from section
				fname = ['/data/project1/data/equatorial_currents_johnson/meanfit_m.cdf'];
				datu  = ncread(fname,'UM');
				depu  = ncread(fname,'ZDEP1_50');
				latu  = ncread(fname,'YLAT11_101');
				lonu  = ncread(fname,'XLON');
				
				% Grab section
				sect_avg = sect;
				datu     = squeeze(datu(find(lonu == sect),:,:));
				datu     = datu.*100;
				diagname = 'Johnson(2002) u-velocity';
			else
				% Load Cravette(2017) mean zonal currents from section
				fname = ['/data/project1/data/Tropical_Currents_Cravatte_2017/Mean_zonal_currents_',sect,'.cdf'];
				datu  = ncread(fname,'U_SADCPD');
				depu  = ncread(fname,'DEPTH');
				if length(depu) ~= size(datu,2)
					depu = ncread(fname,'AXREG1_68');
				end
				latu  = ncread(fname,'LATI');
				
				% Get average longitude based on sect
				if strcmp(sect,'120W_90W');
					sect_avg = [255];
				elseif strcmp(sect,'135E_160E');
					sect_avg = [147.5];
				elseif strcmp(sect,'160E_180');
					sect_avg = [170];
				elseif strcmp(sect,'160W_12OW');
					sect_avg = [220];
				elseif strcmp(sect,'179E_160W');
					sect_avg = [190];
				end

				% Correct 120W_90W
				if sect_avg == 255
					[a,b]        = size(datu);
					extra_depths = [0; -10; -20];
					fill_empty   = ones(a,length(extra_depths));
					fill_data    = nan(size(fill_empty));
					% Add
					depu = [extra_depths;depu];
					datu = [fill_data datu];
				end
				diagname = 'Cravatte(2017) u-velocity';
			end
			% Grid that shit
			[depu,latu] = meshgrid(-depu,latu);
			
			% Get 3D u velocity data
			romsu  = ncread(obj.paths.zavg,'u',[obj.region.lon_lim(1),obj.region.lat_lim(1),1,1],...
											   [diff(obj.region.lon_lim)+1,diff(obj.region.lat_lim)+1,inf,inf]);
			romsln = obj.region.lon_rho;
			romslt = obj.region.lat_rho;
			romsdp = obj.grid.woa1p0.depth;
			
			% Take slice of data for each month
			dmsn = obj.region.ny;
			nl   = obj.region.nt;
			nz   = obj.grid.woa0p25.depth; nzz = length(nz);
			disp(' '); disp(['Slicing ROMS u data']);disp(' ');
			tmpslice = NaN(dmsn,nzz,nl);
			%parfor mnth = 1:nl
			for mnth = 1:nl
				fprintf([num2str(mnth),'...']);
				% Only get data for that month
				tmpdata = squeeze(romsu(:,:,:,mnth));
				% - Interp to lon/lat line
				for i = 1:dmsn
					% - Fill tmpslice			
					for z = 1:nzz;
						try
							tmpslice(i,z,mnth) = interp1(squeeze(romsln(:,i)),squeeze(tmpdata(:,i,z)),sect_avg);
						catch
							disp('Cant interpolate ROMS at this longitude, try again');
							return
						end
					end
				end
			end

			% - Set up grid
			try
				for i = 1:dmsn
					tmpdeg(i) = interp1(squeeze(romsln(:,i)),squeeze(romslt(:,i)),sect_avg);
				end
				tmpdepth = -nz;
				[tmpdepth,tmpdeg] = meshgrid(tmpdepth,tmpdeg);
			catch
				disp('Cant interpolate ROMS at this longitude, try again');
				return
			end
			
			% Interpolate results to validation grid	
			tmpdepth = tmpdepth(:);
			tmpdeg   = tmpdeg(:);
			for mnth = 1:12
				tmpu        = tmpslice(:,:,mnth);
				idx         = find(isnan(tmpu)==0);
				F           = scatteredInterpolant(tmpdeg(idx),tmpdepth(idx),tmpu(idx));
				gridu{mnth} = F(latu,depu);
			end

			% Save results in format for slicePlot 
			idx = find(strcmp(obj.info.var3d,'u')==1);
			obj.romsData.u.slice  = cat(3,gridu{:});
			obj.romsData.u.name  = obj.info.name3d{idx};
			obj.romsData.u.units  = obj.info.unit3d{idx};
			obj.diagData.u.slice  = datu./100;
			obj.diagData.u.name   = diagname;
			obj.diagData.u.units  = obj.info.unit3d{idx};
			obj.slice.deg   = latu;
			obj.slice.depth = depu;
			if nanmean(depu) < 0
				obj.slice.depth = -depu;
			end
			obj.slice.coord = 'longitude';
			obj.slice.sect  = sect_avg;
		end % end method equatorUcmp

		%--------------------------------------------------------------------------------
		function obj = getProfile(obj,vars,lon,lat)
			% ------------------
			% Loads profile data at the nearest lon/lat point 
			%
			% Usage:
			%	- obj = getProfile(obj,vars,lon,lat);
			%
			% Inputs:
			%	- vars  = variables to load, as a cell array
			%	- lon   = vector of longitude points
			%	- lat   = vector of latitude points
			%
			% Inputs (varargin):
			%	- yr_range = range of years to load and average (i.e. 2045:2049)
			%
			% Example:
			%	- obj = getProfile(obj,{'temp','salt','O2','NO3'},[250 250],[-15 -20]);
			% -------------------

			% process inputs
			A.yr_range = [];
			A          = parse_pv_pairs(A,varargin);	

			% Get indices of nearest point
			lon_idx = [];
			lat_idx = [];
			for i = 1:length(lon);
				% Get index of nearest lat/lon grid cell
				all_dist                = distance(lat(i),lon(i),obj.region.lat_rho,obj.region.lon_rho);
				[lon_idx(i),lat_idx(i)] = find(all_dist == min(all_dist(:)));
			end
			   
			% Load data
			if isempty(A.yr_range)
				obj = loadData(obj,vars,'type','z_avg');
			else
				obj = getAvgData(obj,vars,A.yr_range);
			end
				
			% Go through all variables and lon/lats
			for v = 1:length(vars)
				disp(['...grabbing ',vars{v},' profile(s)']);
				tmp.data = obj.romsData.(vars{v}).data;	
				% Save profile data
				for i = 1:length(lon);
					obj.romsData.(vars{v}).profile(i,:,:)  = squeeze(tmp.data(lon_idx(i),lat_idx(i),:,:));
					obj.profile.depth	= obj.grid.z_avg_dep; 
					obj.profile.lon(i)  = obj.region.lon_rho(lon_idx(i),lat_idx(i));
					obj.profile.lat(i)  = obj.region.lat_rho(lon_idx(i),lat_idx(i));
					obj.profile.xidx(i) = lon_idx(i);
					obj.profile.yidx(i) = lat_idx(i);
				end
			end
		end % end method getProfile

		%--------------------------------------------------------------------------------
		function obj = loadDiag(obj,vars,depth,varargin)
			% ------------------
			% Loads and interpolates monthly 2D validation data to the ROMS grid
			% See initDiags for available fields or type (obj.paths.diag)
			%
			% Usage:
			% - obj = loadDiag(obj,vars,depth,varargin);
			%
			% Inputs:
			% - vars  = validation variable(s) to load, as a cell array; paths set in initDiag
			% - depth = depths to interpolate data, will round to nearest z_avg_dep
			%			for surface data, use depth = 0
			%
			% Inputs (varargin)
			% - outer = degrees of lon/lat to pad interpolant, default = 3
			%
			% Example:
			% - obj = loadDiag(obj,vars,depth);
			% -------------------

			% Optional arguments
			A.outer = [3]; % interpolant padding (degrees)
			A=parse_pv_pairs(A,varargin);

			% Check that object is initialized(init) and initDiag has run
			obj = initDiag(obj);
			
			% Check inputs
			diagfields = fieldnames(obj.paths.diag); 
			for i = 1:length(vars)
				if ~strcmp(vars{i},diagfields) & ~strcmp(vars{i},upper(diagfields));
					disp(' ');
					disp([vars{i},' is not a diagnostic variable']);
					disp(' '); return
				end
			end
			if min(depth) < 0 | max(depth) > max(obj.grid.z_avg_dep)
				disp(' ');
				disp('Check depth input');
				disp(' '); return
			end
			for i = 1:length(depth)
				diffd = abs(depth(i) - obj.grid.z_avg_dep);
				ind   = find(diffd == min(diffd));
				depth(i) = obj.grid.z_avg_dep(ind);
			end
			
			% Get temporary output longitude and latitude
			tmp.outlon = obj.region.lon_rho;
			tmp.outlat = obj.region.lat_rho;

			% Process each variable
			for i = 1:length(vars)		
				fprintf(['\n Processing ', vars{i}]);

				% Get path and coords for current variable
				curVar  = obj.paths.diag.(vars{i});

				% Go through multiple files, if they exist
				for j = 1:length(curVar.file)
					
					% Get coordinates
					if strcmp(curVar.type{j},'nc');
						tmp.lon = romsMaster.lon360(ncread(curVar.file{j},curVar.lon{j}));
						tmp.lat = ncread(curVar.file{j},curVar.lat{j});
					elseif strcmp(curVar.type{j},'mat');
						tmp.lon = load(curVar.file{j},curVar.lon{j});
						tmp.lat = load(curVar.file{j},curVar.lat{j});
						tmp.lon = tmp.lon.(curVar.lon{j});
						tmp.lat = tmp.lat.(curVar.lat{j});
						tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
					end

					% Make lon/lat data gridded if it isnt setup that way
					if sum(size(tmp.lon(:,:))==1) > 0
						[tmp.lat,tmp.lon] = meshgrid(tmp.lat, tmp.lon);
					end;
					tmp = romsMaster.struct2double(tmp);
				
					% Get indeces for reduced domain interpolation
					% Note: this breaks if ROMS boundary longitude is close to 0 or 360
					idx = find(tmp.lon(:) > obj.region.minlon_rho-A.outer ...
							 & tmp.lon(:) < obj.region.maxlon_rho+A.outer ...
							 & tmp.lat(:) > obj.region.minlat_rho-A.outer ...
							 & tmp.lat(:) < obj.region.maxlat_rho+A.outer);

					% Initialize interpolated output
					obj.diagData.(vars{i})(j).data = nan(obj.region.nx,obj.region.ny,length(depth),obj.region.nt);

					% Interpolate data on ROMS coords for each month
					lvls = length(depth);
					fprintf('\n month:');
					for k = 1:12
						tmpdata = []; tmpdepth = [];
						fprintf([num2str(k),'...']);

						% Get data
						if strcmp(curVar.type{j},'nc')
							if strcmp(curVar.dim{j},'xyt')
								tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
											[1,1,k],[inf,inf,1]));
								tmpdepth = 0;
							elseif strcmp(curVar.dim{j},'txy');
								tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
											[k,1,1],[1,inf,inf]));
								tmpdepth = 0;
							elseif strcmp(curVar.dim{j},'xyzt')
								tmpdepth = ncread(curVar.file{j},curVar.zvar{j});
								tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
									   [1,1,1,k],[inf,inf,inf,1]));
							end
						elseif strcmp(curVar.type{j},'mat')
							if strcmp(curVar.dim{j},'xyt');
								tmpdata  = load(curVar.file{j},curVar.var{j});
								tmpdata  = tmpdata.(curVar.var{j})(:,:,k);
								tmpdepth = 0;
							elseif strcmp(curVar.dim{j},'xyzt')
								tmpdepth = load(curVar.file{j},curVar.zvar{j});
								tmpdepth = tmpdepth.(curVar.zvar{j});
								tmpdata  = load(curVar.file{j},curVar.var{j});
								tmpdata  = tmpdata.(curVar.var{j})(:,:,:,k);
							end
						end
		
						% Apply any factors (eg nMol to uMol)	
						if isfield(curVar,'factor');
							tmpdata = tmpdata*curVar.factor{j};
						end

						% Build interpolant and interpolate over all depths
						for z = 1:lvls
							data  = tmpdata;
							zdepth = tmpdepth; 
							if nanmean(depth) < 0
								depth = -depth;
							end
							zind  = find(abs([zdepth-depth(z)]) == min(abs([zdepth-depth(z)])));
							if ~strcmp(curVar.dim{j},'xyt')
								data  = data(:,:,zind);
							end
							if strcmp(curVar.dim{j},'xyt') & z > 1
								disp(['No z-data for ',vars{i},', skipping']);
								continue
							end
							% Interpolate
							F = scatteredInterpolant(double(tmp.lon(idx)),double(tmp.lat(idx)),...
													 double(data(idx)),'linear','nearest');
							tmpout{z,k} = F(double(tmp.outlon),double(tmp.outlat));
							tmpout{z,k}(isnan(obj.region.mask_rho)) = nan;
						end % end z-loop
					end % end k-loop
					% Save interpolated data
					for ll = 1:length(depth);
						obj.diagData.(vars{i})(j).data(:,:,ll,:)  = cat(3, tmpout{ll,:});
						obj.diagData.(vars{i})(j).data(:,:,ll,:)  = single(obj.diagData.(vars{i})(j).data(:,:,ll,:));
					end
					obj.diagData.(vars{i})(j).data  = squeeze(obj.diagData.(vars{i})(j).data);
					obj.diagData.(vars{i})(j).depth = depth;
					obj.diagData.(vars{i})(j).units = curVar.units{j};
					obj.diagData.(vars{i})(j).name  = curVar.name{j};
				end % end j-loop
			end % end i-loop
		end % end method loadDiag

		%--------------------------------------------------------------------------------
			function obj = make_zavg(obj,avgfile)
			% ------------------
			% - Vertically interpolates the sigma coordinates ROMS file to 
			% - WOA constant depth levels
			%
			% - Usage:
			% - obj = make_zavg(obj,avgfile)
			%
			% - Inputs:
			% - avgfile - file to perform depth slices on 
			%
			% - Example:
			% - obj = make_zavg(obj,'avg_2009.nc');
			% ------------------	
			
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
				cmd = ['zslice ',ss{:},' ',[gname,gext],' ',avgfile]; 
				system(cmd);
				cmd = ['nczip -1 -vs z_',avgfile];
				system(cmd);
				cmd = ['mv z_',avgfile,' ',path2];
				system(cmd)
			end

			% - return to original directory
			cd(od);

			end % end method make_zavg

		%--------------------------------------------------------------------------------
		function [fig,cb] = mapPlot(obj,dat,varargin);
			% -----------------------
			% - A way to quickly plot a 2D field
			%
			% - Usage:
			%	- [fig,cb] = mapPlot(obj,dat,varargin);
			%
			% - Inputs:
			%	- dat = 2D field to plot (prepare it before using the script)
			%
			% - Varargin:
			%	- lonbounds:  x-boundaries (defaults to whole domain)
			%	- latbounds:  y-boundaries (defaults to whole domain)
			%	- ticks:      2 = fancy, 1 = on, 0 = off
			%	- background: background color (default 'LightGray');
			%	- coastcolor: coast color (default 'DimGray');
			%	- fontsize:   default 10
			%	- figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			%	- figdim:     default 1 (same height as width, its a multiplier)
			%	- prc:        percentage to limit colorbar axes
			%	- bal:        balance colorbar around 0 (1 == yes, 0 == no)
			%	- levels:     hard-coded levels to plot (can't be used with A.prc)
			%	- cmap:       colormap(default = thermal for no balance, balance for balance)
			% -----------------------

			% User-inputs
			A.lonbounds  = [];
			A.latbounds  = [];
			A.ticks      = 0;
			A.background = rgb('LightGray');
			A.coastcolor = rgb('DimGray');
			A.fontsize   = 10;
			A.figtype    = 'mfig';
			A.figdim     = 1;
			A.prc        = 2;
			A.bal        = 0;
			A.levels     = [];
			A.cmap       = cmocean('thermal');
			A = parse_pv_pairs(A,varargin);

			% Balance override
			if A.bal == 1
				A.cmap = cmocean('balance');
			end

			% Make double
			lon = double(obj.region.lon_rho);
			lat = double(obj.region.lat_rho);

			% Get auto-bounds if empty
			if isempty(A.lonbounds) & isempty(A.latbounds);
				lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
				latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
			elseif isempty(A.lonbounds);
				lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
				latbounds = A.latbounds;
			elseif isempty(A.latbounds);
				lonbounds = A.lonbounds;
				latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
			else
				lonbounds = A.lonbounds;
				latbounds = A.latbounds;
			end

			% Set up ticks
			if abs(diff(lonbounds)) > 50
				dx = 10;
			else
				dx = 5;
			end
			if abs(diff(latbounds)) > 50
				dy = 10;
			else
				dy = 5;
			end
			latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
			lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));

			% Initiate figure
			fig = piofigs(A.figtype,A.figdim);

			% Get colormap limits
			if isempty(A.levels)
				clims = romsMaster.prclims(dat,'prc',A.prc,'bal',A.bal);
				clevs = clims(1):(diff(clims)/31):clims(2);
			else
				clevs = A.levels;
				clims = [A.levels(1) A.levels(end)];
			end
			dat(dat<clevs(1))   = clevs(1);
			dat(dat>clevs(end)) = clevs(end);

			% Make map
			set(0,'CurrentFigure',fig);
			m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
			hold on
			if A.ticks == 0
				m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
					   'xticklabels',[],'yticklabels',[],'backgroundcolor',rgb('LightGray')); drawnow
			elseif A.ticks == 1
				m_grid('box','on','linestyle','none','xtick',lonticks,...
					   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
			elseif A.ticks == 2
				m_grid('box','fancy','linestyle','none','xtick',lonticks,...
					   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
			end
			hold on
			m_contourf(lon,lat,dat,clevs,'LineStyle','none');
			cb = colorbar; drawnow
			cb.FontSize = A.fontsize;
			caxis([clims])
			ax = get(gca);
			m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
			colormap(gca,A.cmap);
			
			% Print figure
			if nargout<1
				pltjpg(1);
			end
		end % end method mapPlot

		%--------------------------------------------------------------------------------
		function [fig] = quickMap(obj,varargin);
			% -----------------------
			% - A way to quickly plot a ROMS map
			%
			% - Usage:
			%	- [fig] = quickMap(obj,varargin);
			%
			% - Varargin:
			%	- lonbounds:  x-boundaries (defaults to whole domain)
			%   - latbounds:  y-boundaries (defaults to whole domain)
			%   - ticks:      2 = fancy, 1 = on, 0 = off
			%   - background: background color (default 'LightGray');
			%   - coastcolor: coast color (default 'DimGray');
			%   - fontsize:   default 10
			%   - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			%   - figdim:     default 1 (same height as width, its a multiplier)
			% -----------------------

			% Reject if no output provided
			if nargout < 1
				disp('No point using this without [fig] output');
				return
			end

			% User-inputs
			A.lonbounds  = [];
			A.latbounds  = [];
			A.ticks      = 0;
			A.background = rgb('LightGray');
			A.coastcolor = rgb('DimGray');
			A.fontsize   = 10;
			A.figtype    = 'mfig';
			A.figdim     = 1;
			A = parse_pv_pairs(A,varargin);

			% Make double
			lon = double(obj.region.lon_rho);
			lat = double(obj.region.lat_rho);

			% Get auto-bounds if empty
			if isempty(A.lonbounds) & isempty(A.latbounds);
				lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
				latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
			elseif isempty(A.lonbounds);
				lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
				latbounds = A.latbounds;
			elseif isempty(A.latbounds);
				lonbounds = A.lonbounds;
				latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
			end

			% Set up ticks
			if abs(diff(lonbounds)) > 50
				dx = 10;
			else
				dx = 5;
			end
			if abs(diff(latbounds)) > 50
				dy = 10;
			else
				dy = 5;
			end
			latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
			lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));

			% Make map
			fig = piofigs(A.figtype,A.figdim);
			set(0,'CurrentFigure',fig);
			m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
			hold on
			if A.ticks == 0
				m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
					   'xticklabels',[],'yticklabels',[],'backgroundcolor',rgb('LightGray')); drawnow
			elseif A.ticks == 1
				m_grid('box','on','linestyle','none','xtick',lonticks,...
					   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
			elseif A.ticks == 2
				m_grid('box','fancy','linestyle','none','xtick',lonticks,...
					   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
			end
			m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
		end % end method quickMap

		%--------------------------------------------------------------------------------
		function [figs,cbs] = mapCmp(obj,dat1,dat2,varargin);
			% -----------------------
			% - A way to quickly plot a 2D field comparison
			% - Produces 3 plots: dat1 and dat2 fields with the same colorbar, 
			%   and a 3rd plot of the difference between them (dat1 - dat2)
			%
			% - Usage:
			%	- [figs,cbs] = mapCmp(obj,dat1,dat2,varargin);
			%
			% - Inputs:
			%	- dat1 = 2D field to plot (prepare it before using the script)
			%	- dat2 = same same but different 2D field to plot (prepare it before using the script)
			%
			% - Varargin:
			%	- lonbounds:  x-boundaries (defaults to whole domain)
			%   - latbounds:  y-boundaries (defaults to whole domain)
			%   - ticks:      2 = fancy, 1 = on, 0 = off
			%   - background: background color (default 'LightGray');
			%   - coastcolor: coast color (default 'DimGray');
			%   - fontsize:   default 10
			%   - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			%   - figdim:     default 1 (same height as width, its a multiplier)
			%	- prc:        percentage to limit colorbar axes
			%	- bal:        balance colorbar around 0 (1 == yes, 0 == no, 2 == set min to 0)
			%	- levels:     hard-coded levels to plot (overrides A.prc, A.bal)
			%   - difflevels: hard-coded difference levels to plot
			%	- cmap:       colormap(default = thermal for no balance, balance for balance)
			% -----------------------

			% User-inputs
			A.lonbounds  = [];
			A.latbounds  = [];
			A.ticks      = 0;
			A.background = rgb('LightGray');
			A.coastcolor = rgb('DimGray');
			A.fontsize   = 10;
			A.figtype    = 'mfig';
			A.figdim     = 1;
			A.prc        = 0.5;
			A.bal        = 0;
			A.levels     = [];
			A.difflevels = [];
			A.cmap       = 'thermal';
			A = parse_pv_pairs(A,varargin);

			% Balance override
			if A.bal == 1
				A.cmap = cmocean('balance');
			end		
			
			% Get universal levels
			if isempty(A.levels)
				all_dat = [dat1(:) dat2(:)];
				A.levels = romsMaster.prclims(all_dat,'prc',A.prc,'bal',A.bal);
				A.levels = linspace(A.levels(1),A.levels(2),20);
			end

			% Make figs(1) and figs(2)
			[figs(1),cbs(1)] = mapPlot(obj,dat1,...
				'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
				'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',cmocean(A.cmap,length(A.levels)-1));
			[figs(2),cbs(2)] = mapPlot(obj,dat2,...
				'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
				'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',cmocean(A.cmap,length(A.levels)-1));

			% Get differences
			diff_dat  = dat1 - dat2;
			if isempty(A.difflevels)
				A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1); 
				A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
			end
			
			% Make figs(3)
			[figs(3),cbs(3)] = mapPlot(obj,diff_dat,...
				'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
				'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels,'cmap',cmocean('balance',length(A.difflevels)-1));

			% Auto-print if no output
			if nargout == 0
				set(0,'CurrentFigure',figs(1));
				pltjpg(1);
				
				set(0,'CurrentFigure',figs(2));
				pltjpg(2);

				set(0,'CurrentFigure',figs(3));
				pltjpg(3);
			end
		end % end method mapCmp

		%--------------------------------------------------------------------------------
		function [fig,cb] = slicePlot(obj,dat,varargin)
			% ------------------
			% Loads and interpolates monthly 2D validation data to the ROMS grid
			% See initDiags for available fields or type (obj.paths.diag)
			%
			% Usage:
			% - [fig,cb] = slicePlot(obj,dat,varargin);
			%
			% Inputs:
			% - dat: output structure produced from sliceROMS or sliceDiag
			%
			% Options:
			% - slice: if several slices were made, specify which one you want to plot here
			% - xlims: hard-coded x-limits (degrees)
			% - zlims: hard-coded z-limits (meters)
			% - cmap:  colormaps (if used, must be a cell array of length (vars))
			% - fontsize:   default 10
			% - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			% - figdim:     default 1 (same height as width, its a multiplier)
			% - prc:        percentage to limit colorbar axes
			% - bal:        balance colorbar around 0 (1 == yes, 0 == no)
			% - levels:     hard-coded levels to plot (can't be used with A.prc)
			%
			% Example:
			% - [fig,cb] = slicePlot(obj,dat,'xlims',[140 260],'zlims',[-500 0]);
			% ------------------

			% User-inputs
			A.slice      = [];
			A.xlims      = [];
			A.zlims      = [];
			A.cmap       = [];
			A.fontsize   = 10;
			A.figtype    = 'mfig';
			A.figdim     = 0.33;
			A.prc        = 0.5;
			A.bal        = 0;
			A.levels     = [];
			A = parse_pv_pairs(A,varargin);

			% Check that A.slice has also been reduced
			if ndims(obj.slice.depth)==3 & isempty(A.slice) 
				disp('You forgot to specify which slice! Killing');
				return
			elseif isempty(A.slice);
				A.slice = 1;
			end
			tmpdeg = obj.slice.deg(:,:,A.slice);
			if ndims(obj.slice.depth)==3
				tmpdep = obj.slice.depth(:,:,A.slice);
			else
				tmpdep = obj.slice.depth;
			end
			
			% Reduce data?
			if ~isempty(A.xlims);
				xind = A.xlims(1) <= tmpdeg & tmpdeg <= A.xlims(2);
				dat(xind==0) = NaN;
			end
			if ~isempty(A.zlims);
				zind = A.zlims(1) <= tmpdep & tmpdep <= A.zlims(2);
				dat(zind==0) = NaN;
			end

			% Get universal levels
			if isempty(A.levels)
				A.levels = romsMaster.prclims(dat(:),'prc',A.prc,'bal',A.bal);
				A.levels = linspace(A.levels(1),A.levels(2),20);
			end
			dat(dat<A.levels(1)) = A.levels(1);
			dat(dat>A.levels(end)) = A.levels(end);
			
			% Generate figure	
			fig = piofigs(A.figtype,A.figdim);
			if strcmp(obj.slice.coord,'latitude') | strcmp(obj.slice.coord,'longitude');
				contourf(tmpdeg,tmpdep,dat,A.levels,'linestyle','none');
				set(gca,'YDir','Reverse');
			end
			if ~isempty(A.xlims);
				xlim(A.xlims);
			end
			if ~isempty(A.zlims);
				ylim(A.zlims);
			end
			if strcmp(obj.slice.coord,'latitude')
				xlabel('Longitude','Interpreter','Latex','FontSize',A.fontsize);
			elseif strcmp(obj.slice.coord,'longitude')	
				xlabel('Latitude','Interpreter','Latex','FontSize',A.fontsize);
			end
			ylabel('Depth (m)','Interpreter','Latex','FontSize',A.fontsize);
			cb = colorbar('location','eastoutside');;
			caxis([A.levels(1) A.levels(end)]);
			if ~isempty(A.cmap)
				set(gca,'Colormap',cmocean(A.cmap,length(A.levels)-1));
			end
			set(gca,'FontSize',A.fontsize);

			% Print if no output provided
			if nargout<1
				pltjpg(1);
			end
		end % end method slicePlot

		%--------------------------------------------------------------------------------
		function [figs,cbs] = sliceCmp(obj,dat1,dat2,varargin)
			% ------------------
			% A way to quickly plot slice comparisons
			% Produces 3 plots: dat1 and dat2 fields with the same colorbar, 
			% and a 3rd plot of the difference between them (dat1 - dat2)
			%
			% Usage:
			% - [figs,cbs] = sliceCmp(obj,dat1,dat2,varargin);
			%
			% Inputs:
			% - dat1 = 2D field to plot (prepare it before using the script)
			% - dat2 = same same but different 2D field to plot (prepare it before using the script)
			%
			% Options:
			% - slice: if several slices were made, specify which one you want to plot here
			% - xlims: hard-coded x-limits (degrees)
			% - zlims: hard-coded z-limits (meters)
			% - cmap:  colormaps (if used, must be a cell array of length (vars))
			% - fontsize:   default 10
			% - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			% - figdim:     default 1 (same height as width, its a multiplier)
			% - prc:        percentage to limit colorbar axes
			% - bal:        balance colorbar around 0 (1 == yes, 0 == no)
			% - levels:     hard-coded levels to plot (can't be used with A.prc)
			%
			% Example:
			% - [figs,cbs] = sliceCmp(obj,dat1,dat2,'xlims',[140 260],'zlims',[-500 0]);
			% ------------------

			% User-inputs
			A.slice      = [];
			A.xlims      = [];
			A.zlims      = [];
			A.cmap       = [];
			A.fontsize   = 10;
			A.figtype    = 'mfig';
			A.figdim     = 0.33;
			A.prc        = 2;
			A.bal        = 0;
			A.levels     = [];
			A.difflevels = [];
			A = parse_pv_pairs(A,varargin);

			% Check that A.slice has also been reduced
			if ndims(obj.slice.depth)==3 & isempty(A.slice)
				disp('You forgot to specify which slice! Killing');
				return
			elseif isempty(A.slice);
				A.slice = 1;
			end
			tmpdeg = obj.slice.deg(:,:,A.slice);
			if ndims(obj.slice.depth)==3
				tmpdep = obj.slice.depth(:,:,A.slice);
			elseif ndims(obj.slice.depth)==2
				tmpdep = obj.slice.depth;
			end

			% Reduce data
			if ~isempty(A.xlims);
				xind = A.xlims(1) <= tmpdeg & tmpdeg <= A.xlims(2);
				dat1(xind==0) = NaN;
				dat2(xind==0) = NaN;
			end
			if ~isempty(A.zlims);
				zind = A.zlims(1) <= tmpdep & tmpdep <= A.zlims(2);
				dat1(zind==0) = NaN;
				dat2(zind==0) = NaN;
			end
			if isempty(A.levels);
				alldat  = [dat1(:) dat2(:)];
				lvls     = romsMaster.prclims(alldat,'prc',A.prc,'bal',A.bal);
				A.levels = linspace(lvls(1),lvls(2),20);
			end

			% Make figs(1) and figs(2)
			[figs(1),cbs(1)] = slicePlot(obj,dat1,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);
			[figs(2),cbs(2)] = slicePlot(obj,dat2,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);

			% Get differences
			diff_dat  = dat1 - dat2;
			if isempty(A.difflevels)
				A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1);
				A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
			end

			% Make figs(3), difference
			[figs(3),cbs(3)] = slicePlot(obj,diff_dat,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap','balance','fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels);

		end % end method sliceCmp

		%--------------------------------------------------------------------------------
		function [obj] = OMZthick(obj,omzthresh)
			% ------------------
			% Calculates OMZ thickness for ROMS and diagnostic oxygen products:
			% WOA-18 and Bianchi (2012) objmap2
			%
			% Usage:
			% - [obj] = OMZthick(obj,thresh)  
			%
			% Inputs:
			% - omzthresh: mmol/m3 thresholds in O2 to define 'oxygen minimum zone'
			%
			% Example:
			% - [obj] = OMZthick(obj,10);
			%
			% This will calculate OMZ (defined as O2 < 10 uMol) thickness from ROMS and
			% validation productions (WOA-18, Bianchi objective mapping 2)
			% ------------------

			% Load first validation set
			% WOA-18
			fname = ['/data/project3/data/woa18/oxygen/1p0/RFcorrected_woa18_o2_clim.nc'];
			woavars = {'lat','lon','depth','Hz','o2'};
			for i = 1:length(woavars)
				woa.(woavars{i}) = ncread(fname,woavars{i});
			end
			clear  woavars fname
			woa.lon(woa.lon<0) = woa.lon(woa.lon<0) + 360; % fix lon (0 - 360);
			woa.o2(woa.o2<0) = 0;
			woa.o2  = nanmean(woa.o2,4);
	
			% Load second validation set    
			% Bianchi 2012 objmap2
			fname = ['/data/project1/demccoy/data/Bianchi2012/o2_datasets.mat'];
			load(fname);
			om2.o2    = o2_datasets.o2_gamma_barnes10;
			om2.depth = o2_datasets.depth; 
			om2.lon   = o2_datasets.lon;
			om2.lon(om2.lon<0) = om2.lon(om2.lon<0)+360; % fix lon (0 - 360):
			om2.lat   = o2_datasets.lat;
			tmpdz     = mean([om2.depth(1:end-1) om2.depth(2:end)],2);
			om2.depth_bnds(:,1) = [0;tmpdz(1:end)];
			om2.depth_bnds(:,2) = [tmpdz(1:end);om2.depth(end)];
					
			% Calculate OMZ thickness from first validation set
			woa.dz  = permute(repmat(woa.Hz,[1 360 180]),[2 3 1]);
			for i = 1:length(omzthresh)
				omzind  = find(woa.o2 < omzthresh(i));
				tmpfill = zeros(size(woa.o2));
				tmpfill(omzind) = 1;
				tmpfill = woa.dz .* tmpfill; % blanking dz where O2>thresh
				omz.woa(:,:,i) = nansum(tmpfill,3);
			end

			% Calculate OMZ thickness from second set
			om2.dz  = om2.depth_bnds(:,2) - om2.depth_bnds(:,1);
			om2.dz  = permute(repmat(om2.dz,[1 360 180]),[2 3 1]);
			for i = 1:length(omzthresh)
				omzind  = find(om2.o2 < omzthresh(i));
				tmpfill = zeros(size(om2.o2));
				tmpfill(omzind) = 1;
				tmpfill = om2.dz .* tmpfill; % blanking dz where O2>thresh
				omz.om2(:,:,i) = nansum(tmpfill,3);
			end

			% Get meshgrids
			[woa.LAT,woa.LON]   = meshgrid(woa.lat,woa.lon);
			[om2.LAT,om2.LON]   = meshgrid(om2.lat,om2.lon);
			
			% Blank data outside ROMS grid
			[in,on] = inpolygon(woa.LON,woa.LAT,obj.region.polygon(:,1),obj.region.polygon(:,2));
			woaind = find(in == 1 | on == 1);
			[in,on] = inpolygon(om2.LON,om2.LAT,obj.region.polygon(:,1),obj.region.polygon(:,2));
			om2ind = find(in == 1 | on == 1);
			for i = 1:length(omzthresh)
				tmpdat = squeeze(omz.woa(:,:,i));
				tmpfill = nan(size(tmpdat));
				tmpfill(woaind) = tmpdat(woaind);
				omz.woa(:,:,i) = tmpfill;
				tmpdat = squeeze(omz.om2(:,:,i));
				tmpfill = nan(size(tmpdat));
				tmpfill(om2ind) = tmpdat(om2ind);
				omz.om2(:,:,i) = tmpfill;
			end

			% Now convert obs grids to ROMS grid
			for i = 1:length(omzthresh);
				tmpwoa  = romsMaster.grid_to_grid(woa.LON,woa.LAT,omz.woa(:,:,i),obj.region.lon_rho,obj.region.lat_rho);
				tmp.woa(:,:,i)  = tmpwoa .* obj.region.mask_rho;
				tmpom2  = romsMaster.grid_to_grid(om2.LON,om2.LAT,omz.om2(:,:,i),obj.region.lon_rho,obj.region.lat_rho);
				tmp.om2(:,:,i)  = tmpom2 .* obj.region.mask_rho;
			end
			omz.woa  = tmp.woa;
			omz.om2  = tmp.om2;

			% Load raw O2 data, get thickness of OMZ
			omz.roms = nan(obj.region.nx,obj.region.ny,length(omzthresh));
			if isfield(obj.romsData,'O2');
				clearO2 = 0;
			else
				clearO2 = 1;
			end
			obj = loadData(obj,{'O2'},'type','raw');
			% Find thickness of OMZ over all months
			for i = 1:obj.region.nt
				tmpHz   = squeeze(obj.region.Hz(:,:,:,i));
				tmpO2   = squeeze(obj.romsData.O2.data(:,:,:,i));
				for j = 1:length(omzthresh)
					omzind  = find(tmpO2 <= omzthresh(j));
					tmpfill = zeros(size(tmpHz)); tmpfill(omzind) = 1;
					tmpfill = tmpHz .* tmpfill; % blanks dz where O2>omzthresh
					% Fill thickness for time j
					tmp.roms(:,:,i,j) = nansum(tmpfill,3);
				end
			end
			% Get annual average for current_file
			omz.roms = squeeze(nanmean(tmp.roms,3)) .* obj.region.mask_rho;
			
			% Save results
			% ROMS
			obj.romsData.OMZ.data   = omz.roms;
			obj.romsData.OMZ.name   = 'OMZ Thickness';
			obj.romsData.OMZ.units  = '$m$';
			obj.romsData.OMZ.thresh = omzthresh;
			% DIAG
			obj.diagData.OMZ(1).data   = omz.woa;
			obj.diagData.OMZ(1).name   = 'OMZ Thickness (WOA-18)';
			obj.diagData.OMZ(1).units  = '$m$';
			obj.diagData.OMZ(1).thresh = omzthresh;
			obj.diagData.OMZ(2).data   = omz.om2;
			obj.diagData.OMZ(2).name   = 'OMZ Thickness (Bianchi 2012)';
			obj.diagData.OMZ(2).units  = '$m$';
			obj.diagData.OMZ(2).thresh = omzthresh;
			
			% Clear data
			if clearO2 == 1; obj.romsData.O2 = []; end
		end % end method OMZthick

	end % end methods declarations
	%----------------------------------------------------------------------------------------

	%----------------------------------------------------------------------------------------
	% Utility functions (static methods)
	methods (Static)
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

		%--------------------------------------------------------------------------------
		function dis = gc_dist(lon1,lat1,lon2,lat2); 
			% -------------------
			% - Distance between 2 points along a great circle
			% - Lat/lon in Radians
			% - Jeroen Molemaker, UCLA 2008
			% -------------------
			dlat = lat2-lat1;
			dlon = lon2-lon1;
			dang = 2*asin( sqrt( sin(dlat/2).^2 + cos(lat2).*cos(lat1).*sin(dlon/2).^2 ) );  %% haversine function
			r_earth = 6371315.;
			dis = r_earth*dang;
		end % end static method gc_dist

		%--------------------------------------------------------------------------------
		function [tmpslice,tmpmask,tmpdepth] = lonlat_slice(dims,lonlat,data,mask,deg,depth,fillmat)
			% -------------------
			% - Function to perform slices along lat/lon, but with parfor
			% - ...a separate function should reduce 'broadcast' variables
			%
			% - Usage:
			% [tmpslice,tmpmask,tmpdepth] = lonlat_slice(dims,lonlat,data,mask,deg,depth)
			% -------------------

			% Process inputs		
			dmsn = dims(1);
			nz   = dims(2);
			nl   = dims(3);
			ns   = dims(4);
			nx   = dims(5);
			ny   = dims(6);
			nv   = dims(7);

			% Set up arrays
			tmpslice = NaN(dmsn,nz,nl,ns,nv);
			tmpdepth = NaN(dmsn,nz,nl,ns);
			tmpmask  = NaN(dmsn,nz,nl,ns);
					   
			% - Start parpool
			delete(gcp('nocreate'));
			parpool(12);
			parfor rcrd = 1:nl
				% - Only get data for that month
				tmpdata = squeeze(data(:,:,:,rcrd,:));
				% - Interp to lon/lat line
				for i = 1:dmsn
					for j = 1:ns
						% - Fill tmpslice
						for z = 1:nz
							if dmsn == nx
								for v = 1:nv
									tmpslice(i,z,rcrd,j,v) = interp1(squeeze(lonlat(i,:)),squeeze(tmpdata(i,:,z,v)),deg(j));
								end
								tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(i,:)),squeeze(mask(i,:,z)),deg(j));
							elseif dmsn == ny
								for v = 1:nv
									tmpslice(i,z,rcrd,j,v) = interp1(squeeze(lonlat(:,i)),squeeze(tmpdata(:,i,z,v)),deg(j));
								end
								tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(:,i)),squeeze(mask(:,i,z)),deg(j));
							end
						end
					end
				end
				tmpdepth(:,:,rcrd,:) = depth'.*ones(fillmat);
			end
			delete(gcp);
		end % end static method lonlat_slice

		%--------------------------------------------------------------------------------
		function [salt_abs] = getSA(salt,pres,lon,lat)
			% -------------------
			% Function to call gsw_SA_from_SP, but in parallel
			% 
			% Usage:
			% - [salt_abs] = getSA(salt,pres,lon,lat);
			%
			% Inputs:
			% - salt, pres, lon, lat = ROMS variables...must be the same size!
			% 
			% -------------------

			% Start parpool
			delete(gcp('nocreate'));
			parpool(12);

			% Create filler
			salt_abs = NaN(size(salt));
			dims     = size(salt);
			nt       = dims(4);
			nz       = dims(3);
			tic
			parfor i = 1:nt
				for z = 1:nz;
					tmp_salt = squeeze(salt(:,:,z,i));
					tmp_pres = squeeze(pres(:,:,z,i));
					tmp_lon  = squeeze(lon(:,:,z,i));
					tmp_lat  = squeeze(lat(:,:,z,i));
					tmp_SA   = gsw_SA_from_SP(tmp_salt(:),tmp_pres(:),tmp_lon(:),tmp_lat(:));
					salt_abs(:,:,z,i) = reshape(tmp_SA,dims(1),dims(2));
				end
			end
			toc
			delete(gcp('nocreate'));

		end % end static method getSA

		%--------------------------------------------------------------------------------
		function lims = prclims(data,varargin)
			% -------------------
			% Script to automatically create limits of data to show
			% i.e. xlim( ) or BinLimits( ) 
			%
			% Usage:
			%   lims = prclims(data,varargin)
			%
			% Inputs:
			%   prc -- percentile to limit data by (default 2% to 98%);
			%   bal -- balance limits (i.e. for histograms around 0)
			%          0 = auto-limits, 1 = evenly balanced around 0, 2 = set min to 0
			% -------------------

			% Default arguments
			A.prc = [2]; % default range 2 - 98%
			A.bal = [1]; 
			A     = parse_pv_pairs(A,varargin);

			% Get limits
			low = prctile(data(:),[A.prc]);
			hih = prctile(data(:),[100-A.prc]);

			% If low < 0 and hih > 0, check for balance
			if low < 0 & hih > 0 & A.bal ~= 2
				lims = [-(max(abs([low hih]))) max(abs([low hih]))];
			elseif A.bal == 1
				lims = [-(max(abs([low hih]))) max(abs([low hih]))];
			elseif A.bal == 2
				lims = [0 hih];
			else
				lims = [low hih];
			end
		end % end static method prclims

		%--------------------------------------------------------------------------------
		function [outvar] = grid_to_grid(inlon,inlat,invar,outlon,outlat)
			% --------------------------------------------------------------------
			% Function to grid one meshgrid to another.
			% Only works if 'inlon/inlat' encapsulates 'outlon/outlat'
			% i.e. a subgrid inside a larger grid.
			%
			% Usage:
			% - outvar = grid_to_grid(inlon,inlat,invar,outlon,outlat)
			%
			% Inputs:
			% - inlon  = meshgrid of lon points (m x n)
			% - inlat  = meshgrid of lat points (m x n)
			% - invar  = meshgrid of data on surface (m x n x t)
			% - outlon = meshgrid of lat points (mm x nn)
			% - outlat = meshrgid of lat points (mm x nn)
			%
			% Outputs:
			% - outvar = meshgrid of data on surface (mm x nn x t)
			%
			% NOTES:
			% Cannot be used with 4D grids, reduce your grid to 3D in a loop first
			% --------------------------------------------------------------------

			% get indices for reduced domain interpolation
			minlon_rho = min(outlon(:));
			maxlon_rho = max(outlon(:));
			minlat_rho = min(outlat(:));
			maxlat_rho = max(outlat(:));
			outer      = 3;
			idx = find(inlon(:) > minlon_rho-outer ...
				 & inlon(:) < maxlon_rho+outer ...
				 & inlat(:) > minlat_rho-outer ...
				 & inlat(:) < maxlat_rho+outer);
				
			% interpolate over all time steps
			[x]  = size(invar);
			[xi] = length(x);
			if xi == 2 % no time
				tmpvar = invar;
				F      = scatteredInterpolant(double(inlon(idx)),double(inlat(idx)),...
											  double(tmpvar(idx)),'linear','nearest');
				outvar = F(double(outlon),double(outlat));
			else 
				for t = 1:x(end)
					fprintf([num2str(t),'...']);
					if xi == 3
						tmpvar = squeeze(invar(:,:,t));
					elseif xi == 4
						tmpvar = squeeze(invar(:,:,t));
					end
					% - build interpolant and interpolate
					F 		  = scatteredInterpolant(double(inlon(idx)),double(inlat(idx)),...
													 double(tmpvar(idx)),'linear','nearest');
					tmpout{t} = F(double(outlon),double(outlat));
				end
				outvar = cat(3,tmpout{:});
			end
		end % end static method grid_to_grid
	end % end static methods declarations
	%----------------------------------------------------------------------------------------
end % end class
%------------------------------------------------------------------------------------------------
