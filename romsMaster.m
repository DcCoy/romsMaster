%------------------------------------------------------------------------------------------------
classdef romsMaster
%-------------------
% Class used to load ROMS data and make various plots,
% comparison figures against validation data, and many
% other actions.
%
% To begin, simply initialize romsMaster
% - obj = initROMS(romsMaster,sim)
%
% To view available routines
% - methods(obj)
%
% For help
% - help obj.(method)
%------------------
	%----------------------------------------------------------------------------------------
	properties
		info					  % struct containing simulation and variable information (via romsInfo)
		grid					  % struct containing grid data and dimensions (via loadGrid)
		slice					  % struct containing slice coordinates (via sliceROMS/sliceDiag)
		profile					  % struct containing profile coordinates(via getProfile)
		budget					  % struct containing BGC tracer budget results (via getBudg)
		data					  % struct containing ROMS data for comparisons (via loadData, computeVar, sliceROMS)
		diag					  % struct containing validation data for comparions (via loadDiag, sliceDiag)
		paths = getDiagPaths; 	  % struct containing paths to data, directories, diagnostic data, etc (via initROMS)
		params = getParams;       % struct containing N-cycle parameters
		constants = getConstants; % struct containing constants and conversions
	end % end properties
	%----------------------------------------------------------------------------------------
	methods
		%--------------------------------------------------------------------------------
		function obj = initROMS(obj,simName,varargin)
			% -------------------
			% Initialization method: gathers paths and coordinate variables 
			%
			% Usage: 
			% - obj = initROMS(obj,sim)
			% 
			% Inputs:
			% - simName: ROMS simulation (peru_chile_0p1, peru_chile_0p05, or pacmed_0p25 only)
			%
			% Optional Inputs:
			% - runName:  Simulation name (i.e. spinup, VKV4, specific runs)
			% - domain:   [xi_start xi_end eta_start eta_end] for domain 
			% - coast:    [#] to remove data within #km of coast
			%
			% Example:
			% - obj = initROMS(obj,'peru_chile_0p1','runName','VKV4_tune2')         <-- if obj defined
			% - obj = initROMS(romsMaster,'peru_chile_0p1','runName','VKV4_tune2')  <-- if obj undefined
			% -------------------
			
			%  Begin
			disp(' ');
			disp('---------------------------------');
			disp('---------------------------------');
			disp('            romsMaster           ');
			disp('---------------------------------');
			disp('---------------------------------');
			disp(' ');               	
			disp('Initializing simulation paths');

			% Suppress warning messages and addpath to Danny's scripts
			warning off
			addpath /data/project1/demccoy/matlab_scripts/
			addpath /data/project1/demccoy/ROMS/tools/Roms_tools_MF/Preprocessing_tools/

			% Options	
			% Defaults
			default_settings = 0;
			if nargin < 2
				% Load default settings
				default_settings = 1;
				resolution = 10;
				simName    = 'peru_chile_0p1';	
				A.runName  = 'dccoy_VKV4_tune7';
				A.gridName = [simName,'_grd.nc'];
				A.outFreq  = 'annual';
				A.domain   = [];
				A.coast    = [];
				A.depth    = [];
				A.header   = [];
			end
			% Other run defaults
			if ~default_settings
				if strcmp(simName,'peru_chile_0p1')
					resolution = 10;
					A.runName  = 'dccoy_VKV4_tune7';
					A.gridName = [simName,'_grd.nc'];
					A.outFreq  = 'annual';
					A.domain   = []; %[31 341 301 461];
					A.depth    = []; %[-750 inf];
					A.coast    = []; %[20];
					A.header   = [];
				elseif strcmp(simName,'peru_chile_0p05');
					resolution = 5;
					A.runName  = 'dccoy_VKV4_tune7';
					A.gridName = [simName,'_grd.nc']; 
					A.outFreq  = 'daily';
					A.domain   = [];
					A.coast    = [];
					A.depth    = [];
					A.header   = [];
				elseif strcmp(simName,'pacmed_0p25');
					resolution = 25;
					A.runName  = 'dccoy_VKV4_tune5';
					A.gridName = [simName,'_grd_corrected.nc']; 
					A.outFreq  = 'annual';
					A.domain   = [];
					A.coast    = [];
					A.depth    = [];
					A.header   = []; % monthly often output as pacmed_avg_*nc
				elseif strcmp(simName,'USSW1');
					resolution = 1;
					A.runName  = 'dccoy_files';
					A.gridName = ['usw1_grd.nc'];
					A.outFreq  = 'daily';
					A.domain   = [];
                    A.coast    = [];
                    A.depth    = [];
					A.header   = ['ussw1_'];
				end
			end
			% Set default data to all
			A = parse_pv_pairs(A,varargin);

			% Clear any loaded ROMS data
			obj = clearROMS(obj);
			
			% Set file paths
			obj.paths.simPath   = ['/data/project2/model_output/',simName,'/'];
			obj.paths.runPath   = [obj.paths.simPath,simName,'_',A.runName,'/'];
			obj.paths.config    = ['/data/project2/demccoy/ROMS_configs/',simName,'/'];
			obj.paths.grid      = [obj.paths.config,'grid/',A.gridName];
			obj.paths.avg       = [obj.paths.runPath,'avg/',A.outFreq,'/'];
			obj.paths.his       = [obj.paths.runPath,'his/',A.outFreq,'/'];
			obj.paths.phys_flux = [obj.paths.runPath,'phys_flux/',A.outFreq,'/'];

			% initiate directories if they dont exist
			mkdir([obj.paths.runPath,'Figures']);
			mkdir([obj.paths.runPath,'Figures/']);
			mkdir([obj.paths.runPath,'Figures/Diagnostic']);
			mkdir([obj.paths.runPath,'Figures/Budget']);
		
			% grab plot paths
			obj.paths.plots.diag	= [obj.paths.runPath,'Figures/Diagnostic/'];
			obj.paths.plots.budget  = [obj.paths.runPath,'Figures/Budget/'];
			obj.paths.plots.tmpfigs = ['/data/project1/demccoy/tmpfigs/'];
					
			% Initialize info for subregion
			if ~isempty(A.domain)
				obj.grid.region.lon_lim = [A.domain(1:2)];
				obj.grid.region.lat_lim = [A.domain(3:4)];
			else
				obj.grid.region.lon_lim = [];
				obj.grid.region.lat_lim = [];
			end
			if ~isempty(A.coast)
				obj.grid.region.coast_lim = [A.coast];
			else
				obj.grid.region.coast_lim = [];
			end
			if ~isempty(A.depth)
				obj.grid.region.dep_lim = [A.depth];
			else
				obj.grid.region.dep_lim = [-inf inf];
			end

			% Get LaTex friendly runName (runTit)
			string_to_replace  = {'dccoy_','pdamien_','_spinup','_','KV'};
			replacement_string = {'','','','-','Kv'};
			A.runTit  = A.runName;
			for i = 1:length(string_to_replace)
				A.runTit  = strrep(A.runTit,string_to_replace{i},replacement_string{i});
			end
			
			% Initialize info for ROMS variables
			obj.info = A;
			obj.info.simName = simName;
			obj.info.resolution = resolution;
	
			% Get info for files
			obj = romsInfo(obj);
			
			% Load grid
			obj = loadGrid(obj);
		end % end method initROMS

		%--------------------------------------------------------------------------------
		function [obj] = clearROMS(obj) 
			% ----------------------
			% Clears loaded data and coordinates
			%
			% Usage:
			% - [obj] = clearROMS(obj) 
			% ----------------------
			disp('Clearing loaded data');			

			% Clear fields
			obj.slice = [];
			obj.profile = [];
			obj.diag = [];
			obj.data = [];
			obj.budget = [];

			% Clear file-specific grid data
			file_types = {'avg','his','phys_flux'};
			for ff = 1:length(file_types)
				obj.grid.(file_types{ff}) = [];
			end
			obj.grid.nt = [];
		end % end method clearROMS

		%--------------------------------------------------------------------------------
		function [obj] = romsInfo(obj)
			% ----------------------
			% Obtains info from a sample roms file
			%
			% Usage:
			% - [obj] = romsInfo(obj) 
			% ----------------------
			disp('Loading file meta data (ncinfo) into obj.info');

			rmpath /data/project1/demccoy/ROMS/ROMS_tools/nc_tools/

			% List file_types
			file_types  = {'avg','his','phys_flux'};
			file_string = {'avg','his','phys_flux_avg'};
			if ~isempty(obj.info.header)
				for ff = 1:length(file_string)
					file_string{ff} = [obj.info.header,file_string{ff}];
				end
			end
			empty = [];		
			for ff = 1:length(file_types)
				obj.info.(file_types{ff}) = [];
				% Get files and simplify paths
				tmp = dir(obj.paths.(file_types{ff}));
				idx = [];
				for i = 1:length(tmp)
					if strncmp(tmp(i).name,[file_string{ff}],length(file_string{ff}));
						idx = [idx i];
					end
				end
				tmpfiles = {tmp.name};
				obj.info.(file_types{ff}).files  = tmpfiles(idx); 
				clear tmp tmpfiles

				% Get sample file
				file_path = [obj.paths.(file_types{ff})];

				% list variables
				try
					info_file = [file_path,obj.info.(file_types{ff}).files{1}];
					tmp.info = ncinfo(info_file);
				catch
					disp(['No ',file_types{ff},' files found...skipping']);
					empty = [empty ff];
					obj.info.(file_types{ff}) = [];
					continue
				end
				
				tmpfields = fields(obj.info.(file_types{ff}));
				for i = 1:length(tmpfields)
					tmp.info.(tmpfields{i}) = obj.info.(file_types{ff}).(tmpfields{i});
				end
				obj.info.(file_types{ff}) = tmp.info;

				% Get time-steps
				for i = 1:length(obj.info.(file_types{ff}).files);
					tmp.info = ncinfo([file_path,obj.info.(file_types{ff}).files{i}]);
					tmptime(i) = tmp.info.Dimensions(find(strcmp('time',{obj.info.(file_types{ff}).Dimensions.Name})==1)).Length; 
				end
				obj.info.(file_types{ff}).time = tmptime;
				clear tmp;

				% Separate 2D from 3D variables
				cnt2d = 1;
				cnt3d = 1;
				nt = obj.info.(file_types{ff}).Dimensions(find(strcmp('time',{obj.info.(file_types{ff}).Dimensions.Name})==1)).Length;
				for i = 1:length(obj.info.(file_types{ff}).Variables)
					if length(obj.info.(file_types{ff}).Variables(i).Size)==3 & obj.info.(file_types{ff}).Variables(i).Size(3) == nt;
						obj.info.(file_types{ff}).var2d{cnt2d}  = obj.info.(file_types{ff}).Variables(i).Name;
						for j = 1:length(obj.info.(file_types{ff}).Variables(i).Attributes);
							if strcmp(obj.info.(file_types{ff}).Variables(i).Attributes(j).Name,'long_name');
								nameidx = j;
							elseif strcmp(obj.info.(file_types{ff}).Variables(i).Attributes(j).Name,'units');
								unitidx = j;
							end
						end
						obj.info.(file_types{ff}).name2d{cnt2d} = obj.info.(file_types{ff}).Variables(i).Attributes(nameidx).Value;
						obj.info.(file_types{ff}).unit2d{cnt2d} = obj.info.(file_types{ff}).Variables(i).Attributes(unitidx).Value;
						obj.info.(file_types{ff}).idx2d(cnt2d)  = i;
						cnt2d                  = cnt2d + 1;
					elseif length(obj.info.(file_types{ff}).Variables(i).Size)==4
						obj.info.(file_types{ff}).var3d{cnt3d} = obj.info.(file_types{ff}).Variables(i).Name;
						for j = 1:length(obj.info.(file_types{ff}).Variables(i).Attributes);
							if strcmp(obj.info.(file_types{ff}).Variables(i).Attributes(j).Name,'long_name');
								nameidx = j;
							elseif strcmp(obj.info.(file_types{ff}).Variables(i).Attributes(j).Name,'units');
								unitidx = j;
							end
						end
						obj.info.(file_types{ff}).name3d{cnt3d} = obj.info.(file_types{ff}).Variables(i).Attributes(nameidx).Value;
						obj.info.(file_types{ff}).unit3d{cnt3d} = obj.info.(file_types{ff}).Variables(i).Attributes(unitidx).Value;
						obj.info.(file_types{ff}).idx3d(cnt3d)  = i;
						cnt3d                  = cnt3d + 1;
					end
				end
			end
			
			% Update filetypes
			file_types(empty) = [];
			if isempty(file_types)
				disp('ERROR: No ROMS files found, check initROMS inputs');
				kill	
			end
			
			% Replace Unit Strings with Latex version
			strings_to_replace = {'meter second-1','Celsius','PSU','mMol P','mMol O2','mMol N','mMol',...
								  'mMol C','mg Chl-a','mMol CaCO3','meter2 second-1','W m-2','mmol C/m3/s',...
								  'mmol/m3/s','mg Chl/m3','mmol N/m3/s','mmol N2O/m3/s','meter',...
								  'ppm','mmol/m2/s','mmol/m3','m-3'};

			replacement_string = {'$m$ $s^{-1}$','$^{o}C$','$PSU$','$mmol$ $P$','$mmol$ $O_2$','$mmol$ $N$','$mmol$',...
								  '$mmol$ $C$','$mg$ $Chl-A$','$mmol$ $CaCO_3$','$m^{2}$ $s^{-1}$','$W$ $m^{-2}$','$mmol$ $C$ $m^{-3}$ $s^{-1}$',...
								  '$mmol$ $m^{-3}$ $s^{-1}$','$mg$ $C$ $m^{-3}$','$mmol$ $N$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$m$',...
								  '$ppm$','$mmol$ $m^{-2}$ $s^{-1}$','$mmol$ $m^{-3}$','$m^{-3}$'};

			% Finish
			for ff = 1:length(file_types)
				% Grab time info from average file
				fieldnames = {obj.info.(file_types{ff}).Attributes.Name};
				fieldvalue = {obj.info.(file_types{ff}).Attributes.Value};
				dt         = fieldvalue(strcmp(fieldnames,'dt'));
				dt         = double(dt{1});
				navg       = fieldvalue(strcmp(fieldnames,'navg'));
				navg       = double(navg{1});
				ntimes     = fieldvalue(strcmp(fieldnames,'ntimes'));
				ntimes     = double(ntimes{1});

				% Compute dt according to output frequency
				obj.info.(file_types{ff}).Freq   = dt*navg/86400;
				obj.info.(file_types{ff}).Ntimes = ((dt*ntimes)/86400)/(obj.info.(file_types{ff}).Freq);
				obj.info.(file_types{ff}).dt     = dt*navg;
				if obj.info.(file_types{ff}).Freq > 27 & obj.info.(file_types{ff}).Freq < 32
					obj.info.(file_types{ff}).time_string     = ['Monthly'];
					obj.info.(file_types{ff}).time_string_idv = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
				elseif obj.info.(file_types{ff}).Freq == 1
					obj.info.(file_types{ff}).time_string = ['Daily'];
					obj.info.(file_types{ff}).time_string_idv = num2cell(1:31);
				elseif obj.info.(file_types{ff}).Freq < 1
					obj.info.(file_types{ff}).time_string = ['Hourly'];
					obj.info.(file_types{ff}).time_string_idv = num2cell(1:24);
				end

				% Replace strings
				for i = 1:length(strings_to_replace)
					obj.info.(file_types{ff}).unit3d = strrep(obj.info.(file_types{ff}).unit3d,strings_to_replace{i},replacement_string{i});
					obj.info.(file_types{ff}).unit2d = strrep(obj.info.(file_types{ff}).unit2d,strings_to_replace{i},replacement_string{i});
				end
			
				% Get original grid sizes
				dims = {obj.info.(file_types{ff}).Dimensions.Name};
				for i = 1:length(dims);
					if ~strcmp(dims{i},'time');
						ind = find(strcmp(dims{i},dims)==1);
						obj.info.(file_types{ff}).(dims{i}) = obj.info.(file_types{ff}).Dimensions(ind).Length;
					end
				end

				% Apply any x/y overrides
				if isempty(obj.grid.region.lon_lim); 
					obj.info.(file_types{ff}).xi_rho = [1 obj.info.(file_types{ff}).xi_rho]; 
					obj.info.(file_types{ff}).xi_u   = [1 obj.info.(file_types{ff}).xi_u];
				else
					obj.info.(file_types{ff}).xi_rho = [obj.grid.region.lon_lim];
					obj.info.(file_types{ff}).xi_u   = [obj.grid.region.lon_lim];
				end
				if isempty(obj.grid.region.lat_lim); 
					obj.info.(file_types{ff}).eta_rho = [1 obj.info.(file_types{ff}).eta_rho]; 
					obj.info.(file_types{ff}).eta_v   = [1 obj.info.(file_types{ff}).eta_v];
				else
					obj.info.(file_types{ff}).eta_rho = [obj.grid.region.lat_lim];
					obj.info.(file_types{ff}).eta_v   = [obj.grid.region.lat_lim];
				end			
					
				% get grid coordinates
				obj.info.(file_types{ff}).ind2D = [obj.info.(file_types{ff}).xi_rho(1) obj.info.(file_types{ff}).eta_rho(1); ...
												   diff(obj.info.(file_types{ff}).xi_rho)+1 diff(obj.info.(file_types{ff}).eta_rho)+1];
				obj.info.(file_types{ff}).ind3D = [obj.info.(file_types{ff}).xi_rho(1) obj.info.(file_types{ff}).eta_rho(1) 1; ...
												   diff(obj.info.(file_types{ff}).xi_rho)+1 diff(obj.info.(file_types{ff}).eta_rho)+1 inf];
				obj.info.(file_types{ff}).ind4D = [obj.info.(file_types{ff}).xi_rho(1) obj.info.(file_types{ff}).eta_rho(1) 1 1; ...
												   diff(obj.info.(file_types{ff}).xi_rho)+1 diff(obj.info.(file_types{ff}).eta_rho)+1 inf inf];
			end
		end % end method romsInfo

		%--------------------------------------------------------------------------------
		function [obj] = loadGrid(obj)
			% ----------------------
			% Loads 2D grid information into obj.grid
			%
			%
			% Usage:
			% - obj = loadGrid(obj)
			% ----------------------
			disp('Loading grid data into obj.grid');

			rmpath('/data/project1/demccoy/ROMS/ROMS_tools/nc_tools/');
			addpath('/usr/local/MATLAB/R2019b/toolbox/matlab/imagesci/');

			% Requires romsInfo
			try; obj.info.avg.ind2D;
			catch; obj = romsInfo(obj);
			end

			% Load grid
			u2D = obj.info.avg.ind2D;
			v2D = obj.info.avg.ind2D;
			u2D(2,1) = u2D(2,1)-1;
			v2D(2,2) = v2D(2,2)-1;
			obj.grid.lon_rho  = double(ncread(obj.paths.grid,'lon_rho',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			try
				obj.grid.lon_u    = double(ncread(obj.paths.grid,'lon_u',[u2D(1,:)],[u2D(2,:)]));
				obj.grid.lat_u    = double(ncread(obj.paths.grid,'lat_u',[u2D(1,:)],[u2D(2,:)]));
			catch
			end
			obj.grid.lat_rho  = double(ncread(obj.paths.grid,'lat_rho',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			try
				obj.grid.lon_v    = double(ncread(obj.paths.grid,'lon_v',[v2D(1,:)],[v2D(2,:)]));
				obj.grid.lat_v    = double(ncread(obj.paths.grid,'lat_v',[v2D(1,:)],[v2D(2,:)]));
			catch
			end
			obj.grid.pm       = double(ncread(obj.paths.grid,'pm',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			obj.grid.pn       = double(ncread(obj.paths.grid,'pn',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			obj.grid.angle    = double(ncread(obj.paths.grid,'angle',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			obj.grid.mask_rho = double(ncread(obj.paths.grid,'mask_rho',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));
			try
				obj.grid.mask_u   = double(ncread(obj.paths.grid,'mask_u',[u2D(1,:)],[u2D(2,:)]));
				obj.grid.mask_v   = double(ncread(obj.paths.grid,'mask_v',[v2D(1,:)],[v2D(2,:)]));
			catch
			end
			obj.grid.h        = double(ncread(obj.paths.grid,'h',[obj.info.avg.ind2D(1,:)],[obj.info.avg.ind2D(2,:)]));

			% lon360 fix
			obj.grid.lon_rho(obj.grid.lon_rho<0) = obj.grid.lon_rho(obj.grid.lon_rho<0)+360;

			% Get area dimensions
			obj.grid.area     = (1./(obj.grid.pm .* obj.grid.pn));
			obj.grid.dx       = 1./(obj.grid.pm);
			obj.grid.dy       = 1./(obj.grid.pn);
			obj.grid.area_rho = double(obj.grid.dx.*obj.grid.dy);
			
			% Basic dimensions
			obj.grid.nx         = diff(obj.info.avg.xi_rho)+1;
			obj.grid.ny         = diff(obj.info.avg.eta_rho)+1; 
			obj.grid.nz         = obj.info.avg.s_rho; 
			obj.grid.ndim_xy    = [obj.grid.nx,obj.grid.ny];
			obj.grid.ndim_xyz   = [obj.grid.nx,obj.grid.ny,obj.grid.nz];
			obj.grid.minlon_rho = nanmin(obj.grid.lon_rho(:));
			obj.grid.maxlon_rho = nanmax(obj.grid.lon_rho(:));
			obj.grid.minlat_rho = nanmin(obj.grid.lat_rho(:));
			obj.grid.maxlat_rho = nanmax(obj.grid.lat_rho(:));
			clear tmp

			% Fix masks
			obj.grid.mask_rho(obj.grid.mask_rho==0) = NaN;
			try
				obj.grid.mask_u(obj.grid.mask_u==0) = NaN;
				obj.grid.mask_v(obj.grid.mask_v==0) = NaN;
			catch
			end

			% Get grid polygon
			obj.grid.polygon(:,1) = [obj.grid.lon_rho(1,:) obj.grid.lon_rho(:,end)'...
									 fliplr(obj.grid.lon_rho(end,:)) flipud(obj.grid.lon_rho(:,1))'];
			obj.grid.polygon(:,2) = [obj.grid.lat_rho(1,:) obj.grid.lat_rho(:,end)'...
									 fliplr(obj.grid.lat_rho(end,:)) flipud(obj.grid.lat_rho(:,1))'];

			% Also save z_avg_dep
			tmp = load('/data/project1/demccoy/ROMS/validation/wcoord.mat');
			obj.grid.z_avg_dep = tmp.wcoord0p25.depth;

			% Make coord data all single
			obj.grid = romsMaster.struct2double(obj.grid);

			% (OPTIONAL)
			if obj.grid.region.coast_lim>(-inf) 
				% Calculate distance-from-coast (m)
				disp('Refining mask via Dist2Coast');
				obj = Dist2Coast(obj);
				% Update mask
				obj.grid.mask_rho(obj.grid.coastdist < obj.grid.region.coast_lim*1000) = NaN;
			end
		end % end method loadGrid

		%--------------------------------------------------------------------------------
		function [obj] = loadDepth(obj,file,varargin) 
			% ----------------------
			% Load Z-grid info (z_r,z_w), which are output-dependent.
			% Also updates 3D mask for raw data.
			%
			% Usage:
			% - obj = loadDepth(obj,file,varargin) 
			%
			% Inputs:
			% - file = (#s) load specific files 
			%
			% Optional:
			% - type: file type (his, avg (default), phys_flux)
			%
			% Example:
			% - obj = loadDepth(obj,1,'type','avg');
			% ----------------------
			disp('Grabbing z_r, z_w, Hz');

			% Call loadGrid
			try; obj.grid.h;
			catch; obj = loadGrid(obj);
			end
			
			% Process optional inputs
			A.type = 'avg';
			A = parse_pv_pairs(A,varargin);
			files = obj.info.(A.type).files;

			% Load vertical coordinates (and zlevs4 info)
			atts = {obj.info.(A.type).Attributes.Name};
			atts_to_read = {'theta_s','theta_b','hc'};
			for i = 1:length(atts_to_read);
				ind = find(strcmp(atts_to_read{i},atts)==1);
				obj.grid.(atts_to_read{i}) = obj.info.(A.type).Attributes(ind).Value;
			end

			% Initialize matrices to fill
			nx = obj.info.(A.type).xi_rho;  nx = nx(end)-nx(1)+1;
			ny = obj.info.(A.type).eta_rho; ny = ny(end)-ny(1)+1;
			nz = obj.info.(A.type).s_rho;   nz = nz;
			nt = sum(obj.info.(A.type).time(file));
			files = files(file);
			obj.grid.nt = nt;
			obj.grid.(A.type).Hz  = nan(nx,ny,nz,nt);
			obj.grid.(A.type).z_r = nan(nx,ny,nz,nt);
			obj.grid.(A.type).z_u = nan(nx-1,ny,nz,nt);
			obj.grid.(A.type).z_v = nan(nx,ny-1,nz,nt);
			obj.grid.(A.type).z_w = nan(nx,ny,nz+1,nt);
			
			% Go through file(s) and load data (or calculate it via zlevs4)
			for ff = 1:length(files)
				if length(files)>1
					fprintf([num2str(ff),'..']);
				end
				fname = [obj.paths.(A.type),files{ff}];
				% Load SLA
				tmp.zeta = ncread(fname,'zeta',[obj.info.(A.type).ind3D(1,:)],[obj.info.(A.type).ind3D(2,:)]);
				% Average for more accurate depth calc
				tmp.h = obj.grid.h;
				% Initialize
				tmp.z_r = nan(size(tmp.zeta,1),size(tmp.zeta,2),obj.info.(A.type).s_rho,size(tmp.zeta,3));
				tmp.z_w = nan(size(tmp.zeta,1),size(tmp.zeta,2),obj.info.(A.type).s_rho+1,size(tmp.zeta,3));
				% Get depths
				for i = 1:size(tmp.zeta,3);
					z_r = zlevs4(tmp.h,tmp.zeta(:,:,i),obj.grid.theta_s,obj.grid.theta_b,obj.grid.hc,obj.info.(A.type).s_rho,'r','new2012');
					z_w = zlevs4(tmp.h,tmp.zeta(:,:,i),obj.grid.theta_s,obj.grid.theta_b,obj.grid.hc,obj.info.(A.type).s_rho,'w','new2012');
					z_r = permute(z_r,[2 3 1]);
					z_u = (z_r(1:end-1,:,:) + z_r(2:end,:,:))./2;
					z_v = (z_r(:,1:end-1,:) + z_r(:,2:end,:))./2;
					z_w = permute(z_w,[2 3 1]);
					tmp.z_r(:,:,:,i) = z_r;
					tmp.z_u(:,:,:,i) = z_u;
					tmp.z_v(:,:,:,i) = z_v;
					tmp.z_w(:,:,:,i) = z_w;
				end
				tmp.Hz = diff(tmp.z_w,1,3);
				% Save output
				out.Hz{ff} = tmp.Hz;
				out.z_r{ff} = tmp.z_r;
				out.z_u{ff} = tmp.z_u;
				out.z_v{ff} = tmp.z_v;
				out.z_w{ff} = tmp.z_w;
				clear tmp
			end
			if ff>1
				fprintf(['\n']);
			end
			obj.grid.(A.type).Hz  = cat(4,out.Hz{:});
			obj.grid.(A.type).z_r = cat(4,out.z_r{:});
			obj.grid.(A.type).z_u = cat(4,out.z_u{:});
			obj.grid.(A.type).z_v = cat(4,out.z_v{:});
			obj.grid.(A.type).z_w = cat(4,out.z_w{:});

			% Make 3D area_rho, grid_mask with depth limits applied
			obj.grid.(A.type).mask_rho3d = repmat(obj.grid.mask_rho,[1 1 obj.grid.nz nt]);
			try
				obj.grid.(A.type).mask_u3d = repmat(obj.grid.mask_u,[1 1 obj.grid.nz nt]);
				obj.grid.(A.type).mask_v3d = repmat(obj.grid.mask_v,[1 1 obj.grid.nz nt]);
			catch
			end
			obj.grid.(A.type).area3d = repmat(obj.grid.area_rho,[1 1 obj.grid.nz nt]);
			obj.grid.(A.type).volume = obj.grid.(A.type).area3d .* obj.grid.(A.type).Hz;
			if ~isempty(obj.grid.region.dep_lim);			
				obj.grid.(A.type).mask_rho3d(obj.grid.(A.type).z_r<obj.grid.region.dep_lim(1)) = NaN;
				obj.grid.(A.type).mask_rho3d(obj.grid.(A.type).z_r>obj.grid.region.dep_lim(2)) = NaN;
				try
					obj.grid.(A.type).mask_u3d(obj.grid.(A.type).z_u<obj.grid.region.dep_lim(1)) = NaN;
					obj.grid.(A.type).mask_u3d(obj.grid.(A.type).z_u>obj.grid.region.dep_lim(2)) = NaN;
					obj.grid.(A.type).mask_v3d(obj.grid.(A.type).z_v<obj.grid.region.dep_lim(1)) = NaN;
					obj.grid.(A.type).mask_v3d(obj.grid.(A.type).z_v>obj.grid.region.dep_lim(2)) = NaN;
				catch
				end
				obj.grid.(A.type).area3d(obj.grid.(A.type).z_r<obj.grid.region.dep_lim(1)) = NaN;
				obj.grid.(A.type).area3d(obj.grid.(A.type).z_r>obj.grid.region.dep_lim(2)) = NaN;
			end
		end % end method loadDepth

		%--------------------------------------------------------------------------------
		function [fig,ax] = gridView(obj,varargin)
			% ----------------------
			% Plots the lat/lon indices of a grid file
			%
			% Usage:
			% - [fig,ax] = gridView(obj,varargin)
			%
			% Inputs (varargin):
			% - full  = 1 (default), view entire grid. Use 0 for regional grid
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
			A.full  = 1;
			A.dx    = [20];
			A.dy    = [20];
			A.ticks = [1];
			A.font  = [12];
			A.save  = [0];
			A       = parse_pv_pairs(A,varargin);

			% Get grid
			if A.full == 1
				% Get original lon/lat
				tmp.lon_rho = double(ncread(obj.paths.grid,'lon_rho'));
				tmp.lat_rho = double(ncread(obj.paths.grid,'lat_rho'));
			else
				tmp.lon_rho = obj.grid.lon_rho;
				tmp.lat_rho = obj.grid.lat_rho;
			end

			% Plot lon/lat lines
			fig(1) = piofigs('lfig',1.5);
			ax(1)  = map_plot(fig(1),tmp.lon_rho,tmp.lat_rho,'ticks',A.ticks,'font',A.font);	
			[a,b]  = size(tmp.lon_rho);
			for i = 1:A.dx:a
				m_plot(tmp.lon_rho(i,:),tmp.lat_rho(i,:),'r');
				m_text(tmp.lon_rho(i,end),tmp.lat_rho(i,end),num2str(i),'fontsize',8);
				for j = 1:A.dy:b
					hold on
					m_plot(tmp.lon_rho(:,j),tmp.lat_rho(:,j),'b');
					m_text(tmp.lon_rho(end,j),tmp.lat_rho(end,j),num2str(j),'fontsize',8);
				end
			end
			fname = [obj.info.simName,'_grid'];
			if A.save == 1
				export_fig('-jpg',[obj.paths.runPath,fname]);
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

			% Get original lon/lat
			tmp.lon_rho = ncread(obj.paths.grid,'lon_rho');
			tmp.lat_rho = ncread(obj.paths.grid,'lat_rho');
	
			% Map of region
			if ~isempty(obj.grid.lon_rho)
				% Generate whole map
				fig(1)  = piofigs('lfig',1.5);
				set(0,'CurrentFigure',fig(1));
				[ax(1)] = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks,'font',A.font);
				hold on
				% Plot grid box
				m_plot(tmp.lon_rho(1,:),  tmp.lat_rho(1,:),'k','linewidth',2);
				m_plot(tmp.lon_rho(:,1),  tmp.lat_rho(:,1),'k','linewidth',2);
				m_plot(tmp.lon_rho(end,:),tmp.lat_rho(end,:),'k','linewidth',2);
				m_plot(tmp.lon_rho(:,end),tmp.lat_rho(:,end),'k','linewidth',2);
				% Plot region box
				m_plot(obj.grid.lon_rho(1,:),obj.grid.lat_rho(1,:),'--k','linewidth',2);
				m_plot(obj.grid.lon_rho(:,1),obj.grid.lat_rho(:,1),'--k','linewidth',2);
				m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'--k','linewidth',2);
				m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'--k','linewidth',2);
				fname = [obj.info.simName,'_region'];
				if A.save == 1
					export_fig('-jpg',[obj.paths.runPath,fname]);
				else
					pp = 1;
					pltjpg(1);
				end
				if nargout < 1 & pp ~= 1;
					pltjpg(1);
				end
			end
		end % end method regionView

		%--------------------------------------------------------------------------------
		function obj = Dist2Coast(obj)
			% --------------------
			% Call this function to calculate each grid cell's distance-to-coast
			% --------------------
			disp('Calculating distance from coast');

            % Load grid?
            try; obj.grid.lon_rho;
            catch; obj = loadGrid(obj);
            end

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
		function obj = changeInputs(obj,varargin)
			% ----------------------
			% Easier way to change the input files
			% 
			% Usage:
			% - obj = changeInputs(obj,year)
			%
			% Options:
			% - runName = different extension
			%
			% Example:
			% - obj = changeInputs(obj,year,'ext','_PCref');
			% ----------------------

			% Process inputs		
			A.runName = obj.info.runName;
			A.outFreq = obj.info.outFreq;
			A = parse_pv_pairs(A,varargin);

			% Easy path to avoid script
			if strcmp(A.runName,obj.info.runName) & A.outFreq==obj.info.outFreq;
				disp('Skipping changeInputs');
				return
			else
				disp(['Changing inputs']);			
			end

			% Clear loaded data
			obj = clearROMS(obj);

            % Set file paths
            obj.paths.runPath   = [obj.paths.simPath,obj.info.simName,'_',A.runName,'/'];
            obj.paths.avg       = [obj.paths.runPath,'avg/',A.outFreq,'/'];
            obj.paths.his       = [obj.paths.runPath,'his/',A.outFreq,'/'];
            obj.paths.phys_flux = [obj.paths.runPath,'phys_flux/',A.outFreq,'/'];

			% Update info
			obj.info.runName = A.runName;
			obj.info.outFreq = A.outFreq;
			obj = romsInfo(obj);
			
			% Clear file-dependent data
			tmpfields = {'Hz','z_r','z_w','mask_rho3d'};
			for i = 1:length(tmpfields)
				disp('Warning, dont forget to rerun loadDepth');
				obj.grid.(tmpfields{i}) = [];
			end
		end % end method changeInputs
		
		%--------------------------------------------------------------------------------
		function obj = getBudg(obj,vars,file,varargin)
			% --------------------
			% Main method to perform budget analysis on a ROMS tracer (vars).
			%
			% Usage:
			% - obj = getBudg(obj,vars,file,varargin)
			%
			% Inputs:
			% - vars = cell array of budgets to close 
            % - file = (#s) load specific time-averages 
            %        = 0 load all files and average (for monthly)
			%
			% Optional Inputs:
			% - hisfile  = alternate file#s for history files (in case output is limited) 
			% - physfile = alternate file#s for flux files (in case output is limited) 
			% - srcsnk   = (1) to call sources_sinks (default 0)
			% - int      = (1) t o call intBudg (default 0)
			% - clean    = (0) to keep all data loaded (default 1)
			%
			% Example:
			% - obj = getBudg(obj,'N2O',1)
			%
			%   This will instrcut budget methods on which
			%   tracer to perform budget on. Here, tracking N2O.
			%
			% Available tracer budgets:
			% NO3
			% NO2
			% NH4
			% N2O (+AO1,SODEN,SIDEN,ATM)
			% N2
			% --------------------

			% Toggles
			A.int      = [0];
			A.srcsnk   = [0];
			A.clean    = [1];
			A.hisfile  = [];
			A.physfile = [];
			A          = parse_pv_pairs(A,varargin);

			% Check inputs
			if isempty(vars)
				disp('vars must be defined, see romsMaster.getBudg')
				return
			end
			if isempty(A.hisfile)
				hisfile = file;
			else
				hisfile = A.hisfile;
			end
			if isempty(A.physfile)
				physfile = file;
			else
				physfile = A.physfile;
			end

			% Constants (these may change!!! check ROMS code)
			param.DONrefract = 0.0115;   % Fraction of DON to refractory pool
			param.Q          = 0.137;    % N/C ratio

			% Get budget
			for i = 1:length(vars)
				disp(['Computing budget for ',vars{i}]);
				if strcmp(vars{i},'NO3')
					obj.budget.(vars{i}).info.tits   = {'$NO_3$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N$ $m^{-3}$ $s^{-1}$','$mmol$ $N$ $m^{-2}$ $s^{-1}$','$mmol$ $N$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'NITROX','DENITRIF1','SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE'};
					obj.budget.(vars{i}).info.smseq  = [(1) (-1) (-1) (-1) (-1)];
					obj.budget.(vars{i}).info.fluxes = {'SED_DENITRIF'};
					obj.budget.(vars{i}).info.lvls   = {'sed'};
				elseif strcmp(vars{i},'NO2')
					obj.budget.(vars{i}).info.tits   = {'$NO_2$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N$ $m^{-3}$ $s^{-1}$','$mmol$ $N$ $m^{-2}$ $s^{-1}$','$mmol$ $N$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'AMMOX','N2OAMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2',...
														'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
					obj.budget.(vars{i}).info.smseq  = [(1) (-2) (-1) (-1) (1) (-1) (-1) (-1) (-1)];	
					obj.budget.(vars{i}).info.fluxes = {[]};
					obj.budget.(vars{i}).info.lvls   = {[]};
				elseif strcmp(vars{i},'NH4')
					obj.budget.(vars{i}).info.tits   = {'$NH_4$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N$ $m^{-3}$ $s^{-1}$','$mmol$ $N$ $m^{-2}$ $s^{-1}$','$mmol$ $N$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'SP_NH4_UPTAKE','DIAT_NH4_UPTAKE','DIAZ_NH4_UPTAKE','AMMOX','ANAMMOX',...
														'DON_REMIN','DONr_REMIN','ZOO_LOSS_DIC','SP_LOSS_DIC','DIAT_LOSS_DIC',...
														'DIAZ_LOSS_DIC','SP_GRAZE_DIC','DIAT_GRAZE_DIC','DIAZ_GRAZE_DIC',...
														'POC_REMIN'}; 
					obj.budget.(vars{i}).info.smseq  = [(-1) (-1) (-1) (-1) (-1) (1) (1) (param.Q) (param.Q) (param.Q) (param.Q) ...
														(param.Q) (param.Q) (param.Q) (1-param.DONrefract)];
					obj.budget.(vars{i}).info.fluxes = {[]};
					obj.budget.(vars{i}).info.lvls   = {[]};	
				elseif strcmp(vars{i},'N2')
					obj.budget.(vars{i}).info.tits   = {'$N_2$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2$ $m^{-3}$ $s^{-1}$','$mmol$ $N$ $m^{-3}$ $s^{-1}$'}; 
					obj.budget.(vars{i}).info.rates  = {'DENITRIF3','ANAMMOX'};
					obj.budget.(vars{i}).info.smseq  = [(1) (1)]; % changing to 2, both rates take 2N --> 2N
					%obj.budget.(vars{i}).info.fluxes = {'FG_N2'};
					%obj.budget.(vars{i}).info.lvls   = {  'sfc'};
					obj.budget.(vars{i}).info.fluxes = {[]};
					obj.budget.(vars{i}).info.lvls   = {[]};
				elseif strcmp(vars{i},'N2O')
					obj.budget.(vars{i}).info.tits   = {'$N_2O$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-2}$ $s^{-1}$','$mmol$ $N_2O$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
					obj.budget.(vars{i}).info.smseq  = [(0.5) (1) (-1)]; 
					obj.budget.(vars{i}).info.fluxes = {'FG_N2O'};
					obj.budget.(vars{i}).info.lvls   = {   'sfc'};
				elseif strcmp(vars{i},'N2O_SODEN')
					obj.budget.(vars{i}).info.tits   = {'$N_2O_{den}$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-2}$ $s^{-1}$','$mmol$ $N_2O$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'DENITRIF2','N2OSODEN_CONS'};
					obj.budget.(vars{i}).info.smseq  = [(0.5) (-1)]; 
					obj.budget.(vars{i}).info.fluxes = {'FG_N2O_SODEN'};
					obj.budget.(vars{i}).info.lvls   = {   'sfc'};
				elseif strcmp(vars{i},'N2O_AO1')
					obj.budget.(vars{i}).info.tits   = {'$N_2O_{nit}$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-2}$ $s^{-1}$','$mmol$ $N_2O$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'N2OAMMOX','N2OAO1_CONS'};
					obj.budget.(vars{i}).info.smseq  = [(1) (-1)]; 
					obj.budget.(vars{i}).info.fluxes = {'FG_N2O_AO1'};
					obj.budget.(vars{i}).info.lvls   = {   'sfc'};		
				elseif strcmp(vars{i},'N2O_ATM')
					obj.budget.(vars{i}).info.tits   = {'$N_2O_{atm}$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-2}$ $s^{-1}$','$mmol$ $N_2O$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'N2OATM_CONS'};
					obj.budget.(vars{i}).info.smseq  = [(-1)]; 
					obj.budget.(vars{i}).info.fluxes = {'FG_N2O_ATM'};
					obj.budget.(vars{i}).info.lvls   = {   'sfc'};
				elseif strcmp(vars{i},'N2O_SIDEN')
					obj.budget.(vars{i}).info.tits   = {'$N_2O_{bou}$'};
					obj.budget.(vars{i}).info.units  = {'$mmol$ $N_2O$ $m^{-3}$ $s^{-1}$','$mmol$ $N_2O$ $m^{-2}$ $s^{-1}$','$mmol$ $N_2O$ $s^{-1}$'};; 
					obj.budget.(vars{i}).info.rates  = {'N2OSIDEN_CONS'};
					obj.budget.(vars{i}).info.smseq  = [(-1)]; 
					obj.budget.(vars{i}).info.fluxes = {'FG_N2O_SIDEN'};
					obj.budget.(vars{i}).info.lvls   = {   'sfc'};
				end
				
				% Load all 3D output terms
				terms = [vars{i},obj.budget.(vars{i}).info.rates,obj.budget.(vars{i}).info.fluxes];
				terms = [terms(find(~cellfun(@isempty,terms)))];
				obj = loadData(obj,terms,file,'type','avg');

				% Integrate 3D variables, rates vertically
				terms = [vars{i},obj.budget.(vars{i}).info.rates];
				terms = [terms(find(~cellfun(@isempty,terms)))];
				obj = intVar(obj,terms);

				% Load and process 2D fluxes
				obj = getFluxes(obj,vars(i),obj.budget.(vars{i}).info.fluxes,obj.budget.(vars{i}).info.lvls);

				% Get dCdt, advection, sms, net terms
				obj = computeDcDt(obj,vars(i),hisfile);
				obj = computeXYZflux(obj,vars(i),physfile);
				obj = computeSMS(obj,vars(i),obj.budget.(vars{i}).info.rates,obj.budget.(vars{i}).info.smseq);
				obj = computeNet(obj,vars(i));

				% OPTIONAL SWITCHES
				% Integrate vertically (and horizontally)
				if A.int == 1
					obj = intBudg(obj,vars(i)); 
				end
				% Split sources and sinks
				if A.srcsnk == 1
					obj = sources_sinks(obj,vars(i));
				end
				% Save average, clear loaded data
				if A.clean == 1
					tmp = obj.data.avg.(vars{i});
					obj.data = [];
					obj.data.avg.(vars{i}) = tmp;
				end
			end
		end % end method getBudg

		%--------------------------------------------------------------------------------
		function obj = getFluxes(obj,vars,fluxes,lvls,varargin);
			% -------------------
			% Grab 2D fluxes for budget, convert to 3D 
			% Called in getBudg
			% Fluxes in mmol/m2/s
			%
			% Usage:
			% - obj = getFluxes(obj,vars,fluxes,lvls)
			%
			% Inputs:
			% - vars = BGC budget that you are closing (i.e. 'NO2')
			% - fluxes  = 2D fluxes (air-sea, sediment, etc)
			% - lvls    = levels to apply 2D flux ('sfc','sed', as a cell array)
			%
			% Optional inputs:
			% - type    = file type (avg, his, rst)
			%
			% Example:
			% - obj = getFluxes(obj,'N2O',{'FG_N2O'},{'sfc'})
			% -------------------
			disp('---------------------------------');
			disp('Get 2D fluxes, apply to 3D grid')

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);
			
			% Kill if no fluxes
			if isempty(fluxes{1});
				obj.budget.(vars{1}).fg = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
				obj.budget.(vars{1}).sed = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
				return
			end

			% Convert 2D flux to 3D based on lvls
			for i = 1:length(fluxes)
			
				% Initialize matrices-to-fill
				obj.budget.(vars{1}).fg = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
				obj.budget.(vars{1}).sed = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
	
				% Apply 2D mask to 2D data
				obj.data.(A.type).(fluxes{i}).data = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
				
				% Apply 2D flux to correct z-level to make 3D
				if strcmp(lvls{i},'sfc')
					tmpfg = zeros(size(obj.grid.(A.type).mask_rho3d));
					% Apply value into 3D grid
					tmpfg(:,:,obj.grid.nz,:) = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
					% Divide by z, save as 3D rate
					obj.budget.(vars{1}).fg = tmpfg ./ obj.grid.(A.type).Hz;
					% Mask padding
					obj.budget.(vars{1}).fg(end,:,:,:) = nan;
					obj.budget.(vars{1}).fg(:,end,:,:) = nan;
				elseif strcmp(lvls{i},'sed')
					tmpsed = zeros(size(obj.grid.(A.type).mask_rho3d));
					% Apply value into 3D grid
					tmpsed(:,:,1,:) = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
					% Divide by z, save as 3D rate
					obj.budget.(vars{1}).sed = tmpsed ./ obj.grid.(A.type).Hz;
					% Mask padding
					obj.budget.(vars{1}).sed(end,:,:,:) = nan;
					obj.budget.(vars{1}).sed(:,end,:,:) = nan;
				end
			end
		end % end method getFluxes

		%--------------------------------------------------------------------------------
		function obj = computeDcDt(obj,vars,file)
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
            % - file = (#s) load specific time-averages 
            %        = 0 load all files and average (for monthly)
			%
			% Example:
			% - obj = computeDcDt(obj,'NO2',1);
			% -------------------
			disp('---------------------------------');
			disp('Computing dC/dt')
		
			% Get file
			fname = [obj.paths.his,obj.info.his.files{file}];
			% Load data on either side of 'history' (first,last snapshot)
			for v = 1:length(vars)
				% Get start/finish from history
				his  = double(ncread(fname,vars{v},[obj.info.his.ind4D(1,:)],[obj.info.his.ind4D(2,:)]));
				% Divide by dt, apply 3d mask
				obj.budget.(vars{v}).dcdt = ((his(:,:,:,2:end) - his(:,:,:,1:end-1))/obj.info.his.dt) .* obj.grid.avg.mask_rho3d;

				% Mask padding
				obj.budget.(vars{v}).dcdt(end,:,:,:) = nan;
				obj.budget.(vars{v}).dcdt(:,end,:,:) = nan;
			end 
		end % end method computeDcDt

		%--------------------------------------------------------------------------------
		function obj = computeXYZflux(obj,vars,file)
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
			% - vars = variable to get flux terms for, as a cell array
			%
			% Example:
			% - obj = computeXYZflux(obj,{'NO2'});
			% -------------------
			disp('---------------------------------');
			disp('Get advective/diffusion terms')
		
            % Process optional inputs

            % Get file
            fname = [obj.paths.phys_flux,obj.info.phys_flux.files{file}];

			% Cycle through and load XYfluxes
			for i = 1:length(vars)
				% Load XYZ advection
				temp.X  = double(ncread(fname,['HorXAdvFlux_',vars{i}],[obj.info.phys_flux.ind4D(1,:)], [obj.info.phys_flux.ind4D(2,:)]));
				temp.Y  = double(ncread(fname,['HorYAdvFlux_',vars{i}],[obj.info.phys_flux.ind4D(1,:)], [obj.info.phys_flux.ind4D(2,:)]));
				temp.Z  = double(ncread(fname,['VertAdvFlux_',vars{i}],[obj.info.phys_flux.ind4D(1,:)], [obj.info.phys_flux.ind4D(2,:)]));
				temp.Zd = double(ncread(fname,['VertDiffFlux_',vars{i}],[obj.info.phys_flux.ind4D(1,:)],[obj.info.phys_flux.ind4D(2,:)]));
				
				% Add extra lons/lats for calc
				pad_x  = nan(1,size(temp.X,2),size(temp.X,3),size(temp.X,4));
				pad_y  = nan(size(temp.Y,1),1,size(temp.Y,3),size(temp.Y,4));
				temp.X = cat(1,temp.X,pad_x);
				temp.Y = cat(2,temp.Y,pad_y);

				% Initiate matrices-to-fill, compute adv X and Y
				% Converts from mmol/s to mmol/m3/s by dividing by grid area and height of cell (volume)
				% Get dimensions
				nx = obj.grid.nx;
				ny = obj.grid.ny;
				nz = obj.grid.nz;
				nt = obj.info.phys_flux.time(file);
				
				% X advection
				adx(1:nx,1:ny,1:nz,1:nt) = NaN;
				adx(1:nx,:,:,1:nt)       = (temp.X(1:nx,:,:,:) - temp.X(2:nx+1,:,:,:)) ./ ...
											obj.grid.avg.area3d(1:nx,:,:,:) ./ obj.grid.avg.Hz(1:nx,:,:,:);
				
				% Y advection
				ady(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				ady(:,1:ny,:,:)          = (temp.Y(:,1:ny,:,:) - temp.Y(:,2:ny+1,:,:)) ./ ...
											obj.grid.avg.area3d(:,1:ny,:,:) ./ obj.grid.avg.Hz(:,1:ny,:,:);
				
				% Z advection
				adz(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				adz(:,:,:,:)             = (temp.Z(:,:,1:nz,:) - temp.Z(:,:,2:nz+1,:)) ./ obj.grid.avg.Hz(:,:,:,:);
				
				% Z diffusion
				dfz(1:nx,1:ny,1:nz,1:nt) = NaN; 		
				dfz(:,:,:,:)             = (temp.Zd(:,:,1:nz,:) - temp.Zd(:,:,2:nz+1,:)) ./ obj.grid.avg.Hz(:,:,:,:);
				
				% Apply 3D mask
				obj.budget.(vars{i}).adx = adx .* obj.grid.avg.mask_rho3d(:,:,:,:);
				obj.budget.(vars{i}).ady = ady .* obj.grid.avg.mask_rho3d(:,:,:,:);
				obj.budget.(vars{i}).adz = adz .* obj.grid.avg.mask_rho3d(:,:,:,:);
				obj.budget.(vars{i}).dfz = dfz .* obj.grid.avg.mask_rho3d(:,:,:,:);

				% Mask padding
				obj.budget.(vars{i}).adx(end,:,:,:) = nan;
				obj.budget.(vars{i}).ady(end,:,:,:) = nan;
				obj.budget.(vars{i}).adz(end,:,:,:) = nan;
				obj.budget.(vars{i}).dfz(end,:,:,:) = nan;
				obj.budget.(vars{i}).adx(:,end,:,:) = nan;
				obj.budget.(vars{i}).ady(:,end,:,:) = nan;
				obj.budget.(vars{i}).adz(:,end,:,:) = nan;
				obj.budget.(vars{i}).dfz(:,end,:,:) = nan;
			end
		end % end method computeXYZFlux

		%--------------------------------------------------------------------------------
		function obj = computeSMS(obj,vars,rates,smseq,varargin)
			% ---------------------
			% Gathers sources and sinks
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeSMS(obj,rates,smseq)
			%
			% Inputs:
			% - vars = budget variable (i.e. 'NO2')
			% - rates = BGC rates, set in getBudg
			% - smseq = S-M-S equation, set in getBudg
			%
			% Optional:
			% - type = file type, defaults to 'avg' for budget calcs
			%
			% Example:
			% - vars = 'N2O_AO1';
			% - rates = {'N2OAMMOX','N2OAO1_CONS');
			% - smseq = [(1) (-1)];
			% - obj = computeSMS(obj,vars,rates,smseq);
			% ---------------------
			disp('---------------------------------');
			disp('Computing sources-minus-sinks (SMS)');

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

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
			obj.budget.(vars{1}).info.sms_eq = eq;
			obj.budget.(vars{1}).info.sms_factors = smseq;

			% Get SMS 
			dims = ndims(obj.data.(A.type).(rates{1}).data);
			eq = [];
			for i = 1:length(rates)
				eq{i} = [(smseq(i).*obj.data.(A.type).(rates{i}).data)];
			end
			tmpsms = sum(cat(dims+1,eq{:}),dims+1);
			obj.budget.(vars{1}).sms = tmpsms;

			% Get production 
			ind = find(smseq>0);
			if ~isempty(ind);
				tmprates = rates(ind);
				obj.budget.(vars{1}).info.prod_sources = tmprates;
				obj.budget.(vars{1}).info.prod_factors = smseq(ind);
				eq = [];
				for i = 1:length(tmprates)
					eq{i} = [(smseq(ind(i)).*obj.data.(A.type).(tmprates{i}).data)];
				end
				tmpprod = sum(cat(dims+1,eq{:}),dims+1);
				obj.budget.(vars{1}).prod = tmpprod;
			else
				obj.budget.(vars{1}).info.prod_sources = [];
				obj.budget.(vars{1}).info.prod_factors = []; 
				obj.budget.(vars{1}).prod = zeros([obj.grid.ndim_xyz obj.grid.nt]);
			end

			% Get consumption 
			ind = find(smseq<0);
			if ~isempty(ind)
				tmprates = rates(ind);
				obj.budget.(vars{1}).info.cons_sources = tmprates;
				obj.budget.(vars{1}).info.cons_factors = smseq(ind);
				eq = [];
				for i = 1:length(tmprates)
					eq{i} = [(smseq(ind(i)).*obj.data.(A.type).(tmprates{i}).data)];
				end
				tmpcons = sum(cat(dims+1,eq{:}),dims+1);
				obj.budget.(vars{1}).cons = tmpcons;
			else
				obj.budget.(vars{1}).info.cons_sources = [];
				obj.budget.(vars{1}).info.cons_factors = []; 
				obj.budget.(vars{1}).cons = zeros([obj.grid.ndim_xyz obj.grid.nt]);
			end

			% mask padding
			obj.budget.(vars{1}).sms(end,:,:,:) = nan;
			obj.budget.(vars{1}).prod(end,:,:,:) = nan;
			obj.budget.(vars{1}).cons(end,:,:,:) = nan;
			obj.budget.(vars{1}).sms(:,end,:,:) = nan;
			obj.budget.(vars{1}).prod(:,end,:,:) = nan;
			obj.budget.(vars{1}).cons(:,end,:,:) = nan;
		end % end method getSMS

		%--------------------------------------------------------------------------------
		function obj = computeNet(obj,vars)
			% ---------------------
			% Computes remainder (net) from budget equation
			% dC/dt = adv + dfz + sms + fg
			% Called in getBudg
			% End result is mmol/m3/s
			%
			% Usage:
			% - obj = computeNet(obj,vars)
			%
			% Inputs:
			% - vars = budget to close (i.e 'NO2');
			%
			% Example:
			% - obj = computeNet(obj,'NO2');
			% ---------------------
			disp('---------------------------------');
			disp('Computing budget remainder (net)');

			% Calculate remainder (net)
			obj.budget.(vars{1}).net = obj.budget.(vars{1}).dcdt - (obj.budget.(vars{1}).adx + ...
										    obj.budget.(vars{1}).ady  +  obj.budget.(vars{1}).adz + ...
										    obj.budget.(vars{1}).dfz  +  obj.budget.(vars{1}).sms + ...
										    obj.budget.(vars{1}).fg   +  obj.budget.(vars{1}).sed);
		end % end method computeNet
	
		%--------------------------------------------------------------------------------
		function obj = intVar(obj,vars,varargin)
			% ------------------
			% Vertically integrate 3D variable(s) 
			%
			% Usage:
			% - obj = intVar(obj,vars,varargin)
			%
			% Inputs:
			% - vars = 3D variables to integrate, as a cell array
			%
			% Optional:
			% - type = file type ('avg','his', etc)
			%
			% Example:
			% - obj = loadData(obj,{'NO3'},file);
			% - obj = intVar(obj,{'NO3'});
			% ------------------
			disp('---------------------------------');
			disp('Integrating 3D variables');

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

			% Load depths?
			try; obj.grid.(A.type).Hz;
			catch; disp('call loadDepth first'); return
			end

			% Go through each 3D rate, integrate vertically (and totally)
			for i = 1:length(vars)		
				tmpdata = obj.data.(A.type).(vars{i}).data .* obj.grid.(A.type).Hz .* obj.grid.(A.type).mask_rho3d; % mmol/m3/s --> mmol/m2/s
				tmpdata = squeeze(nansum(tmpdata,3));
				tmpdata = tmpdata .* obj.grid.mask_rho;
				obj.data.(A.type).(vars{i}).int = tmpdata; 
				obj.data.(A.type).(vars{i}).intunits = strrep(obj.data.(A.type).(vars{i}).units,'$m^{-3}$','$m^{-2}$');
				for t = 1:size(obj.grid.(A.type).Hz,4);
					obj.data.(A.type).(vars{i}).tot(t) = nansum(obj.data.(A.type).(vars{i}).int(:,:,t) .*obj.grid.area_rho,'all');
				end
				obj.data.(A.type).(vars{i}).totunits = strrep(obj.data.(A.type).(vars{i}).units,' $m^{-3}$','');
			end
		end % end method intVar

		%--------------------------------------------------------------------------------
		function obj = intFlux(obj,vars,varargin)
			% ------------------
			% Totally integrate 2D variable(s) 
			%
			% Usage:
			% - obj = intFlux(obj,vars,varargin)
			%
			% Inputs:
			% - vars = 2D variables to integrate, as a cell array
			%
			% Optional:
			% - type = file type ('avg','his',etc)
			%
			% Example:
			% - obj = loadData(obj,{'FG_N2O'},file);
			% - obj = intVar(obj,{'FG_N2O'});
			% ------------------
			disp('---------------------------------');
			disp('Integrating 2D variables');
		
            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);
	
			% Go through each 3D rate, integrate vertically (and totally)
			for i = 1:length(vars)		
				obj.data.(A.type).(vars{i}).int = obj.data.(A.type).(vars{i}).data;
				for t = 1:size(obj.data.(A.type).(vars{i}).data,3);
					obj.data.(A.type).(vars{i}).tot(t) = nansum(obj.data.(A.type).(vars{i}).data(:,:,t) .*obj.grid.area_rho,'all');
				end
			end
		end % end method intVar

		%--------------------------------------------------------------------------------
		function obj = intBudg(obj,vars,terms)
			% ------------------
			% Vertically integrate budget terms
			% Called in getBudg
			% Terms are all in mmol/m3/s, get it in mmol/m2/s by using dz
			%
			% Usage:
			% - obj = intBudg(obj)
			%
			% Inputs:
			% - vars = budget variable (i.e. 'NO2');
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
				disp('Must supply vars');
				return
			end

			disp('---------------------------------');
			disp('Computing vertical integration...');
			disp('...and getting grid area fluxes');
	
			% Grab terms, including remainder
			for i = 1:length(terms)
				if ~isempty(obj.budget.(vars{1}).(terms{i}))
					eval([terms{i},' = obj.budget.(vars{1}).',terms{i},' .* obj.grid.avg.mask_rho3d;']);
				end
			end

			% Integrate vertically (fg term should match real flux)
			% ...mmol/m3/s to mmol/m2/s
			for i = 1:length(terms)
				if ~isempty(obj.budget.(vars{1}).(terms{i}))
					eval(['obj.budget.(vars{1}).int',terms{i},' = squeeze(nansum(',terms{i},'.*obj.grid.avg.Hz,3));']);
				end
			end
			
			% ...mmol/m2/s to mmol/s
			for i = 1:length(terms);
				if ~isempty(obj.budget.(vars{1}).(terms{i}))
					for t = 1:size(obj.budget.(vars{1}).(terms{i}),4);
						eval(['obj.budget.(vars{1}).tot',terms{i},'(t)  = nansum(obj.budget.(vars{1}).int',terms{i},...
							  '(:,:,t) .*obj.grid.area_rho,''all'');']); 
					end
				end
			end

			% Check for total advection
			if sum(ismember({'adx','ady','adz'},terms)) == 3
				obj.budget.(vars{1}).intadv  =  obj.budget.(vars{1}).intadx + ...
												obj.budget.(vars{1}).intady + ...
										    	obj.budget.(vars{1}).intadz;
				obj.budget.(vars{1}).totadv  =	obj.budget.(vars{1}).totadx + ...
												obj.budget.(vars{1}).totady + ...
												obj.budget.(vars{1}).totadz;
			end
		end

		function [obj] = sources_sinks(obj,vars)
			% ------------------------------------------------
			% Function to extract sources and sinks
			%
			% Usage: 
			% - [obj] = sources_sinks(obj,vars) 
			% 
			% Inputs:
			% - obj:    romsMaster object
			% - vars:   budget variable
			%
			% Example:
			% - [obj] = sources_sinks(obj,{'NH4'});
			% ------------------------------------------------

			% Get sources
			for i = 1:length(vars)
				for j = 1:length(obj.budget.(vars{i}).info.prod_sources);
					this_var = obj.budget.(vars{i}).info.prod_sources{j};
					this_factor = obj.budget.(vars{i}).info.prod_factors(j);
					tmp.sources.(this_var) = obj.data.avg.(this_var).data.*this_factor;
				end
				for j = 1:length(obj.budget.(vars{i}).info.cons_sources);
					this_var = obj.budget.(vars{i}).info.cons_sources{j};
					this_factor = obj.budget.(vars{i}).info.cons_factors(j);
					tmp.sinks.(this_var) = obj.data.avg.(this_var).data.*this_factor;
				end
				obj.budget.(vars{i}).sources = tmp.sources;
				obj.budget.(vars{i}).sinks   = tmp.sinks;
			end
		end % end function sources_and_sinks
	
		%--------------------------------------------------------------------------------
		function plotBudg(obj,vars,varargin);
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
			A.lonbounds = [floor(min(obj.grid.lon_rho(:))) ceil(max(obj.grid.lon_rho(:)))];
			A.latbounds = [floor(min(obj.grid.lat_rho(:))) ceil(max(obj.grid.lat_rho(:)))];
			A.tfont     = 18;
			A.cbfont    = 14;
			A           = parse_pv_pairs(A,varargin); % parse method arguments to A
			
			% Process inputs
			if isempty(A.time)
				A.time = 1:size(obj.budget.(vars{1}).intdcdt,3);
			end

			% First generate uniform color bars for each axis	
			% Go through each variables
			terms = {'dcdt','adv','dfz','sms','fg','sed','net'};
			for i = 1:length(vars)
				for j = 1:length(terms)
					if ~strcmp(terms{j},'adv')
						tmp.(terms{j}) = obj.budget.(vars{i}).(['int',terms{j}]);
					else
						tmp.adv = [obj.budget.(vars{i}).intadx + ...
								   obj.budget.(vars{i}).intady + ... 
								   obj.budget.(vars{i}).intadz];	
					end
					lims = [];
					for t = 1:length(A.time)
						lims{t} = prclims(tmp.(terms{j})(:,:,t),'prc',A.prc,'bal',1);
					end
					lims = cat(1,lims{:});
					lims = max(abs(lims(:)));
					lims = [-lims lims];
					% Plot results
					for t = 1:length(A.time)
						% Plot results
						dat  = tmp.(terms{j})(:,:,t);
						[fig,cb] = mapPlot(obj,dat,'levels',linspace(lims(1),lims(end),31),'cmap',cmocean('-balance',31));
						hold on
						m_plot(obj.grid.lon_rho(1,:),obj.grid.lat_rho(1,:),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(:,1),obj.grid.lat_rho(:,1),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'--k','linewidth',2);
						title([terms{j},': time ',num2str(t)],...
							  'Interpreter', 'Latex','FontSize',A.tfont);
						ylabel(cb,obj.budget.(vars{1}).info.units{2},'Interpreter','Latex');
						fname = [vars{i},'_',terms{j}];
						if t == 1
							if exist([obj.paths.plots.budget,fname,'.pdf']) == 2
								cmd = ['rm ',obj.paths.plots.budget,fname,'.pdf'];
								system(cmd);
							end
							export_fig('-pdf',[obj.paths.plots.budget,fname]);
						else
							export_fig('-pdf',[obj.paths.plots.budget,fname],'-append');
						end
						close(fig)	
					end
				end
			end
		end % end method plotBudg

		%--------------------------------------------------------------------------------
		function obj = loadData(obj,vars,file,varargin)
			% -------------------
			% Method to Load ROMS data.
			%
			% Usage:
			% - obj = loadData(obj,vars,file,varargin)
			%
			% Inputs:
			% - vars   = variable(s) to load, as a cell array
			% - file   = (#s) time-average file to load
			%          = 0 for load all and average
			%
			% Inputs (varargin):
			% - type    = 'avg','his','rst','flux' 
			% - time    = option to load speciific time stamp from file
			% 
			% Examples:
			% - obj = loadData(obj)
			% - obj = loadData(obj,{'temp','salt'},1,'type','avg');
			% -------------------
			disp(['Loading ROMS variables']);
			
			% defaults for optional  arguments
			A.type   = ['avg']; % input file type
			A.time   = [];      % optional time
			A        = parse_pv_pairs(A,varargin); % parse method arguments to A
	
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
		
            % Initialize matrices to fill
            nx = obj.info.(A.type).xi_rho;  nx = nx(end)-nx(1)+1;
            ny = obj.info.(A.type).eta_rho; ny = ny(end)-ny(1)+1;
			nz = obj.info.(A.type).s_rho;
			nt = sum(obj.info.(A.type).time(file));
			if strcmp(vars{i},'u');
				nx = nx-1;
			elseif strcmp(vars{i},'v');
				ny = ny-1;
			end

			% Check for depth
			if isempty(obj.grid.(A.type)) | ~nt == size(obj.grid.(A.type).z_r) 
				obj = loadDepth(obj,file,'type',A.type);
			end

			% Get files
			files = obj.info.(A.type).files(file);
			
			% Blank
			for i = 1:length(vars)
				obj.data.(A.type).(vars{i}).data = nan(nx,ny,nz,nt);
			end
		
			% Start variable loop	
			for i = 1:length(vars)
				% 2D check
				idx = [];
				idx = find(ismember(obj.info.(A.type).var2d,vars{i})==1);
				if ~isempty(idx)
					for ff = 1:length(files)
						% Load data
						tmp.data = ncread([obj.paths.(A.type),files{ff}],vars{i},...
										  [obj.info.(A.type).ind3D(1,:)],[obj.info.(A.type).ind3D(2,:)]);
						if strcmp(vars{i},'u');
							out.data{ff} = tmp.data .* obj.grid.mask_u;
						elseif strcmp(vars{i},'v');
							out.data{ff} = tmp.data .* obj.grid.mask_v;
						else
							out.data{ff} = tmp.data .* obj.grid.mask_rho;
						end
					end
					obj.data.(A.type).(vars{i}).data = cat(3,out.data{:});
					clear out
					% Grab variable info
					obj.data.(A.type).(vars{i}).name = obj.info.(A.type).name2d{idx};
					obj.data.(A.type).(vars{i}).units = obj.info.(A.type).unit2d{idx};
					continue
				end
				% 3D check
				idx = find(ismember(obj.info.(A.type).var3d,vars{i})==1);
				if ~isempty(idx);
                    for ff = 1:length(files)
                        % Load data
                        tmp.data = ncread([obj.paths.(A.type),files{ff}],vars{i},...
                                          [obj.info.(A.type).ind4D(1,:)],[obj.info.(A.type).ind4D(2,:)]);
						if strcmp(vars{i},'u');
							out.data{ff} = tmp.data .* obj.grid.mask_u;
						elseif strcmp(vars{i},'v');
							out.data{ff} = tmp.data .* obj.grid.mask_v;
						else
							out.data{ff} = tmp.data .* obj.grid.mask_rho;
						end
                    end
					obj.data.(A.type).(vars{i}).data = cat(4,out.data{:});
					clear out
					% Grab variable info
					obj.data.(A.type).(vars{i}).name = obj.info.(A.type).name3d{idx};
					obj.data.(A.type).(vars{i}).units = obj.info.(A.type).unit3d{idx};
					continue
				end
				% Derived check
				if isempty(idx)
					try
						obj = computeVar(obj,vars(i),file,'type',A.type);
					catch
						disp(['ERROR: Unreckognizable variable ',vars{i}]);
					end
					% Force end of loop
					continue
				end
			end
		end % end method loadData
	
		%--------------------------------------------------------------------------------
		function obj = sliceROMS(obj,vars,choice,deg,file,varargin);
			% -------------------
			% Takes 2D depth slice of ROMS along a given latitude or longitude
			% Options are set by user
			% 
			% Usage:
			% - obj = sliceROMS(obj,vars,choice,deg,file,varargin);
			% 
			% Inputs:
			% - vars   = ROMS variable(s) to slice, as a cell array
			% - choice = 'lon','lat','xi','eta' 
			%			  (lat/lon slices along a given lat/lon degree)
			%			  (xi/eta slices along a given xi or eta index, use gridView)
			% - deg    = lon/lat degree(s) or x/y indicies
			% - file   = file number to slice
			%
			% Inputs (varargin):
			% - type     = 'avg','his','rst' (or z_avg etc. for obj.grid.z_avg_dep levels) 
			% - zlim     = if using 'z_avg' or others, limit max obj.grid.z_avg_dep levels) 
			% 
			% Example
			% - obj = sliceROMS(obj,{'temp','salt'},'lon',0,1);
			% 
			% This will slice temp and salt data at 0 degrees longitude
			% -------------------

			% Clear slice struct
			obj.slice = [];
	
			% Check inputs
			if nargin<5
				disp('Incorrect number of inputs');
				help sliceROMS
				return
			end	

			% Grab user inputs
			A.type     = ['avg'];
			A.zlim     = [inf];
			A          = parse_pv_pairs(A,varargin);

			% If calling lat/lon, force z_avg/z_his
			if strcmp(choice,'lat') | strcmp(choice,'lon');
				if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
					A.type = ['z_',A.type];
				end
			end

			% Load raw data, or slice it first
			if strcmp(A.type,'z_avg') | strcmp(A.type,'z_his') | strcmp(A.type,'z_phys_flux');
				% If requesting z_avg/z_his, zslice data first
				obj = zslice(obj,vars,obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim),file,'type',A.type(3:end));
				nz  = length(obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim));
			else
				% Use raw output data
				nz = obj.grid.nz;
				% Check if grid is loaded
				try; obj.grid.(A.type).z_r;
				catch; obj = loadDepth(obj,file,'type',A.type);
				end
				% Load variables
				for i = 1:length(vars)
					if ismember(vars(i),obj.info.(A.type).var3d) | ismember(vars(i),obj.info.(A.type).var2d)
						obj = loadData(obj,vars(i),file,'type',A.type);
					else
						obj = computeVar(obj,vars(i),file,'type',A.type);
					end
				end
			end

			% Get dimensions
			if strcmp(choice,'lat') | strcmp(choice,'eta');
				dmsn = obj.grid.nx;
			elseif strcmp(choice,'lon') | strcmp(choice,'xi');
				dmsn = obj.grid.ny;
			end
			if strcmp(A.type,'z_avg') | strcmp(A.type,'z_his');
				nt = size(obj.grid.(A.type(3:end)).Hz,4);
			else
				nt = size(obj.grid.(A.type).Hz,4);
			end
			ns = length(deg);

			% Choose latitude, longitude, x-idx or y-idx
			% Latitude slice
			if strcmp(choice,'lat')
				lonlat = obj.grid.lat_rho;
				dstr   = ['latitude'];
			% Longitude slice
			elseif strcmp(choice,'lon')
				lonlat = obj.grid.lon_rho; 
				dstr   = ['longitude'];
			% X slice
			elseif strcmp(choice,'xi')
				if length(deg)>1
					order = [2 3 1 4];
				else
					order = [1 2 4 3];
				end
				for v = 1:length(vars)
					if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
						obj.data.(A.type).(vars{v}).slice  = squeeze(obj.data.(A.type).(vars{v}).data(deg,:,:,:));
						slicedepth = squeeze(obj.grid.(A.type).z_r(deg,:,:,:));
					else
						obj.data.(A.type(3:end)).(vars{v}).slice  = squeeze(obj.data.(A.type(3:end)).(vars{v}).slice(deg,:,:,:));
						slicedepth = squeeze(obj.slice.depth(deg,:,:,:)); 
					end
				end
				slicelon = squeeze(obj.grid.lon_rho(deg,:))';
				slicelat = squeeze(obj.grid.lat_rho(deg,:))';
			% Y slice
			elseif strcmp(choice,'eta');
				if length(deg)>1
					order = [1 3 2 4];
				else
					order = [1 2 4 3];
				end
                for v = 1:length(vars)
                    if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
                        obj.data.(A.type).(vars{v}).slice  = squeeze(obj.data.(A.type).(vars{v}).data(:,deg,:,:));
                        slicedepth = squeeze(obj.grid.(A.type).z_r(:,deg,:,:));
                    else
                        obj.data.(A.type(3:end)).(vars{v}).slice  = squeeze(obj.data.(A.type(3:end)).(vars{v}).slice(:,deg,:,:));
                        slicedepth = squeeze(obj.slice.depth(:,deg,:,:));
                    end
				end
				slicelon = squeeze(obj.grid.lon_rho(:,deg));
				slicelat = squeeze(obj.grid.lat_rho(:,deg));
			end

			% If xi-idx or eta-idx, organize slice output and end routine
			if strcmp(choice,'xi') | strcmp(choice,'eta');
				if ns>1
					for v = 1:length(vars)
						if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
							obj.data.(A.type).(vars{v}).slice = permute(obj.data.(A.type).(vars{v}).slice,order);
						else
							obj.data.(A.type(3:end)).(vars{v}).slice = permute(obj.data.(A.type(3:end)).(vars{v}).slice,order);
						end
					end
					obj.slice.depth = permute(slicedepth,order);
					obj.slice.lat   = permute(repmat(slicelat,[1 1 nz nt]),[1 3 2 4]);
					obj.slice.lon   = permute(repmat(slicelon,[1 1 nz nt]),[1 3 2 4]);
				else
					for v = 1:length(vars)
						if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
							obj.data.(A.type).(vars{v}).slice = permute(obj.data.(A.type).(vars{v}).slice,order);
						else
							obj.data.(A.type(3:end)).(vars{v}).slice = permute(obj.data.(A.type(3:end)).(vars{v}).slice,order);
						end
					end
					obj.slice.depth = permute(slicedepth,order);
					obj.slice.lat   = repmat(slicelat,[1 nz ns nt]);
					obj.slice.lon   = repmat(slicelon,[1 nz ns nt]);
				end
				if strcmp(choice,'xi');
					obj.slice.deg = obj.slice.lat;
					obj.slice.coord = 'latitude';
				elseif strcmp(choice,'eta');
					obj.slice.deg = obj.slice.lon;
					obj.slice.coord = 'longitude';
				end
				obj.slice.mask  = double(~isnan(obj.slice.depth));
				obj.slice.mask(obj.slice.mask==0) = NaN;
				% Kill routine	
				if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
					obj.data.(A.type).(vars{v}).data = [];
				else
					obj.data.(A.type(3:end)).(vars{v}).data = [];
				end
				return
			end

			% Dimensions to fill
			fillmat = [dmsn nz ns];	

			% Copy mask
			mask = obj.grid.mask_rho;
			mask(isnan(mask)) = 0;
			mask = repmat(mask,1,1,nz);
		
			% Organize 3D ROMS data and slice
			for ff = 1:length(vars)
				if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
					tmpdata{ff} = obj.data.(A.type).(vars{ff}).data;			
				else
					tmpdata{ff} = obj.data.(A.type(3:end)).(vars{ff}).slice;			
				end
			end
			data = cat(5,tmpdata{:});
			clear tmpdata;

			% Take slice of data for each month
			disp(' '); disp(['Slicing ROMS data @ ',num2str(deg),'deg ',dstr,'...']);disp(' ');
			dims = [dmsn nz nt ns obj.grid.nx obj.grid.ny 1];
			tmpslice = [];
			for i = 1:length(vars)
				if length(vars)==1
					tmpdata = data;
				else
					tmpdata = squeeze(data(:,:,:,:,i));	
				end
				[tmpslice{i},tmpmask{i},tmpdepth{i}] = romsMaster.lonlat_slice(dims,lonlat,tmpdata,mask,deg,...
					obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim),fillmat,i,length(vars));
			end 
			tmpdim = ndims(tmpslice{1});
			tmpslice = cat(tmpdim+1,tmpslice{:});
			tmpmask = tmpmask{1};
			tmpdepth = tmpdepth{1};

			% Get lon/lat data (outside parfor)
			tmpdeg  = NaN(dmsn,ns);
			for i = 1:dmsn
				for j = 1:ns	
					if dmsn == obj.grid.nx;
						tmpdeg(i,j)  = interp1(squeeze(obj.grid.lat_rho(i,:)),squeeze(obj.grid.lon_rho(i,:)),deg(j));
					elseif dmsn == obj.grid.ny;
						tmpdeg(i,j)  = interp1(squeeze(obj.grid.lon_rho(:,i)),squeeze(obj.grid.lat_rho(:,i)),deg(j));
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
					for rcrd = 1:nt
						if ns == 1
							tmp = squeeze(tmpslice(:,:,rcrd,ff));
							tmp(masktmp==0) = NaN;
							tmpslice(:,:,rcrd,ff) = tmp;
						else
							tmp = squeeze(tmpslice(:,:,rcrd,j,ff));
							tmp(masktmp==0) = NaN;
							tmpslice(:,:,rcrd,j,ff) = tmp;
						end
					end
				end
			end

			% Save results
			for ff = 1:length(vars);
				if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
					obj.data.(A.type).(vars{ff}).slice = [];
				else
					obj.data.(A.type(3:end)).(vars{ff}).slice = [];
				end
				if ns == 1
					tmpdat = squeeze(tmpslice(:,:,:,ff));
					if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
						obj.data.(A.type).(vars{ff}).slice  = permute(repmat(tmpdat,[1 1 1 ns]),[1 2 4 3]); 
					else
						obj.data.(A.type(3:end)).(vars{ff}).slice  = permute(repmat(tmpdat,[1 1 1 ns]),[1 2 4 3]); 
					end
					obj.slice.deg = permute(repmat(tmpdeg,[1 1 nt ns]),[1 2 4 3]);
				else
					if ~strcmp(A.type,'z_avg') & ~strcmp(A.type,'z_his');
						obj.data.(A.type).(vars{ff}).slice  = permute(squeeze(tmpslice(:,:,:,:,ff)),[1 2 4 3]);
					else
						obj.data.(A.type(3:end)).(vars{ff}).slice  = permute(squeeze(tmpslice(:,:,:,:,ff)),[1 2 4 3]);
					end
					obj.slice.deg = repmat(tmpdeg,[1 1 1 nt]); 
				end
				if dmsn == obj.grid.nx;
						obj.slice.coord = 'latitude';
				elseif dmsn == obj.grid.ny;
						obj.slice.coord = 'longitude';
				end
				obj.slice.depth = permute(repmat(tmpdepth,[1 1 nt ns]),[1 2 4 3]);
				obj.slice.sect  = deg;
				fprintf(['\n']);
			end
			obj.slice.mask  = double(~isnan(obj.slice.depth));
			obj.slice.mask(obj.slice.mask==0) = NaN;
		end % end method sliceROMS

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
			% - choice = 'lon','lat','xi', or 'eta' 
			%			  (lat/lon slices along a given lat/lon degree)
			%			  (x/y slices along a given x or y index (lon_rho or lat_rho))
			% - deg    = lon/lat degree or x/y index
			%
			% Inputs (varargin):
			% - zlim = depth limit to restrict output
			% 
			% Example
			% - obj = sliceDiag(obj,{'temp','salt'},'lon',0);
			% 
			% This will slice temp and salt data at 0 degrees longitude, and will also return
			% sliced temperature (but not salinity) data from the validation dataset
			% -------------------

			% Get optional inputs
			A.zlim = inf;
			A = parse_pv_pairs(A,varargin);

			% Grab longitude/latitude data
			lon  = obj.grid.lon_rho;
			lat  = obj.grid.lat_rho;

			% Choose latitude, longitude
			if strcmp(choice,'lat')
				lonlat = lat;
				dmsn   = obj.grid.nx; 
				nz     = length(obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim)); 
				nt     = 12; 
				ns     = length(deg); 
				dstr   = ['latitude'];
			elseif strcmp(choice,'lon')
				lonlat = lon; 
				dmsn   = obj.grid.ny; 
				nz     = length(obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim)); 
				nt     = 12; 
				ns     = length(deg); 
				dstr   = ['longitude'];
			end
			
			% Get lon/lat data (outside parfor)
			tmpdeg  = NaN(dmsn,ns);
			for i = 1:dmsn
				for j = 1:ns	
					if dmsn == obj.grid.nx;
						tmpdeg(i,j)  = interp1(squeeze(lat(i,:)),squeeze(lon(i,:)),deg(j));
					elseif dmsn == obj.grid.ny;
						tmpdeg(i,j)  = interp1(squeeze(lon(:,i)),squeeze(lat(:,i)),deg(j));
					end
				end
			end
			tmpdeg = repmat(tmpdeg,1,1,nz); tmpdeg = permute(tmpdeg,[1 3 2]);

			% Dimensions to fill
			fillmat  = [dmsn nz ns];
			tmpdepth = obj.grid.z_avg_dep(obj.grid.z_avg_dep<A.zlim)'.*ones(fillmat);

			% Perform slices of each variable
			for ff = 1:length(vars);
				obj.diag.(vars{ff}).slice = [];
				disp(' '); disp(['Slicing validation ',vars{ff},' data...']);disp(' ');
				% get path and coords for current variable
				if ~isfield(obj.paths.diag,(vars{ff}));
					disp([vars{ff},' is not a diagnostic variable']);
					disp(['Filling with NaN']);
					obj.diag.(vars{ff}).slice = nan(dmsn,nz,ns,nt);
                    obj.diag.(vars{ff}).name = 'null';
                    obj.diag.(vars{ff}).units = 'null';
					continue
				end
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
					else
						tmp.latg = tmp.lat;	
						tmp.long = tmp.lon;
						tmp.depthg = tmp.depth;
						tmp.lat = squeeze(tmp.lat(1,:,1,1));
						tmp.lon = squeeze(tmp.lon(:,1,1,1));
						tmp.depth = squeeze(tmp.depth(1,1,:,1));
					end
					
					% go through each requested slice
					for j = 1:ns
						% go through each time record
						for rcrd = 1:nt
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
						obj.diag.(vars{ff})(i).slice(:,:,:,j)  = cat(ndims(tmp.out{1})+1,tmp.out{:});
						obj.diag.(vars{ff})(i).units = curVar.units{i};
						obj.diag.(vars{ff})(i).name  = curVar.name{i};
					end % end j-loop
					% Check if slice is already filled
					if isempty(obj.slice)
						if dmsn == obj.grid.nx;
							obj.slice.coord = 'latitude';
						elseif dmsn == obj.grid.ny;
							obj.slice.coord = 'longitude';
						end
						fill_slice = 1;
					else
						fill_slice = 0;
					end
					
					if fill_slice	
						% Organize slice output
						obj.slice.deg(:,:,j,:) = permute(repmat(tmpdeg,[1 1 1 nt]),[1 2 4 3]);
						obj.slice.depth(:,:,j,:) = permute(repmat(tmpdepth,[1 1 1 nt]),[1 2 4 3]);
						obj.slice.sect(j) = deg;
					end
					% Finalize output
					obj.diag.(vars{ff})(i).slice = permute(obj.diag.(vars{ff})(i).slice,[1 2 4 3]);
				end % end i-loop
			end % end var-loop
		end % end method sliceDiag

		%--------------------------------------------------------------------------------
		function obj = computeVar(obj,vars,file,varargin)
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
			% - OW       (Okubo-Weiss)
			% - vort     (relative vorticity)
			% - sN       (normal component of strain)
			% - sS       (shear component of strain)
			% - Jden_N2O (SMS of N2O via denitrification)
			% - Jnit_N2O (SMS of N2O via nitrification)
			% - AOU      (Apparent oxygen utilization)
			% - DeltaN2O (Supersaturated N2O)
			% - OMZ      (OMZ thickness)
			% - NH4vNO2  (NH4 vs NO2 ratio)
			% - NO2vNH4  (NO2 vs NH4 ratio)
			% - NOX      (NO2+NO3)
			%
			% Usage:
			% - obj = computeVar(obj,vars,varargin)
			%
			% Inputs:
			% - var  = variable(s) to compute, as a cell array	
			% - file = file number to extract data from
			%
			% Inputs (varargin)
			% - type   = 'avg' (default), 'his'
			% - dep/ip = depths/isopycnal to slice Okubo-Weiss fields on
			% - thresh = thresholds, if needed
			%
			% Example:
			% - obj = computeVar(obj,{'sigma0'},1,'type','avg');
			% ------------------

			% process inputs
			A.type    = 'avg';
			A.ip      = [];
			A.dep     = [];
			A.thresh  = 20; %uM O2 for OMZ thickness
			A.zlim    = inf;
			A         = parse_pv_pairs(A,varargin);

			% dims for vertical computations
			nx = obj.grid.nx;
			ny = obj.grid.ny;
			nz = obj.grid.nz;
			nt = sum(obj.info.(A.type).time(file));

			% Check for depth
			if isempty(obj.grid.(A.type))
				obj = loadDepth(obj,file,'type',A.type);
			elseif ~size(obj.grid.(A.type).z_r,4) == nt
				obj = loadDepth(obj,file,'type',A.type);
			end
			zgrid = obj.grid.(A.type).z_r;

			% process requests
			for i = 1:length(vars)
				% pressure calc
				if strcmp(vars{i},'pres') % & ~isfield(obj.data,'pres');
					disp('Calculating sea pressure');
					obj.data.(A.type).pres.data = sw_pres(-obj.grid.(A.type).z_r,repmat(obj.grid.lat_rho,[1 1 nz nt]));
					obj.data.(A.type).pres.data = obj.data.(A.type).pres.data .* obj.grid.(A.type).mask_rho3d;
					obj.data.(A.type).pres.name  = 'averaged Pressure';
					obj.data.(A.type).pres.units = '$dbar$';
				% rho calc
				elseif strcmp(vars{i},'sigma')
					disp('Calculating sigma');
					obj = computeVar(obj,{'pres'},file,'type',A.type);
					if isfield(obj.data.(A.type),'temp');
						tt = 1;
						origtemp = obj.data.(A.type).temp;
					else
						tt = 0;
					end
					if isfield(obj.data.(A.type),'salt');
						ss = 1;
						origsalt = obj.data.(A.type).salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},file,'type',A.type);
					tmptemp = sw_temp(obj.data.(A.type).salt.data(:),obj.data.(A.type).temp.data(:),obj.data.(A.type).pres.data(:),0); 
					tmpsigma = sw_dens(obj.data.(A.type).salt.data(:),tmptemp,obj.data.(A.type).pres.data(:));
					obj.data.(A.type).sigma.data = reshape(tmpsigma,size(obj.data.(A.type).temp.data))-1000;
					obj.data.(A.type).sigma.data = real(obj.data.(A.type).sigma.data) .* obj.grid.(A.type).mask_rho3d;
					obj.data.(A.type).sigma.name  = 'averaged Density';
					obj.data.(A.type).sigma.units = '$kg$ $m^{-3}$';
					if tt == 1
						obj.data.(A.type).temp = origtemp;
					else
						obj.data.(A.type).temp = [];
					end
					if ss == 1
						obj.data.(A.type).salt = origsalt;
					else
						obj.data.(A.type).salt = [];
					end
				% bvf calc
				elseif strcmp(vars{i},'bvf') | strcmp(vars{i},'pv');
					disp('Calculating bvf (Brunt Vaisala)');
					obj = computeVar(obj,{'pres'},file,'type',A.type);
					if isfield(obj.data.(A.type),'temp');
						tt = 1;
						origtemp = obj.data.(A.type).temp;
					else
						tt = 0;
					end
					if isfield(obj.data.(A.type),'salt');
						ss = 1;
						origsalt = obj.data.(A.type).salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},file,'type',A.type);
					tmptemp = reshape(sw_temp(obj.data.(A.type).salt.data(:),obj.data.(A.type).temp.data(:),...
											  obj.data.(A.type).pres.data(:),0),size(obj.data.(A.type).salt.data)); 
					tmplat  = repmat(obj.grid.lat_rho,1,1,nz,nt);
					% Permute to get depth in front
					tmptemp = permute(obj.data.(A.type).temp.data,[3 1 2 4]);
					tmpsalt = permute(obj.data.(A.type).salt.data,[3 1 2 4]);
					tmppres = permute(obj.data.(A.type).pres.data,[3 1 2 4]);
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
					obj.data.(A.type).bvf.data  = bvf_rho;
					obj.data.(A.type).bvf.name  = 'averaged Brunt Vaisala frequency';
					obj.data.(A.type).bvf.units = '$s^{-2}$';
					obj.data.(A.type).pv.data   = pv_rho;
					obj.data.(A.type).pv.name   = 'averaged Potential Vorticity';
					obj.data.(A.type).pv.units  = '$m^{-1}$ $s^{-1}$';
					if tt == 1
						obj.data.(A.type).temp = origtemp;
					else
						obj.data.(A.type).temp = [];
					end
					if ss == 1
						obj.data.(A.type).salt = origsalt;
					else
						obj.data.(A.type).salt = [];
					end
					obj.data.(A.type).bvf.data = obj.data.(A.type).bvf.data .* obj.grid.(A.type).mask_rho3d;
					obj.data.(A.type).pv.data  = obj.data.(A.type).pv.data .* obj.grid.(A.type).mask_rho3d;
				% sigma0
				elseif strcmp(vars{i},'sigma0')
					disp('Calculating sigma0');
					obj = computeVar(obj,{'pres'},file,'type',A.type);
					if isfield(obj.data.(A.type),'temp');
						tt = 1;
						origtemp = obj.data.(A.type).temp;
					else
						tt = 0;
					end
					if isfield(obj.data.(A.type),'salt');
						ss = 1;
						origsalt = obj.data.(A.type).salt;
					else
						ss = 0;
					end
					obj = loadData(obj,{'temp','salt'},file,'type',A.type);
					tmp_lon      = repmat(obj.grid.lon_rho,1,1,nz,nt);
					tmp_lat      = repmat(obj.grid.lat_rho,1,1,nz,nt);
					tmp_salt     = obj.data.(A.type).salt.data;
					tmp_pres     = obj.data.(A.type).pres.data;
					tmp_salt_abs = romsMaster.getSA(tmp_salt,tmp_pres,tmp_lon,tmp_lat);
					tmp_theta    = gsw_CT_from_pt(tmp_salt_abs(:),obj.data.(A.type).temp.data(:));
					tmp_sigma0   = gsw_sigma0(tmp_salt_abs(:),tmp_theta(:));
					tmp_sigma0   = reshape(tmp_sigma0,nx,ny,nz,nt);
					obj.data.(A.type).sigma0.data  = tmp_sigma0;
					obj.data.(A.type).sigma0.data = obj.data.(A.type).sigma0.data .* obj.grid.(A.type).mask_rho3d;
					obj.data.(A.type).sigma0.name  = 'averaged Potential Density';
					obj.data.(A.type).sigma0.units = '$kg$ $m^{-3}$';
					if tt == 1
						obj.data.(A.type).temp = origtemp;
					else
						obj.data.(A.type).temp = [];
					end
					if ss == 1
						obj.data.(A.type).salt = origsalt;
					else
						obj.data.(A.type).salt = [];
					end
				% nstar
				elseif strcmp(vars{i},'nstar');
					obj = loadData(obj,{'NO3','NO2','PO4'},file,'type',A.type);
					tmpnstar = (obj.data.(A.type).NO3.data + obj.data.(A.type).NO2.data) - 16.*obj.data.(A.type).PO4.data + 2.9;
					obj.data.(A.type).nstar.data = tmpnstar;
					obj.data.(A.type).nstar.data = obj.data.(A.type).nstar.data .* obj.grid.(A.type).mask_rho3d;
					obj.data.(A.type).nstar.name = 'averaged N*';
					obj.data.(A.type).nstar.units = '$mmol$ $N$ $m^{-3}$';
				% NPP
				elseif strcmp(vars{i},'NPP') | strcmp(vars{i},'npp');
					obj     = loadData(obj,{'TOT_PROD'},file,'type',A.type);
					tmpprod = obj.data.(A.type).TOT_PROD.data;
					tmpHz   = obj.grid.(A.type).Hz; 
					tmpNPP  = squeeze(nansum(tmpprod.*tmpHz,3)).*3600*24*12;
					obj.data.(A.type).NPP.data  = tmpNPP;
					obj.data.(A.type).NPP.data  = obj.data.(A.type).NPP.data .* obj.grid.mask_rho;
					obj.data.(A.type).NPP.name  = 'averaged Net Primary Production (NPP)';
					obj.data.(A.type).NPP.units = '$mg$ $C$ $m^{-2}$ $d^{-1}$'; 
				% SSH	
				elseif strcmp(vars{i},'SSH') | strcmp(vars{i},'ssh');
					obj = loadData(obj,{'zeta'},file,'type',A.type);
					obj = loadDiag(obj,{'SSH'},0);
					slacorr = nanmedian(obj.diag.SSH.data(:)) - nanmedian(obj.data.(A.type).zeta.data(:));
					disp(' '); disp(['Adding correction of ',num2str(slacorr),'m to ROMS SSH']);
					obj.data.(A.type).SSH.data = obj.data.(A.type).zeta.data + slacorr;
					obj.data.(A.type).SSH.name = 'averaged sea-surface height';
					obj.data.(A.type).SSH.units = '$m$'; 
				% Wind Stress or Wind Stress Curl
				elseif strcmp(vars{i},'WS') | strcmp(vars{i},'ws') | strcmp(vars{i},'WSC') | strcmp(vars{i},'wsc');
					obj    = loadData(obj,{'sustr','svstr'},file,'type',A.type);
					tmpu   = obj.data.(A.type).sustr.data;
					tmpv   = obj.data.(A.type).svstr.data;
					tmpang = obj.grid.angle;
					[tmpws,tmpwsc] = romsMaster.WindStress(tmpu,tmpv,obj.grid.lon_rho,obj.grid.lat_rho,tmpang);
					obj.data.(A.type).ws.data  = tmpws;
					obj.data.(A.type).wsc.data = tmpwsc;
					obj.data.(A.type).ws.data  = obj.data.(A.type).ws.data .* obj.grid.mask_rho;
					obj.data.(A.type).wsc.data = obj.data.(A.type).wsc.data .* obj.grid.mask_rho;
					obj.data.(A.type).ws.name  = 'averaged wind-stress';
					obj.data.(A.type).wsc.name = 'averaged wind-stress curl';
					obj.data.(A.type).ws.units = '$N$ $m^{-2}$';
					obj.data.(A.type).wsc.units = '$N$ $m^{-2}$';
				% Mixed-layer depth (MLD)
				elseif strcmp(vars{i},'MLD') | strcmp(vars{i},'mld');
					obj  = zslice(obj,{'sigma0'},10,file,'type',A.type,'clear','off');
					dat  = obj.data.(A.type).sigma0.data;
					[nx,ny,nz,nt] = size(dat);
					N    = nx*ny*nt;
					d10  = squeeze(obj.data.(A.type).sigma0.slice);
					d10  = d10 + 0.03; % new threshold (10m density + 0.03 kg/m3)
					d10  = permute(repmat(d10,[1 1 1 nz]),[1 2 4 3]);
					dep  = obj.grid.(A.type).z_r;
					dat  = reshape(dat,[nx*ny*nt nz]);
					d10  = reshape(d10,[nx*ny*nt nz]);
					dep  = reshape(dep,[nx*ny*nt nz]);
					lvls = zeros(size(dat));
					lvls(dat<d10)=lvls(dat<d10)+1;
					lvls  = squeeze(nansum(lvls,2));
					lvls  = nz - lvls;
					for i = 1:(nx*ny*nt)
						if ~lvls(i)==nz
							outdep(i) = dep(i,lvls(i));	
						else
							outdep(i) = NaN;
						end
					end
					outdep = reshape(outdep,[nx ny nt]);
					for t = 1:nt
						outdep(:,:,t) = outdep(:,:,t).*obj.grid.mask_rho;
					end
					% Save ROMS data (make positive)
					obj.data.(A.type).MLD.data  = outdep;	
					obj.data.(A.type).MLD.data  = obj.data.(A.type).MLD.data;
					obj.data.(A.type).MLD.name  = 'averaged mixed-layer depth';
					obj.data.(A.type).MLD.units = '$m$';
				% ChlA
				elseif strcmp(vars{i},'SFC_CHL') | strcmp(vars{i},'sfc_chl');
					obj     = zslice(obj,{'TOT_CHL'},obj.grid.z_avg_dep(obj.grid.z_avg_dep<50),file);
					tmpchla = obj.data.(A.type).TOT_CHL.slice;
					tmpchla = squeeze(nanmean(tmpchla,3));
					obj.data.(A.type).SFC_CHL.data  = tmpchla .* obj.grid.mask_rho;
					obj.data.(A.type).SFC_CHL.name  = 'averaged surface chlA';
					obj.data.(A.type).SFC_CHL.units = '$mg$ $chlA$ $m^{-3}$';
                % Okubo-Weiss, vorticity, sN, or sS (normal, shear of strain)
                elseif strcmp(vars{i},'OW') | strcmp(vars{i},'vort') | strcmp(vars{i},'sN') | strcmp(vars{i},'sS');
					if isempty(A.ip) & isempty(A.dep)
						disp('Supply a depth input via ''ip'' or ''dep''');
						return
					elseif ~isempty(A.ip);
						obj = ipslice(obj,{'u','v'},A.ip,file);
					elseif ~isempty(A.dep);
						if ismember(A.dep,obj.grid.z_avg_dep);
							disp('Choose a depth from obj.grid.z_avg_dep only');
							return
						end
						obj = zslice(obj,{'u','v'},A.dep,file);
					end	
					[OW,vort,sS,sN] = romsMaster.okubo_weiss(obj.data.(A.type).u.data,obj.data.(A.type).v.data,obj.grid.pm,obj.grid.pn);
					obj.data.(A.type).OW.slice   = OW;
					obj.data.(A.type).OW.name    = 'averaged Okubo-Weiss parameter';
					obj.data.(A.type).OW.units   = '$s^{-2}$';
					obj.data.(A.type).vort.slice = vort;
					obj.data.(A.type).vort.name  = 'averaged relative vorticity';
					obj.data.(A.type).vort.units = '$s^{-1}$';
					obj.data.(A.type).sN.slice   = sN;
					obj.data.(A.type).sN.name    = 'averaged normal component of strain';
					obj.data.(A.type).sN.units   = '$s^{-1}$';
					obj.data.(A.type).sS.slice   = sS;
					obj.data.(A.type).sS.name    = 'averaged shear component of strain';
					obj.data.(A.type).sS.units   = '$s^{-1}$';
				% Jden_N2O
				elseif strcmp(vars{i},'Jden_N2O');
					obj = loadData(obj,{'DENITRIF2','N2OSODEN_CONS'},file,'type',A.type);
					obj.data.(A.type).Jden_N2O.data = (obj.data.(A.type).DENITRIF2.data ./ 2) - obj.data.(A.type).N2OSODEN_CONS.data; 
					obj.data.(A.type).Jden_N2O.name = 'Net $N_2O$ production from denitrification';
					obj.data.(A.type).Jden_N2O.units = '$mmol$ $N$ $m^{-3}$ $s^{-1}$';
				% Jnit_N2O
				elseif strcmp(vars{i},'Jnit_N2O');
					obj = loadData(obj,{'N2OAMMOX'},file,'type',A.type);
					obj.data.(A.type).Jnit_N2O.data = obj.data.(A.type).N2OAMMOX.data - obj.data.(A.type).N2OAO1_CONS.data;	
					obj.data.(A.type).Jnit_N2O.name = 'Net $N_2O$ production from nitrification';	
					obj.data.(A.type).Jnit_N2O.units = '$mmol$ $N$ $m^{-3}$ $s^{-1}$';
				% Apparent Oxygen Utiliziation (AOU)
				elseif strcmp(vars{i},'AOU');
					obj = loadData(obj,{'temp','salt','O2'},file,'type',A.type);
					o2_sat = romsMaster.o2_sat(obj.data.(A.type).temp.data,obj.data.(A.type).salt.data);
					obj.data.(A.type).AOU.data  = o2_sat - obj.data.(A.type).O2.data;
					obj.data.(A.type).AOU.name  = 'averaged AOU';
					obj.data.(A.type).AOU.units = '$mmol$ $O_2$ $m^{-3}$';
				% Sat/Delta N2O
				elseif strcmp(vars{i},'DeltaN2O');
					obj = loadData(obj,{'temp','salt','N2O'},file,'type',A.type);
					n2o_sat = romsMaster.n2o_sat(obj.data.(A.type).temp.data,obj.data.(A.type).salt.data);
					obj.data.(A.type).satN2O.data    = n2o_sat;
					obj.data.(A.type).satN2O.name    = 'averaged $N_2O_{sat}$';
					obj.data.(A.type).satN2O.units   = '$mmol$ $N_2O$ $m^{-3}$';
					obj.data.(A.type).deltaN2O.data  = obj.data.(A.type).N2O.data - n2o_sat;
					obj.data.(A.type).deltaN2O.name  = 'averaged $\Delta N_2O$';
					obj.data.(A.type).deltaN2O.units = '$mmol$ $N_2O$ $m^{-3}$';
				% OMZ thickness
				elseif strcmp(vars{i},'OMZ');
					obj = loadDepth(obj,file);
					obj = loadData(obj,{'O2'},file,'type',A.type);
					tmp.OMZ = nan(obj.grid.nx,obj.grid.ny,length(A.thresh),obj.grid.nt);
					for th = 1:length(A.thresh)
						tmp.O2 = obj.data.(A.type).O2.data;	
						tmp.Hz = obj.grid.(A.type).Hz;
						tmp.Hz(tmp.O2>A.thresh(th)) = 0;
						tmp.OMZ(:,:,th,:) = squeeze(nansum(tmp.Hz,3));
					end
					obj.data.(A.type).OMZ.int    = tmp.OMZ;
					for t = 1:obj.info.avg.time(file);
						for th = 1:length(A.thresh)
							obj.data.(A.type).OMZ.tot(th,t) = ...
								nansum(obj.data.(A.type).OMZ.int(:,:,th,t).*obj.grid.mask_rho.* obj.grid.area_rho,'all');
						end
					end
					obj.data.(A.type).OMZ.name   = 'OMZ thickness';
					obj.data.(A.type).OMZ.units  = '$m$'; 
					obj.data.(A.type).OMZ.thresh = A.thresh;
				% NH4 / NO2 ratio
				elseif strcmp(vars{i},'NH4vNO2') | strcmp(vars{i},'NO2vNH4');
					obj = loadData(obj,{'NH4','NO2'},file,'type',A.type);
					tmpnh4 = obj.data.(A.type).NH4.data;
					tmpnh4(tmpnh4<0) = 0;
					tmpno2 = obj.data.(A.type).NO2.data;
					tmpno2(tmpno2<0) = 0;
					tmpratio = tmpnh4 ./ tmpno2;
					tmpratio(tmpratio==0)   = NaN;
					tmpratio(tmpratio==Inf) = NaN;
					obj.data.(A.type).NH4vNO2.data   = tmpratio;
					obj.data.(A.type).NH4vNO2.name   = '$NH_4$ / $NO_2$';
                    obj.data.(A.type).NH4vNO2.units  = 'ratio ($\mu M$/$\mu M$)';
					tmpratio = tmpno2 ./ tmpnh4;
					tmpratio(tmpratio==0)   = NaN;
					tmpratio(tmpratio==Inf) = NaN;
					obj.data.(A.type).NO2vNH4.data   = tmpratio;
					obj.data.(A.type).NO2vNH4.name   = '$NO_2$ / $NH_4$';
                    obj.data.(A.type).NO2vNH4.units  = 'ratio ($\mu M$/$\mu M$)';
				elseif strcmp(vars{i},'NOX')
					obj = loadData(obj,{'NO3','NO2'},file,'type',A.type);
					tmpno3 = obj.data.(A.type).NO3.data;
					tmpno2 = obj.data.(A.type).NO2.data;
					tmpno3(tmpno3<0) = 0;
					tmpno2(tmpno2<0) = 0;
					tmpnox = tmpno3 + tmpno2;
					obj.data.(A.type).NOX.data  = tmpnox;
					obj.data.(A.type).NOX.name  = 'averaged $NO_x$';
					obj.data.(A.type).NOX.units = obj.data.(A.type).NO3.units;
				else
					disp([vars{i},' is already, or cant be, calculated']);
				end
			end
		end % end method computeVar
	
		%--------------------------------------------------------------------------------
		function obj = ipslice(obj,vars,ip,file,varargin)
			% ------------------
			% Slices ROMS along a given isopycnal
			% 
			% Usage:
			% - obj = ipslice(obj,vars,ip,file);
			% 
			% Inputs:
			% - vars = cell array of ROMS variables to slice
			% - ip   = isopycnal(s) to slice along
			% - file = file to slice
			%
			% Optional:
			% - type = 'avg','his','rst',etc
			%
			% Example:
			% - obj = ipslice(obj,{'O2'},26.5,file);
			% -------------------	

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

			% Compute sigma0
			try; obj.data.(A.type).sigma0.data;
			catch; obj = computeVar(obj,{'sigma0'},file,'type',A.type);
			end

			% Load variables
			for i = 1:length(vars)
				if ismember(vars{i},obj.info.(A.type).var2d) | ismember(vars{i},obj.info.(A.type).var3d);
					obj = loadData(obj,vars(i),file,'type',A.type);
				else
					obj = computeVar(obj,vars(i),file,'type',A.type);
				end
			end

			% Interpolate variable to ip level
			for i = 1:length(vars)
				vnew   = nan(obj.grid.nx,obj.grid.ny,length(ip),sum(obj.info.(A.type).time(file)));
				tmpvar = obj.data.(A.type).(vars{i}).data;
				tmpip  = obj.data.(A.type).sigma0.data;
				for z = 1:length(ip);
                    % Get sizes
                    [lx,ly,lz,lt] = size(tmpvar);
                    N = lx*ly*lt;

                    % Convert to 2D
                    tmp.ip = permute(tmpip,[1 2 4 3]);
                    tmp.var = permute(tmpvar,[1 2 4 3]);
                    tmp.ip = reshape(tmp.up,[lx*ly*lt lz]);
                    tmp.var = reshape(tmp.var,[lx*ly*lt lz]);

                    % Find level where z_r > z
                    a=tmp.ip>ip(z);
                    levs=squeeze(sum(a,2));
                    levs(levs==lz)=lz;
                    levs(levs==0)=1;
                    mask=levs./levs;
                    mask(isnan(mask)) = 0;
                    mask(levs==1) = 0;
                    mask(levs==lz) = 0;

                    % Index
                    zidx1  = [1:N]'+((levs(1:N)-1).*N);
                    zidx2  = [1:N]'+((levs(1:N)).*N);
                    z1 = reshape(tmp.ip(zidx1),[lx,ly,lt]);
                    z2 = reshape(tmp.ip(zidx2),[lx,ly,lt]);
                    v1 = reshape(tmp.var(zidx1),[lx,ly,lt]);
                    v2 = reshape(tmp.var(zidx2),[lx,ly,lt]);
                    z1(mask==0) = NaN;
                    z2(mask==0) = NaN;
                    v1(mask==0) = NaN;
                    v2(mask==0) = NaN;

                    % Perform linear interpolation
                    vnew(:,:,z,:) = [v2.*(z1-zdep(z)) + v1.*(zdep(z)-z2)]./[(z1-z2)];
				end
				% Replace data
				obj.data.(A.type).(vars{i}).data  = [];
				obj.data.(A.type).(vars{i}).slice = squeeze(vnew);
			end
            % Save into grid
            tmpip = repmat(ip,[1 lx ly t]);
            obj.slice.sigma0 = permute(tmpip,[2 3 1 4]);
            obj.slice.ip     = ip;
		end % end method ipslice

		%--------------------------------------------------------------------------------
		function obj = zslice(obj,vars,zdep,file,varargin)
			% ------------------
			% Slices ROMS along a given depth 
			% Use obj.grid.z_avg_dep for WOA18 depths
			% 
			% Usage:
			% - obj = zslice(obj,vars,zdep,file);
			% 
			% Inputs:
			% - vars = cell array of ROMS variables to slice
			% - zdep = depth(s) to slice along
			% - file = file to slice
			%
			% Optional:
			% - type  = 'avg','his','rst',etc
			% - clear = 'on' to clear raw data after slicing ('off' otherwise) 
			%
			% Example:
			% - obj = zslice(obj,{'O2'},[0 100 200],file);
			% -------------------	
            
			% Process optional inputs
            A.type = 'avg';
			A.clear = 'on';
            A = parse_pv_pairs(A,varargin);

			% Compute depth 
			try; obj.grid.(A.type).z_r;
			catch; obj = loadDepth(obj,file,'type',A.type);
			end

			% Change '0' to 5
			zdep(zdep==0) = 5;

			% Load variables
			for i = 1:length(vars)
				if ismember(vars{i},obj.info.(A.type).var2d) | ismember(vars{i},obj.info.(A.type).var3d);
					obj = loadData(obj,vars(i),file,'type',A.type);
				else
					obj = computeVar(obj,vars(i),file,'type',A.type);
				end
			end

			% Interpolate variable to zdep level
			for i = 1:length(vars)
				vnew   = nan(obj.grid.nx,obj.grid.ny,length(zdep),sum(obj.info.(A.type).time(file)));
				tmpvar = obj.data.(A.type).(vars{i}).data;
				if strcmp(vars{i},'u');
					tmpz = -obj.grid.(A.type).z_u;
					vnew = vnew(1:end-1,:,:,:);
				elseif strcmp(vars{i},'v');
					tmpz = -obj.grid.(A.type).z_v;
					vnew = vnew(:,1:end-1,:,:);
				else
					tmpz   = -obj.grid.(A.type).z_r;
				end
				for z = 1:length(zdep)		
					% Get sizes
					[lx,ly,lz,lt] = size(tmpvar);
					N = lx*ly*lt;
		
					% Convert to 2D				
					tmp.z = permute(tmpz,[1 2 4 3]);
					tmp.var = permute(tmpvar,[1 2 4 3]);
					tmp.z = reshape(tmp.z,[lx*ly*lt lz]);
					tmp.var = reshape(tmp.var,[lx*ly*lt lz]);

					% Find level where z_r > z
					a=tmp.z>zdep(z);
					levs=squeeze(sum(a,2));
					levs(levs==lz)=lz;
					levs(levs==0)=1;
					mask=levs./levs;
					mask(isnan(mask)) = 0;
					mask(levs==1) = 0;
					mask(levs==lz) = 0;

					% Index
					zidx1  = [1:N]'+((levs(1:N)-1).*N);
					zidx2  = [1:N]'+((levs(1:N)).*N);
					z1 = reshape(tmp.z(zidx1),[lx,ly,lt]);
					z2 = reshape(tmp.z(zidx2),[lx,ly,lt]);
					v1 = reshape(tmp.var(zidx1),[lx,ly,lt]);
					v2 = reshape(tmp.var(zidx2),[lx,ly,lt]);
					z1(mask==0) = NaN;
					z2(mask==0) = NaN;
					v1(mask==0) = NaN;
					v2(mask==0) = NaN;

                    % Perform linear interpolation
					try
                    vnew(:,:,z,:) = [v2.*(z1-zdep(z)) + v1.*(zdep(z)-z2)]./[(z1-z2)];
					catch
						disp('debug zslice');
						keyboard
					end
				end
                % Replace data
				if strcmp(A.clear,'on');
					obj.data.(A.type).(vars{i}).data   = [];
				end
                obj.data.(A.type).(vars{i}).slice  = squeeze(vnew);
            end
            % Save into grid
            tmpdepth = repmat(zdep,[1 lx ly lt]);
            obj.slice.depth = permute(tmpdepth,[2 3 1 4]);
            obj.slice.zdep  = zdep;
		end % end method zslice

		%--------------------------------------------------------------------------------
		function obj = getProfile(obj,vars,lon,lat,file,varargin)
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
			%   - file  = input file number
			%
			% Inputs (varargin):
			%   - type     = 'avg' or 'z_avg' (raw cant be used with yr_range)
			%	- yr_range = range of years to load and average (i.e. 2045:2049)
			%
			% Example:
			%	- obj = getProfile(obj,{'temp','salt','O2','NO3'},[250 250],[-15 -20]);
			% -------------------

			% process inputs
			A.type     = 'avg';
			A          = parse_pv_pairs(A,varargin);	

			% Get indices of nearest point
			lon_idx = [];
			lat_idx = [];
			for i = 1:length(lon);
				% Get index of nearest lat/lon grid cell
				all_dist                = distance(lat(i),lon(i),obj.grid.lat_rho,obj.grid.lon_rho);
				[lon_idx(i),lat_idx(i)] = find(all_dist == min(all_dist(:)));
			end
			   
			% Load data
			for i = 1:length(vars)
				try
					obj = loadData(obj,vars(i),file,'type',A.type);
				catch
					obj = computeVar(obj,vars(i),file,'type',A.type);
				end
			end
				
			% Go through all variables and lon/lats
			for v = 1:length(vars)
				disp(['...grabbing ',vars{v},' profile(s)']);
				tmp.data = obj.data.(A.type).(vars{v}).data;	
				% Save profile data
				for i = 1:length(lon);
					obj.data.(A.type).(vars{v}).profile(i,:,:)  = squeeze(tmp.data(lon_idx(i),lat_idx(i),:,:));
					if strcmp(A.type,'z_avg');
						obj.profile.depth	= repmat(obj.grid.z_avg_dep,[1 length(lon)]);; 
					else
						obj.profile.depth(i,:,:) = squeeze(obj.grid.(A.type).z_r(lon_idx(i),lat_idx(i),:,:));
					end
					obj.profile.lon(i)  = obj.grid.lon_rho(lon_idx(i),lat_idx(i));
					obj.profile.lat(i)  = obj.grid.lat_rho(lon_idx(i),lat_idx(i));
					obj.profile.xidx(i) = lon_idx(i);
					obj.profile.yidx(i) = lat_idx(i);
				end
			end
		end % end method getProfile

		%--------------------------------------------------------------------------------
		function obj = loadDiag(obj,vars,depth,varargin)
			% ------------------
			% Loads and interpolates monthly 2D validation data to the ROMS grid
			% See obj.paths.diag for available fields
			%
			% Usage:
			% - obj = loadDiag(obj,vars,depth,varargin);
			%
			% Inputs:
			% - vars  = validation variable(s) to load, as a cell array; paths set in getDiagpaths 
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

			% Check inputs
			diagfields = fieldnames(obj.paths.diag); 
			for i = 1:length(vars)
				if ~strcmp(vars{i},diagfields) & ~strcmp(vars{i},upper(diagfields));
					disp(' ');
					disp([vars{i},' is not a diagnostic variable']);
					disp(['Filling with NaN']);
					obj.diag.(vars{i}).slice = nan(obj.grid.nx,obj.grid.ny,length(depth),12);
					obj.diag.(vars{i}).name = 'null';
					obj.diag.(vars{i}).units = 'null';
					obj.diag.(vars{i}).depth = depth;
					disp(' ');
					skip(i) = 1;
				else
					skip(i) = 0;
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
			tmp.outlon = obj.grid.lon_rho;
			tmp.outlat = obj.grid.lat_rho;

			% Process each variable
			for i = 1:length(vars)		
				% Skip empty data
				if skip(i) == 1
					continue
				end
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
						if ndims(tmp.lon)==3
							tmp.lat = squeeze(tmp.lat(:,:,1));
							tmp.lon = squeeze(tmp.lon(:,:,1));
						end
					end

					% Make lon/lat data gridded if it isnt setup that way
					if sum(size(tmp.lon(:,:))==1) > 0
						[tmp.lat,tmp.lon] = meshgrid(tmp.lat, tmp.lon);
					end;
					tmp = romsMaster.struct2double(tmp);
				
					% Get indeces for reduced domain interpolation
					% Note: this breaks if ROMS boundary longitude is close to 0 or 360
					idx = find(tmp.lon(:) > obj.grid.minlon_rho-A.outer ...
							 & tmp.lon(:) < obj.grid.maxlon_rho+A.outer ...
							 & tmp.lat(:) > obj.grid.minlat_rho-A.outer ...
							 & tmp.lat(:) < obj.grid.maxlat_rho+A.outer);

					% Initialize interpolated output
					obj.diag.(vars{i})(j).slice = nan(obj.grid.nx,obj.grid.ny,length(depth),12);

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
				
						% Fix tmpdepth?
						if ndims(tmpdepth)==3
							tmpdepth = squeeze(tmpdepth(1,1,:));
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
							tmpout{z,k}(isnan(obj.grid.mask_rho)) = nan;
						end % end z-loop
					end % end k-loop
					% Save interpolated data
					for ll = 1:length(depth);
						obj.diag.(vars{i})(j).slice(:,:,ll,:)  = cat(3, tmpout{ll,:});
						obj.diag.(vars{i})(j).slice(:,:,ll,:)  = single(obj.diag.(vars{i})(j).slice(:,:,ll,:));
					end
					obj.diag.(vars{i})(j).slice  = squeeze(obj.diag.(vars{i})(j).slice);
					obj.diag.(vars{i})(j).depth = depth;
					obj.diag.(vars{i})(j).units = curVar.units{j};
					obj.diag.(vars{i})(j).name  = curVar.name{j};
				end % end j-loop
			end % end i-loop
		end % end method loadDiag

		%--------------------------------------------------------------------------------
		function obj = make_z_avg(obj)
			% ------------------
			% - Vertically interpolates the sigma coordinates ROMS file to 
			% - WOA constant depth levels
			%
			% - Usage:
			% - obj = make_z_avg(obj,avgfile)
			%
			% - Inputs:
			% - avgfile - file to perform depth slices on 
			%
			% - Example:
			% - obj = make_z_avg(obj,'avg_2009.nc');
			% ------------------	
			
			% Get initial directory
			od = pwd;

			% Set avgfile name
			[pathstr, avgname, avgext] = fileparts(obj.paths.monthly);
			avgfile = [avgname,avgext];

			% Create z_avg_dir if it doesnt exist
			z_avg_dir = [obj.paths.runPath,'z_avg'];
			if exist(z_avg_dir) ~= 7
				mkdir(z_avg_dir);
			end
			
			% Create z_avg file (woa depths) 
			path1 = [obj.paths.runPath,'avg'];
			path2 = [obj.paths.runPath,'z_avg'];
			cd(path1)
			cmd = ['ln -s ',obj.paths.grid,' .'];
			system(cmd);
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

			% Return to original directory
			cd(od);
		end % end method make_z_avg

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
			%   - caxis       colorbar limits 
			%   - log         log-scale (1), use with caxis to set limits
			% -----------------------
			
			% User-inputs
			A.lonbounds  = [];
			A.latbounds  = [];
			A.lonticks   = [];
			A.latticks   = [];
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
			A.caxis      = [];
			A.log        = 0;
			A = parse_pv_pairs(A,varargin);

			% Balance override
			if A.bal == 1
				A.cmap = cmocean('-balance');
			end

			% Make double
			lon = double(obj.grid.lon_rho);
			lat = double(obj.grid.lat_rho);

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
			if isempty(A.latticks);
				latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
			else
				latticks = A.latticks;
			end
			if isempty(A.lonticks);
				lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));
			else	
				lonticks = A.lonticks;
			end

			% Initiate figure
			fig = piofigs(A.figtype,A.figdim);

			% Get colormap limits
			if isempty(A.levels)
				clims = romsMaster.prclims(dat,'prc',A.prc,'bal',A.bal);
				clevs = linspace(clims(1),clims(2),31); 
			else
				clevs = A.levels;
				clims = [A.levels(1) A.levels(end)];
			end
			if ~isempty(A.caxis);
				clevs(1) = A.caxis(1);
				clevs(end) = A.caxis(2);
				clims = A.caxis;
			end
			if max(dat(:)) == 0 & min(dat(:)) == 0 | isnan(max(dat)) ==1 & isnan(min(dat(:)) == 1);
				dat = nan(size(dat));
				clevs = [0 1];
				clims = linspace(0,1,11);
			else
				dat(dat<clevs(1))   = clevs(1);
				dat(dat>clevs(end)) = clevs(end);
			end

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
			try
				m_contourf(lon,lat,dat,clevs,'LineStyle','none');
			catch
			end
			cb = colorbar; drawnow
			cb.FontSize = A.fontsize;
			try
				caxis([clims])
			catch
				caxis([0 1]);
			end
			ax = get(gca);
			m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
			colormap(gca,A.cmap);
			if A.log == 1 & ~isempty(A.caxis)
				set(gca,'ColorScale','log');
				caxis(A.caxis);
			end
			
			% Print figure
			if nargout<1
				pltjpg(1);
			end
		end % end method mapPlot

		%--------------------------------------------------------------------------------
		function [fig,ax,cb] = axPlot(obj,fig,ax,dat,varargin);
			% -----------------------
			% - A way to quickly plot a 2D field into a specific axis (i.e. subplot)
			%
			% - Usage:
			%	- [fig,ax,cb] = axPlot(obj,fig,ax,dat,varargin);
			%
			% - Inputs:
			%	- fig = figure handle
			%	- ax  = axes handle
			%	- dat = 2D field to plot (prepare it before using the script)
			%
			% - Varargin:
			%	- lonbounds:  x-boundaries (defaults to whole domain)
			%	- latbounds:  y-boundaries (defaults to whole domain)
			%	- ticks:      2 = fancy, 1 = on, 0 = off
			%   - xtick/ytick 0 = off (to toggle only one tick on/off)
			%	- background: background color (default 'LightGray');
			%	- coastcolor: coast color (default 'DimGray');
			%	- fontsize:   default 10
			%	- prc:        percentage to limit colorbar axes
			%	- bal:        balance colorbar around 0 (1 == yes, 0 == no)
			%	- levels:     hard-coded levels to plot (can't be used with A.prc)
			%	- cmap:       colormap(default = thermal for no balance, balance for balance)
			%	- caxis       colorbar limits 
			%	- log         log-scale (1), use with caxis to set limits
			% -----------------------
			
			% User-inputs
			A.lonbounds  = [];
			A.latbounds  = [];
			A.lonticks   = [];
			A.latticks   = [];
			A.ticks      = 0;
			A.ytick      = 1;
			A.xtick      = 1;
			A.background = rgb('LightGray');
			A.coastcolor = rgb('DimGray');
			A.fontsize   = 10;
			A.prc        = 2;
			A.bal        = 0;
			A.levels     = [];
			A.cmap       = cmocean('thermal');
			A.caxis      = [];
			A.log        = 0;
			A = parse_pv_pairs(A,varargin);

			% Balance override
			if A.bal == 1
				A.cmap = cmocean('-balance');
			end

			% Make double
			lon = double(obj.grid.lon_rho);
			lat = double(obj.grid.lat_rho);

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
			if isempty(A.latticks);
				latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
			else
				latticks = A.latticks;
			end
			if isempty(A.lonticks);
				lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));
			else	
				lonticks = A.lonticks;
			end

			% Set current axes
			set(0,'CurrentFigure',fig);
			set(fig,'CurrentAxes',ax);

			% Get colormap limits
			if isempty(A.levels)
				clims = romsMaster.prclims(dat,'prc',A.prc,'bal',A.bal);
				clevs = linspace(clims(1),clims(2),31); 
			else
				clevs = A.levels;
				clims = [A.levels(1) A.levels(end)];
			end
			if ~isempty(A.caxis);
				clevs(1) = A.caxis(1);
				clevs(end) = A.caxis(2);
				clims = A.caxis;
			end
			if max(dat(:)) == 0 & min(dat(:)) == 0 | isnan(max(dat)) ==1 & isnan(min(dat(:)) == 1);
				dat = nan(size(dat));
				clevs = [0 1];
				clims = linspace(0,1,11);
			else
				dat(dat<clevs(1))   = clevs(1);
				dat(dat>clevs(end)) = clevs(end);
			end

			% Make map
			m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
			hold on
			if A.ticks == 0
				m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
					   'xticklabels',[],'yticklabels',[],'backgroundcolor',rgb('LightGray')); drawnow
			elseif A.ticks == 1
				if A.ytick == 1 & A.xtick == 1
					m_grid('box','on','linestyle','none','xtick',lonticks,...
						   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
				elseif A.ytick ==1 & A.xtick == 0
					m_grid('box','on','linestyle','none','xtick',0,...
						   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',[]);
				elseif A.ytick == 0 & A.xtick == 1
					m_grid('box','on','linestyle','none','xtick',lonticks,...
						   'ytick',0,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',[],'xticklabels',lonticks);
				end
			elseif A.ticks == 2
				m_grid('box','fancy','linestyle','none','xtick',lonticks,...
					   'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
			end
			hold on
			try
				m_contourf(lon,lat,dat,clevs,'LineStyle','none');
			catch
			end
			cb = colorbar; drawnow
			cb.FontSize = A.fontsize;
			try
				caxis([clims])
			catch
				caxis([0 1]);
			end
			m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
			colormap(gca,A.cmap);
			if A.log == 1 & ~isempty(A.caxis)
				set(gca,'ColorScale','log');
				caxis(A.caxis);
			end
		end % end method axPlot

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
			A.lonbounds  = [99 300];
			A.latbounds  = [-46 62];
			A.ticks      = 0;
			A.box        = 'on';
			A.background = rgb('LightGray');
			A.coastcolor = rgb('DimGray');
			A.fontsize   = 7;
			A.figtype    = 'mfig';
			A.figdim     = 1;
			A = parse_pv_pairs(A,varargin);

			% Make double
			lon = double(obj.grid.lon_rho);
			lat = double(obj.grid.lat_rho);

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

			% Make map
			fig = piofigs(A.figtype,A.figdim);
			set(0,'CurrentFigure',fig);
			m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
			hold on
			if A.ticks == 0
				m_grid('box',A.box,'linestyle','none','xtick',0,'ytick',0,...
					   'xticklabels',[],'yticklabels',[],'backgroundcolor',A.background); drawnow
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
			%   - units:      string containing units for colorbar
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
			A.prc        = 0.1;
			A.bal        = 0;
			A.levels     = [];
			A.difflevels = [];
			A.cmap       = 'thermal';
			A.units      = [];
			A = parse_pv_pairs(A,varargin);

			% Balance override
			if A.bal == 1
				A.cmap = cmocean('-balance');
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
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels,'cmap',cmocean('-balance',length(A.difflevels)));

			% Add units?
			if ~isempty(A.units)
				ylabel(cbs(1),A.units,'Interpreter','Latex','FontSize',A.fontsize);
				ylabel(cbs(2),A.units,'Interpreter','Latex','FontSize',A.fontsize);
				ylabel(cbs(3),A.units,'Interpreter','Latex','FontSize',A.fontsize);
			end

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
		function [fig,cb] = slicePlot(obj,slicedata,t,varargin)
			% ------------------
			% Plots 2D sliced variables obtained from sliceROMS or sliceDiag
			%
			% Usage:
			% - [fig,cb] = slicePlot(obj,vars,t,varargin);
			%
			% Inputs:
			% - slicedata = sliced data from sliceROMS or sliceDiag (must be size(obj.slice.deg)) 
			% - t         = time to plot (0 == average)
			%
			% Optional inputs (varargin):
			% - slice:      if several slices were made, specify which one you want to plot here
			% - xlims:      hard-coded x-limits (degrees)
			% - zlims:      hard-coded z-limits (meters)
			% - cmap:       colormaps (if used, must be a cell array of length (vars))
			% - fontsize:   default 10
			% - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
			% - figdim:     default 1 (same height as width, its a multiplier)
			% - prc:        percentage to limit colorbar axes
			% - bal:        balance colorbar around 0 (1 == yes, 0 == no)
			% - levels:     hard-coded levels to plot (can't be used with A.prc)
			% - background: background color (i.e. coast). Default = rgb('DimGray');
			%
			% Example:
			% - [fig,cb] = slicePlot(obj,{'O2'},0,'xlims',[140 260],'zlims',[-500 0]);
			% This will averaged O2 from 0 -500m with lon limits of 140/260
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
			A.background = rgb('DimGray');
			A = parse_pv_pairs(A,varargin);

			% Check for slice
			if isempty(obj.slice)
				disp('Slice your data first using sliceROMS or sliceDiag');
				return
			end
		
			% Check that A.slice has also been reduced
			if ndims(obj.slice.depth)>3  & isempty(A.slice) 
				disp('You forgot to specify which slice! Killing');
				return
			elseif ndims(obj.slice.depth)>3
				if t == 0
					tmpdeg = nanmean(squeeze(obj.slice.deg(:,:,A.slice,:)),3);
					tmpdep = nanmean(squeeze(obj.slice.depth(:,:,A.slice,:)),3);
				else
					tmpdeg = squeeze(obj.slice.deg(:,:,A.slice,t));
					tmpdep = squeeze(obj.slice.depth(:,:,A.slice,t));
				end
			elseif ndims(obj.slice.depth)==3
				tmpdeg = squeeze(obj.slice.deg(:,:,A.slice));
				tmpdep = squeeze(obj.slice.depth(:,:,A.slice));
			elseif ndims(obj.slice.depth)==2
				tmpdeg = obj.slice.deg;
				tmpdep = obj.slice.depth;
			end

			% Check for xlim or zlims
			if ~isempty(A.xlims);
				xind = A.xlims(1) <= tmpdeg & tmpdeg <= A.xlims(2);
			end
			if ~isempty(A.zlims);
                zind = A.zlims(1) <= tmpdep & tmpdep <= A.zlims(2);
            end

			% Grab data, reduce to spatial and data limits	
			if t == 0
				slicedata = nanmean(squeeze(slicedata(:,:,A.slice,:)),3);
			else
				slicedata = squeeze(slicedata(:,:,A.slice,t));
			end
			
			% Get universal levels
			if isempty(A.levels)
				tmpdata = slicedata;
				% Reduce data?
				if ~isempty(A.xlims);
					tmpdata(xind==0) = NaN;
				end
				if ~isempty(A.zlims);
					tmpdata(zind==0) = NaN;
				end
				A.levels = romsMaster.prclims(tmpdata(:),'prc',A.prc,'bal',A.bal);
				A.levels = linspace(A.levels(1),A.levels(2),20);
			end
			slicedata(slicedata<A.levels(1)) = A.levels(1);
			slicedata(slicedata>A.levels(end)) = A.levels(end);
			
			% Generate figure(s)	
			fig = piofigs(A.figtype,A.figdim);
			set(0,'CurrentFigure',fig);
			if strcmp(obj.slice.coord,'latitude') | strcmp(obj.slice.coord,'longitude');
				contourf(tmpdeg,tmpdep,slicedata,A.levels,'linestyle','none');
			end
			if ~isempty(A.xlims);
				xlim(A.xlims);
			end
			if ~isempty(A.zlims);
				ylim(A.zlims);
			end
			if strcmp(obj.slice.coord,'longitude')
				xlabel('Latitude','Interpreter','Latex','FontSize',A.fontsize);
			    xlbl = get(gca,'XTickLabel');
				for i = 1:length(xlbl)
					if str2num(xlbl{i})<0
						newlbl{i} = [num2str(-str2num(xlbl{i})),char(176),'S'];
					else
						newlbl{i} = [num2str(str2num(xlbl{i})),char(176),'N'];
					end
				end
				set(gca,'XTickLabel',newlbl);	
			elseif strcmp(obj.slice.coord,'latitude')	
				xlabel('Longitude','Interpreter','Latex','FontSize',A.fontsize);
                xlbl = get(gca,'XTickLabel');
                for i = 1:length(xlbl)
					if str2num(xlbl{i})>180
                        newlbl{i} = [num2str(str2num(xlbl{i})-360),char(176),'W'];
                    else
                        newlbl{i} = [num2str(str2num(xlbl{i})),char(176),'E'];
					end
                end
				set(gca,'XTickLabel',newlbl);	
			end
			ylabel('Depth (m)','Interpreter','Latex','FontSize',A.fontsize);
			cb = colorbar('location','eastoutside');
			caxis([A.levels(1) A.levels(end)]);
			if ~isempty(A.cmap)
				if ischar(A.cmap)
					set(gca,'Colormap',cmocean(A.cmap,length(A.levels)));
				else
					set(gca,'Colormap',A.cmap);
				end
			end
			set(gca,'FontSize',A.fontsize);
			if nanmean(tmpdep(:)) > 0
				set(gca,'YDir','Reverse');
			end
			set(gca,'Color',A.background);
			set(gcf,'inverthardcopy','off');

			% Print if no output provided
			if nargout<1
				pltjpg(1);
			end
		end % end method slicePlot

		%--------------------------------------------------------------------------------
		function [figs,cbs] = sliceCmp(obj,dat1,dat2,t,varargin)
			% ------------------
			% A way to quickly plot slice comparisons
			% Produces 3 plots: dat1 and dat2 fields with the same colorbar, 
			% and a 3rd plot of the difference between them (dat1 - dat2)
			%
			% Usage:
			% - [figs,cbs] = sliceCmp(obj,dat1,dat2,t,varargin);
			%
			% Inputs:
			% - dat1: 2D field to plot (prepare it before using the script)
			% - dat2: same same but different 2D field to plot (prepare it before using the script)
			% - t:    time dimension to plot (0 == average)    
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
			% - [figs,cbs] = sliceCmp(obj,dat1,dat2,0,'xlims',[140 260],'zlims',[-500 0]);
			% This will take the average of dat1 + dat2 and plot them as well as their differences
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

			% Set levels if not supplied
			if isempty(A.levels);
				alldat  = [dat1(:) dat2(:)];
				lvls     = romsMaster.prclims(alldat,'prc',A.prc,'bal',A.bal);
				A.levels = linspace(lvls(1),lvls(2),20);
			end

			% Make figs(1) and figs(2)
			[figs(1),cbs(1)] = slicePlot(obj,dat1,t,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);
			[figs(2),cbs(2)] = slicePlot(obj,dat2,t,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);

			% Get differences
			diff_dat  = dat1 - dat2;
			if isempty(A.difflevels)
				A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1);
				A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
			end

			% Make figs(3), difference
			[figs(3),cbs(3)] = slicePlot(obj,diff_dat,t,...
				'slice',A.slice,'xlims',A.xlims,'zlims',A.zlims,'cmap','-balance','fontsize',A.fontsize,...
				'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels);

		end % end method sliceCmp

		%--------------------------------------------------------------------------------
		function [obj] = OMZthick(obj,omzthresh,varargin)
			% ------------------
			% Calculates OMZ thickness diagnostic oxygen products:
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
			% This will calculate OMZ (defined as O2 < 10 uMol) thickness from 
			% validation productions (WOA-18, Bianchi objective mapping 2)
			%
			% ------------------

           % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

			% Set defaults
			if nargin < 3
				diag = 0;
			end

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
			[in,on] = inpolygon(woa.LON,woa.LAT,obj.grid.polygon(:,1),obj.grid.polygon(:,2));
			woaind = find(in == 1 | on == 1);
			[in,on] = inpolygon(om2.LON,om2.LAT,obj.grid.polygon(:,1),obj.grid.polygon(:,2));
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
				tmpwoa  = romsMaster.grid_to_grid(woa.LON,woa.LAT,omz.woa(:,:,i),obj.grid.lon_rho,obj.grid.lat_rho);
				tmp.woa(:,:,i)  = tmpwoa .* obj.grid.mask_rho;
				tmpom2  = romsMaster.grid_to_grid(om2.LON,om2.LAT,omz.om2(:,:,i),obj.grid.lon_rho,obj.grid.lat_rho);
				tmp.om2(:,:,i)  = tmpom2 .* obj.grid.mask_rho;
			end
			omz.woa  = tmp.woa;
			omz.om2  = tmp.om2;

			% Save output
			obj.diag.OMZ(1).int    = omz.woa;
			obj.diag.OMZ(1).name   = 'OMZ Thickness (WOA-18)';
			obj.diag.OMZ(1).units  = '$m$';
			obj.diag.OMZ(1).thresh = omzthresh;
			obj.diag.OMZ(2).int    = omz.om2;
			obj.diag.OMZ(2).name   = 'OMZ Thickness (Bianchi 2012)';
			obj.diag.OMZ(2).units  = '$m$';
			obj.diag.OMZ(2).thresh = omzthresh;
		end % end method OMZthick

		%--------------------------------------------------------------------------------
		function [obj] = getAvgData(obj,vars,files,varargin)
			% -------------------
			% Get an average field across multiple files
			% 
			% Usage:
			% - obj = getAvgData(obj,vars,files,varargin)
			%
			% Inputs:
			% - vars     = cell array of variables to load and average
			% - files    = file #s to average across
			%
			% Optional Inputs (varargin):
			% - zlvl     = zlevel(s) to grab, if zslicing
			% - iplvl    = iplevel(s) to grab, if ipslicing
			%
			% Example:
			% - obj = getAvgData(obj,{'O2','NO3'},2000:2009);
			% - obj = getAvgData(obj,{'O2','NO3'},2000:2009,'lvl',[300 500]);
			% -------------------

			% defaults for optional  arguments
			A.zlvl  = [];
			A.iplvl = [];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

			% Loop through years and load variables
			for i = 1:length(vars)
				for j = 1:length(files)
					if ~isempty(A.zlvl)
						obj = zslice(obj,vars(i),A.zlvl,files(j));
						if j == 1
							out = obj.data.avg.(vars{i}).slice;
						else
							out = out + obj.data.avg.(vars{i}).slice;
						end
					elseif ~isempty(A.iplvl)
						obj = ipslice(obj,vars(i),A.iplvl,files(j));
						if j == 1
							out = obj.data.avg.(vars{i}).slice;
						else
							out = out + obj.data.avg.(vars{i}).slice;
						end	
					else
						obj = loadDepth(obj,files(j));
						obj = loadData(obj,vars(i),files(j));
						if j == 1
							out = obj.data.avg.(vars{i}).data;
							z_r = obj.grid.avg.z_r;
							z_u = obj.grid.avg.z_u;
							z_v = obj.grid.avg.z_v;
							z_w = obj.grid.avg.z_w;
							Hz  = obj.grid.avg.Hz;
						else
							out = out + obj.data.avg.(vars{i}).data;
							z_r = z_r + obj.grid.avg.z_r;
							z_u = z_u + obj.grid.avg.z_u;
							z_v = z_v + obj.grid.avg.z_v;
							z_w = z_w + obj.grid.avg.z_w;
							Hz  = Hz  + obj.grid.avg.Hz;
						end
					end
				end
				if ~isempty(A.zlvl) | ~isempty(A.iplvl)
					obj.data.avg.(vars{i}).slice = out./length(files);
				else
					obj.data.avg.(vars{i}).data = out./length(files);
					obj.grid.avg.z_r = z_r./length(files);
					obj.grid.avg.z_u = z_u./length(files);
					obj.grid.avg.z_v = z_v./length(files);
					obj.grid.avg.z_w = z_w./length(files);
					obj.grid.avg.Hz  = Hz./length(files);
				end
			end
		end % end method getAvgData	

		%--------------------------------------------------------------------------------
		function [obj] = getAvgSlice(obj,vars,choice,deg,files,varargin)
			% -------------------
			% Get an average depth slice across multiple files
			% 
			% Usage:
			% - obj = getAvgSlice(obj,vars,files,varargin)
			%
			% Inputs:
			% - vars     = cell array of variables to load and average
			% - choice   = 'lon','lat','xi','eta' 
            %              (lat/lon slices along a given lat/lon degree)
            %              (xi/eta slices along a given xi or eta index, use gridView)
			% - deg      = lon/lat degree(s) or x/y indicies
			% - files    = file #s to average across
			%
			% Example:
			% - obj = getAvgSlice(obj,{'O2','NO3'},'lon',0,1);
			% -------------------

			% defaults for optional  arguments
			A.type  = ['avg'];
			A.zlim  = [];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

			% Loop through years and load variables
			for i = 1:length(vars)
				for j = 1:length(files)
					obj = sliceROMS(obj,vars(i),choice,deg,files(j),'type',A.type,'zlim',A.zlim);
					if j == 1
						out = obj.data.avg.(vars{i}).slice;
						dep = obj.slice.depth;
					else
						out = out + obj.data.avg.(vars{i}).slice;
						dep = dep + obj.slice.depth;
					end
				end
				obj.data.avg.(vars{i}).data = out./length(files);
				obj.slice.depth = dep./length(files);
			end
		end % end method getAvgSlice	

		%--------------------------------------------------------------------------------
		function obj = equatorUcmp(obj,sect,file,varargin)
			% -----------------------
			% Compares equatorial current structure with Cravette2017 results or Johnson2002
			% Returns depth slices of u-velocity along specific longitude bands
			%
			% Usage:
			% - obj = equatorUcmp(obj,file,sect,varargin)
			%
			% Inputs:
			% - sect: Section to compare
			%		  choose from the following (as a string) for Cravatte 2017
			%		  '120W_90W', '135E_160E', '160E_180', '160W_120W', '179W_160W'
			%		  choose from the following (as an integer) for Johnson 2002
			%		  143, 156, 165, 180, 190, 205, 220, 235, 250, 265
			% - file: file(s) to load data from
			%
			% Inputs (optional):
			% - deps = depths to slice at, default = 0:10:500
			%
			% -----------------------

			% Process inputs
			A.deps = [0:10:500];
			
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
			% Grid 
			[depu,latu] = meshgrid(-depu,latu);
			
			% Get 3D u velocity data
			for i = 1:length(file)
				obj = zslice(obj,{'u'},A.deps,file(i));
				tmp{i} = obj.data.avg.u.slice;
			end
			if length(file)>1
				romsu = cat(4,tmp{i});
				romsu = nanmean(romsu,4); 
			else
				romsu  = nanmean(obj.data.avg.u.slice,4);
			end
			romsln = obj.grid.lon_u;
			romslt = obj.grid.lat_u;
			romsdp = A.deps;
			
			% Take slice of data for each time
			dmsn = obj.grid.ny;
			nz   = A.deps; nzz = length(nz);
			disp(' '); disp(['Slicing ROMS u data']);disp(' ');
			tmpslice = NaN(dmsn,nzz);
			tmpdata = squeeze(romsu);
			% - Interp to lon/lat line
			for i = 1:dmsn
				% - Fill tmpslice			
				for z = 1:nzz;
					try
						tmpslice(i,z) = interp1(squeeze(romsln(:,i)),squeeze(tmpdata(:,i,z)),sect_avg);
					catch
						disp('Cant interpolate ROMS at this longitude, try again');
						return
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
			tmpu     = tmpslice(:,:);
			idx      = find(isnan(tmpu)==0);
			F        = scatteredInterpolant(tmpdeg(idx),tmpdepth(idx),tmpu(idx));
			gridu{1} = F(latu,depu);

			% Save results in format for slicePlot 
			idx = find(strcmp(obj.info.avg.var3d,'u')==1);
			obj.data.avg.u.slice  = cat(3,gridu{:});
			obj.data.avg.u.name   = obj.info.avg.name3d{idx};
			obj.data.avg.u.units  = obj.info.avg.unit3d{idx};
			obj.diag.u.slice  = datu./100;
			obj.diag.u.name   = diagname;
			obj.diag.u.units  = obj.info.avg.unit3d{idx};
			obj.slice.deg   = latu;
			obj.slice.depth = depu;
			if nanmean(depu) < 0
				obj.slice.depth = -depu;
			end
			obj.slice.coord = 'longitude';
			obj.slice.sect  = sect_avg;
		end % end method equatorUcmp

        %--------------------------------------------------------------------------------
        function obj = zslice_data(obj,vars,zdep); 
            % -----------------------
			% Shortcut to zslice loaded data
			%
			% Inputs:
			% - obj  = romsObject
			% - vars = variable to slice 
			% - zdep = depth to slice
			%
			% Steps:
			%	Load z_r into obj.grid.avg.z_r;
			%   Load data into obj.grid.avg.(vars).data      
            % -----------------------

			% Interpolate variable to zdep level
			for i = 1:length(vars)
				tmpvar = obj.data.avg.(vars{i}).data;
				tmpz   = -obj.grid.avg.z_r;
				[lx,ly,lz,lt] = size(tmpvar);
				N = lx*ly*lt;
				vnew   = nan(lx,ly,length(zdep),lt);
				for z = 1:length(zdep)
					% Convert to 2D
					tmp.z = permute(tmpz,[1 2 4 3]);
					tmp.var = permute(tmpvar,[1 2 4 3]);
					tmp.z = reshape(tmp.z,[lx*ly*lt lz]);
					tmp.var = reshape(tmp.var,[lx*ly*lt lz]);

					% Find level where z_r > z
					a=tmp.z>zdep(z);
					levs=squeeze(sum(a,2));
					levs(levs==lz)=lz;
					levs(levs==0)=1;
					mask=levs./levs;
					mask(isnan(mask)) = 0;
					mask(levs==1) = 0;
					mask(levs==lz) = 0;

					% Index
					zidx1  = [1:N]'+((levs(1:N)-1).*N);
					zidx2  = [1:N]'+((levs(1:N)).*N);
					z1 = reshape(tmp.z(zidx1),[lx,ly,lt]);
					z2 = reshape(tmp.z(zidx2),[lx,ly,lt]);
					v1 = reshape(tmp.var(zidx1),[lx,ly,lt]);
					v2 = reshape(tmp.var(zidx2),[lx,ly,lt]);
					z1(mask==0) = NaN;
					z2(mask==0) = NaN;
					v1(mask==0) = NaN;
					v2(mask==0) = NaN;

					% Perform linear interpolation
					vnew(:,:,z,:) = [v2.*(z1-zdep(z)) + v1.*(zdep(z)-z2)]./[(z1-z2)];
				end
				% Save slices
				obj.data.avg.(vars{i}).slice  = squeeze(vnew);
			end
		end % end method zslice_data
	end % end methods declarations
	%----------------------------------------------------------------------------------------

	%----------------------------------------------------------------------------------------
	% Utility functions (static methods)
	methods (Static)
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
		function [tmpslice,tmpmask,tmpdepth] = lonlat_slice(dims,lonlat,data,mask,deg,depth,fillmat,I,LI)
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
			nt   = dims(3);
			ns   = dims(4);
			nx   = dims(5);
			ny   = dims(6);

			% Set up arrays
			tmpslice = NaN(dmsn,nz,nt,ns);
			tmpdepth = NaN(dmsn,nz,nt,ns);
			tmpmask  = NaN(dmsn,nz,nt,ns);
			
			% Start parpool
			if I == 1;
				tic
				delete(gcp('nocreate'));
				parpool(12);
			end

			% Display progress
			disp([num2str(I),'/',num2str(LI)]);

			% Slice!
			%parfor rcrd = 1:nl
			for rcrd = 1:nt
				% - Only get data for that month
				tmpdata = squeeze(data(:,:,:,rcrd,:));
				% - Interp to lon/lat line
				for i = 1:dmsn
					for j = 1:ns
						% - Fill tmpslice
						for z = 1:nz
							if dmsn == nx
								tmpslice(i,z,rcrd,j) = interp1(squeeze(lonlat(i,:)),squeeze(tmpdata(i,:,z)),deg(j));
								tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(i,:)),squeeze(mask(i,:,z)),deg(j));
							elseif dmsn == ny
								tmpslice(i,z,rcrd,j) = interp1(squeeze(lonlat(:,i)),squeeze(tmpdata(:,i,z)),deg(j));
								tmpmask(i,z,rcrd,j)  = interp1(squeeze(lonlat(:,i)),squeeze(mask(:,i,z)),deg(j));
							end
						end
					end
				end
			end
			tmpdepth = permute(repmat(depth,[1,dmsn,nt,ns]),[2 1 3 4]); 
			if I == LI
				delete(gcp('nocreate'));
				toc
			end
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

		%--------------------------------------------------------------------------------
		function [lambda2,xi,ST,SN]=okubo_weiss(u,v,pm,pn); 
			% --------------------------------------------------------------------
			% Compute the OkuboWeiss parameter
			%
			% Usage:
			%   [lambda,x,ST,SN] = okubo_weiss(u,v,pm,pn);
			%
			% Inputs:
			%   u = ROMS u velocity
			%   v = ROMS v velocity
			%   pm = ROMS grid 'pm', or curvilinear coordinate metrix in XI
			%   pn = ROMS grid 'pn', or curvilinear coordinate metrix in ETA
			%
			% Outputs:
			%	lambda2 = Okubo-Weiss (lambda^2) in s^-2
			%   xi      = Relative vorticity in s^-1
			%   ST      = Shear strain in s^-1
			%   SN      = Normal strain in s^-1
			% --------------------------------------------------------------------
			
			% Get grid dimensions		
			[Mp,Lp]=size(pm);
			L=Lp-1; M=Mp-1;
			Lm=L-1; Mm=M-1;
			
			% Get u/v dimensions
			nt = size(u,3);

			% Initialize output matrices
			xi      = zeros(Mp,Lp,nt);  
			mn_p    = zeros(M,L,nt);
			ST      = zeros(Mp,Lp,nt);	
			SN      = zeros(Mp,Lp,nt);
			lambda2 = zeros(Mp,Lp,nt);
			SN      = zeros(Mp,Lp,nt);

			% Initialize output matrices
			for i = 1:nt;
				uom  =zeros(M,Lp); 
				von  =zeros(Mp,L);
				uom=2*u(:,1:L,i)./(pm(:,1:L)+pm(:,2:Lp));
				uon=2*u(:,1:L,i)./(pn(:,1:L)+pn(:,2:Lp));
				von=2*v(1:M,:,i)./(pn(1:M,:)+pn(2:Mp,:));
				vom=2*v(1:M,:,i)./(pm(1:M,:)+pm(2:Mp,:));
				mn=pm.*pn;
				mn_p=(mn(1:M,1:L)+mn(1:M,2:Lp)+...
					  mn(2:Mp,2:Lp)+mn(2:Mp,1:L))/4;

				% relative vorticity
				xi(:,:,i) = mn.*romsMaster.psi2rho(von(:,2:Lp)-von(:,1:L)-uom(2:Mp,:)+uom(1:M,:));

				% Sigma_T
				ST(:,:,i) = mn.*romsMaster.psi2rho(von(:,2:Lp)-von(:,1:L)+uom(2:Mp,:)-uom(1:M,:));

				% Sigma_N
				SN(2:end-1,2:end-1,i) = mn(2:end-1,2:end-1).*(uon(2:end-1,2:end)...
									  -uon(2:end-1,1:end-1)...
									  -vom(2:end,2:end-1)...
									  +vom(1:end-1,2:end-1));
				% Lambda^2
				lambda2(:,:,i) = SN(:,:,i).^2 + ST(:,:,i).^2 - xi(:,:,i).^2;
			end
		end % end static method okubo_weiss

        %--------------------------------------------------------------------------------
		function [var_rho] = psi2rho(var_psi)
			% --------------------------------------------------------------------
			% Transfert a field at psi points to the rho points
			%
			% Usage:
			% - [var_rho] = psi2rho(var_psi)
            % --------------------------------------------------------------------

			% Convert
			[M,L]=size(var_psi);
			Mp=M+1;
			Lp=L+1;
			Mm=M-1;
			Lm=L-1;
			var_rho=zeros(Mp,Lp);
			var_rho(2:M,2:L)=0.25*(var_psi(1:Mm,1:Lm)+var_psi(1:Mm,2:L)+...
								   var_psi(2:M,1:Lm)+var_psi(2:M,2:L));
			var_rho(1,:)=var_rho(2,:);
			var_rho(Mp,:)=var_rho(M,:);
			var_rho(:,1)=var_rho(:,2);
			var_rho(:,Lp)=var_rho(:,L);
		end % end static method psi2rho

        %--------------------------------------------------------------------------------
        function [matrixOut] = spatial_filter(matrixIn,Nr,Nc)
            % --------------------------------------------------------------------
			% Smooths 2D array data, ignores NaNs.
			%	
			% This function smooths the data in matrixIn using a mean filter over a
			% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
			% element "i" by the mean of the rectange centered on "i".  Any NaN
			% elements are ignored in the averaging.  If element "i" is a NaN, then it
			% will be preserved as NaN in the output.  At the edges of the matrix,
			% where you cannot build a full rectangle, as much of the rectangle that
			% fits on your matrix is used (similar to the default on Matlab's builtin
			% function "smooth").
			%
			% Usage:
			% - matrixOut = spatial_filter(matrixIn,Nr,Nc)
			% 
			% Inputs:
			% - matrixIn: original matrix
			% - Nr      : number of points used to smooth rows
			% - Nc      : number of points to smooth columns.  If not specified, Nc = Nr.
			% 
			% Outputs:
			% - matrixOut: smoothed version of original matrix

			% Initial error statements and definitions
			if nargin < 2, error('Not enough input arguments!'), end

			N(1) = Nr;
			if nargin < 3, N(2) = N(1); else N(2) = Nc; end

			if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
			if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

			% Building matrices that will compute running sums.  The left-matrix, eL,
			% smooths along the rows.  The right-matrix, eR, smooths along the
			% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
			% (2*Nc+1) rectangle centered on element "i".
			[row,col] = size(matrixIn);
			eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
			eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

			% Setting all "NaN" elements of "matrixIn" to zero so that these will not
			% affect the summation.  (If this isn't done, any sum that includes a NaN
			% will also become NaN.)
			A = isnan(matrixIn);
			matrixIn(A) = 0;

			% For each element, we have to count how many non-NaN elements went into
			% the sums.  This is so we can divide by that number to get a mean.  We use
			% the same matrices to do this (ie, "eL" and "eR").
			nrmlize = eL*(~A)*eR;
			nrmlize(A) = NaN;

			% Actually taking the mean.
			matrixOut = eL*matrixIn*eR;
			matrixOut = matrixOut./nrmlize;
		end % end static method spatial_filter
	
        %--------------------------------------------------------------------------------
        function [out,grd] = scatter_density(xdata,ydata,xbounds,ybounds)
            % --------------------------------------------------------------------
			% Probability heatmap given xbounds and ybounds, based on
			% independent X variable and dependent Y variable
			% Usage: 
			% - [out,grd] = scatter_density(xdata,ydata,xbounds,ybounds)  
			% 
			% Inputs:
			% - xdata, ydata     = independent, dependent variables (matrices ok)
			% - xbounds, ybounds = increasing xbin/ybin limits
			%
			% Outputs:
			% - out = count, probability, mean, and median results for each xbin
			% - grd = meshgrid for plotting
            % --------------------------------------------------------------------

			% Initialize matrices
			try
				dat = xdata(:)+ydata(:);
			catch
				disp('debug scatter_density');
				keyboard
			end
			xdata = xdata(~isnan(dat)); 
			ydata = ydata(~isnan(dat));
			xdata = xdata(:);
			ydata = ydata(:);
			xbounds = xbounds(:); ybounds = ybounds(:);
			out.count  = nan(length(ybounds)-1,length(xbounds)-1);
			out.prob   = out.count;
			out.mean   = nan(1,length(xbounds)-1);
			out.median = out.mean;

			% Get output grid
			grd.xgrid = 0.5*(xbounds(1:end-1)+xbounds(2:end));
			grd.ygrid = 0.5*(ybounds(1:end-1)+ybounds(2:end));
			[grd.X grd.Y]  = meshgrid(grd.xgrid,grd.ygrid);

			% Get 2D histogram
			for i = 1:length(xbounds)-1
				xi = xbounds(i);
				xf = xbounds(i+1);
				indx = find(xi<=xdata & xdata<xf);
				if isempty(indx)
					out.count(:,i) = zeros(length(ybounds)-1,1);
					out.prob(:,i)  = zeros(length(ybounds)-1,1);;
					out.mean(1,i)  = NaN;
					out.median(1,i) = NaN;
					continue
				else
					tmpdat = ydata(indx);
					for j = 1:length(ybounds)-1
						yi = ybounds(j);
						yf = ybounds(j+1);
						if j == 1
							indy = find(tmpdat<yf);
						elseif j == length(ybounds)-1
							indy = find(tmpdat>=yi);
						else
							indy = find(yi<=tmpdat & tmpdat <yf); 
						end
						if isempty(indy)
							out.count(j,i) = 0;
							out.prob(j,i)  = 0;
						else
							out.count(j,i) = length(indy);
							out.prob(j,i)  = 100*(length(indy)./length(tmpdat));
						end
					end
					out.mean(1,i)   = nanmean(tmpdat);
					out.median(1,i) = nanmedian(tmpdat); 
				end
			end
		end % end static method scatter_density

        %--------------------------------------------------------------------------------
		function [o2_sat] = o2_sat(pT,s)
			% --------------------------------------------------------------------
			% Saturated O2 in mmol/m3
			%
			% Usage:
			%	[o2_sat] = o2_sat(pT,s)
			%
			% Inputs:
			%	pT = potential temperature
			%	s  = salinity;
			% --------------------------------------------------------------------

			% Constants
			a_0 = 2.00907;
			a_1 = 3.22014;
			a_2 = 4.05010;
			a_3 = 4.94457;
			a_4 = -2.56847E-1;
			a_5 = 3.88767;
			b_0 = -6.24523E-3;
			b_1 = -7.37614E-3;
			b_2 = -1.03410E-2;
			b_3 = -8.17083E-3;
			c_0 = -4.88682E-7;

			TS = log(((273.16+25.0)-pT)./(273.16+pT));
			O2SAT = exp(a_0+TS.*(a_1+TS.*(a_2+TS.*(a_3+TS.*(a_4+TS.*a_5)))) + ...
					s.*((b_0+TS.*(b_1+TS.*(b_2+TS.*b_3))) + s.*c_0));
			o2_sat = O2SAT.*44.6596;
		end % end static method o2_sat

        %--------------------------------------------------------------------------------
		function [n2o_sat] = n2o_sat(pT,s)
			% --------------------------------------------------------------------
			% Saturated N2O in mmol/m3
			%
			% Usage:
			%	[n2o_sat] = n2o_sat(pT,s)
			%
			% Inputs:
			%	pT = potential temperature
			%	s  = salinity;
			% --------------------------------------------------------------------

			% Constants
			a_1 = -165.8802;
			a_2 = 222.8743;
			a_3 = 92.0792;
			a_4 = -1.48425;
			b_1 = -0.056235;
			b_2 = 0.031619;
			b_3 = -0.0048472;

			TS = 273.16 + pT;
			N2OSAT =  exp(a_1+a_2.*(100.0./TS)+a_3.*log(TS./100.0)+ ...
						  a_4.*(TS./100.0).^2+s.*(b_1+b_2.*(TS./100.0)+b_3.*(TS./100.0).^2));
			n2o_sat = N2OSAT.*1000.*1000.*(300.0*1e-9);
		end % end static method n2o_sat

        %--------------------------------------------------------------------------------
		function [ammox,nitrox] = oxic_offline(vars_in,PAR,params,mask,dims);
			% --------------------------------------------------------------------
			% Calculates AMMOX and NITROX offline
			% NOTE:
			%	Must call par_offline beforehand, since it is used as an input here
			%
			% Usage:
			%	[ammox,nitrox] = oxic_offline(vars_in,PAR,params,mask,dims)
			%
			% Inputs:
			%	vars_in  = Structure containing O2, NH4, and NO2 data (i.e. vars_in.O2) 
			%	PAR      = Offline PAR calculation (see par_offline)
			%	params   = NitrOMZ parameter scheme (i.e. params.kao, params.nitrif_par_lim)
			%	mask     = Land-mask for 3D data (usually in obj.grid.mask_rho)
			%   dims     = [nx,ny,nz,nt] (use NaN for empty fields)
			%
			% Outputs:
			%	ammox  = NH4 oxidation (offline)
			%	nitrox = NO2 oxidation (offline)
			% --------------------------------------------------------------------
	
			% Process dims
			nx  = dims(1);
			ny  = dims(2);
			nz  = dims(3);
			nt  = dims(4);
			ind = isnan(dims);
			dims(ind==1) = [];
			dimsname = {'nx','ny','nz','nt'};
			dimsname(ind==1) = [];		

			% Loads O2, NH4, and NO2
			% Load 3D filtered fields
			o2  = vars_in.O2.*mask;
			nh4 = vars_in.NH4.*mask;
			no2 = vars_in.NO2.*mask;

			% Filter bad tracer values
			o2(o2<0) = 0;
			nh4(nh4<0) = 0; 
			no2(no2<0) = 0;

			% idx1, PAR.out suggests no light limitation at bottom of cell, nh4 and o2>0
			% idx2, PAR.in suggests there is light limitation at top of cell 
			idx1 = PAR.out<params.nitrif_par_lim & nh4>0 & o2>0;
			idx2 = PAR.out<params.nitrif_par_lim & nh4>0 & o2>0 & PAR.in>params.nitrif_par_lim;
			ind1 = find(idx1==1);
			ind2 = find(idx1==1 & idx2==1);
			ammox = zeros(size(o2)).*mask;
			ammox(ind1) = params.kao.*(o2(ind1)./(o2(ind1) + params.ko2_ao)).*(nh4(ind1)./(nh4(ind1) + params.knh4_ao));
			ammox(ind2) = ammox(ind2).*log(PAR.out(ind2)./params.nitrif_par_lim)./(-PAR.kPARdz(ind2));

			% idx1, PAR.out suggests no light limitation at bottom of cell, no2 and o2>0
			% idx2, PAR.in suggests there is light limitation at top of cell 
			nitrox = zeros(size(o2)).*mask;
			idx1 = PAR.out<params.nitrif_par_lim & no2>0 & o2>0;
			idx2 = PAR.out<params.nitrif_par_lim & no2>0 & o2>0 & PAR.in>params.nitrif_par_lim;
			ind1 = find(idx1==1);
			ind2 = find(idx1==1 & idx2==1);
			nitrox(ind1) = params.kno.*(o2(ind1)./(o2(ind1) + params.ko2_no)).*(no2(ind1)./(no2(ind1) + params.kno2_no));
			nitrox(ind2) = nitrox(ind2).*log(PAR.out(ind2)./params.nitrif_par_lim)./(-PAR.kPARdz(ind2));
		end % end static method oxic_offline

        %--------------------------------------------------------------------------------
		function [denitrif1,denitrif2,denitrif3,anammox] = anoxic_offline(vars_in,params,mask,dims);
			% --------------------------------------------------------------------
			% Calculates DENITRIF1, DENITRIF2, DENITRIF3, and ANAMMOX offline 
			%
			% Usage:
			%	[denitrif1,denitrif2,denitrif3,anammox] = anoxic_offline(vars_in,params,mask);
			%
			% Inputs:
			%	vars_in  = Structure containing O2, NH4, NO3, NO2, N2O, POC_REMIN, 
			%	           DOC_REMIN, SED_DENITRIF, and OTHER_REMIN (i.e. mean_out.O2) 
			%	params   = NitrOMZ parameter scheme (i.e. params.kao, params.nitrif_par_lim)
			%	mask     = Land-mask for 3D data (usually in obj.grid.mask_rho)
			%   dims     = [nx,ny,nz,nt] (use NaN for empty fields)
			%
			% Outputs:
			%	denitrif1 = NO3 reduction (offline) 
			%	denitrif2 = NO2 reduction (offline)
			%	denitrif3 = N2O reduction (offline)
			%	anammox   = anaerboic NH4 oxidation (offline)
			% --------------------------------------------------------------------

			% Process dims
			nx  = dims(1);
			ny  = dims(2);
			nz  = dims(3);
			nt  = dims(4);
			ind = isnan(dims);
			dims(ind==1) = [];
			dimsname = {'nx','ny','nz','nt'};
			dimsname(ind==1) = [];

			% Very small number
			epsN = 10^(-30);

			% Load 3D filtered fields
			nh4 = vars_in.NH4.*mask;
			no3 = vars_in.NO3.*mask;
			no2 = vars_in.NO2.*mask;
			n2o = vars_in.N2O.*mask;
			o2  = vars_in.O2.*mask;
			poc = vars_in.POC_REMIN.*mask;
			doc = vars_in.DOC_REMIN.*mask;	

			% Filter bad tracer values
			o2(o2<0) = 0;
			nh4(nh4<0) = 0;
			no2(no2<0) = 0;
			no3(no3<0) = 0;
			n2o(n2o<0) = 0;

			% Load 2D filtered fields
			sed = vars_in.SED_DENITRIF;
			oth = vars_in.OTHER_REMIN;

			% Get offline calculations
			Rden = (poc + doc); 

			% Add sed denitrif and other remin processes from sediment
			% Only do this if depth is supplied
			if ~isnan(nz)
				ind = find(strcmp(dimsname,'nz')==1);
				if ind == 1
					Rden(1,:) = Rden(1,:)-(sed.*params.denitrif_C_N)-oth;
				elseif ind == 2
					Rden(:,1,:) = Rden(:,1,:)-(sed.*params.denitrif_C_N)-oth;
				elseif ind == 3
					Rden(:,:,1) = Rden(:,:,1)-(sed.*params.denitrif_C_N)-oth;
				else
					Rden(:,:,:,1) = Rden(:,:,:,1)-(sed.*params.denitrif_C_N)-oth;
				end
			end

			% Offline calculation
			d0 = zeros(size(o2));
			R_oxic = params.koxic.*(o2./(params.ko2_oxic+o2)) + epsN;
			R_den1 = params.kden1.*(no3./(params.kno3_den1+no3)).*exp(-max(o2,d0)./params.ko2_den1) + epsN;
			R_den2 = params.kden2.*(no2./(params.kno2_den2+no2)).*exp(-max(o2,d0)./params.ko2_den2) + epsN;
			R_den3 = params.kden3.*(n2o./(params.kn2o_den3+n2o)).*exp(-max(o2,d0)./params.ko2_den3) + epsN;
			roxic = R_oxic./(R_oxic+R_den1+R_den2+R_den3);
			rden1 = R_den1./(R_oxic+R_den1+R_den2+R_den3);
			rden2 = R_den2./(R_oxic+R_den1+R_den2+R_den3);
			rden3 = R_den3./(R_oxic+R_den1+R_den2+R_den3);

			% Anammox
			idx1 = nh4>0 & no2>0;
			anammox = zeros(size(o2));
			anammox(idx1) = params.kax.*(nh4(idx1)./(params.knh4_ax+nh4(idx1))).*...
										(no2(idx1)./(params.kno2_ax+no2(idx1))).*...
										exp(-max(o2(idx1),d0(idx1))./params.ko2_ax) ;

			% Denitrif1
			idx1 = no3>0; 
			denitrif1 = zeros(size(o2));
			denitrif1(idx1) = rden1(idx1).*Rden(idx1).*params.denitrif_N_C;

			% Denitrif2
			idx1 = no2>0;
			denitrif2 = zeros(size(o2));
			denitrif2(idx1) = rden2(idx1).*Rden(idx1).*params.denitrif_N_C;

			% Denitrif3
			idx1 = n2o>0;
			denitrif3 = zeros(size(o2));
			denitrif3(idx1) = rden3(idx1).*Rden(idx1).*params.denitrif_N_C;
		end % end static method denitrif_offline

        %--------------------------------------------------------------------------------
		function [n2oammox,no2ammox] = n2oammox_offline(vars_in,params,ammox,mask);
			% --------------------------------------------------------------------
			% Calculates N2O and NO2 partitioning from AMMOX
			% i.e. some AMMOX goes to NO2, the rest to N2O following Nevison et al., 2003	
			% NOTE:
			%	Usually paired with oxic_offline beforehand, since ammox is used as an input
			%	However, it can also be applied to the online AMMOX calculation...
			% 
			% Usage:
			%	[n2oammox,no2ammox] = n2oammox_offline(vars_in,params,ammox,mask)
			%
			% Inputs:
			%	vars_in  = Structure containing O2 data (i.e. vars_in.O2) 
			%	params   = NitrOMZ parameter scheme (i.e. params.kao, params.nitrif_par_lim)
			%	mask     = Land-mask for 3D data (usually in obj.grid.mask_rho)
			%
			% Outputs:
			%	n2oammox = N2O production from NH4 oxidation
			%	no2ammox = NO2 production from NO2 oxidation
			% --------------------------------------------------------------------

			% Very small number
			epsN = 10^(-30);

			% Load 3D filtered fields
			o2 = vars_in.O2.*mask;

			% Filter bad tracer values
			o2(o2<0) = 0;

			% Offline calculation
			% Get O2-dependent fractions to N2O/NO2 (accounting for 2Ns in N2O)
			idx1 = o2>- 1e-5;
			n2o_nh4_yield1 = 0.5*ones(size(o2));
			no2_nh4_yield1 = 0.0*ones(size(o2));
			n2o_nh4_yield1(idx1) = 1./(2.0+2.0./((params.n2o_ji_a./o2(idx1)+params.n2o_ji_b)./100));
			no2_nh4_yield1(idx1) = 1./(1.*((params.n2o_ji_a./o2(idx1)+params.n2o_ji_b)./100)+1);

			% Get rates
			n2oammox = ammox.*n2o_nh4_yield1;
			no2ammox = ammox.*no2_nh4_yield1;
		end % end static method n2oammox_offline

        %--------------------------------------------------------------------------------
		function [n2oao1_cons,n2osiden_cons,n2osoden_cons,n2oatm_cons] = n2ocons_offline(vars_in,denitrif3,mask);
			% --------------------------------------------------------------------
			% Calculates the decomposed N2O reduction (since we occassionally decompose N2O i.e. McCoy et al, 2023) 
			% NOTE:
			%	Usually paired with anoxic_offline beforehand, since denitrif3 is used as an input
			%	However, it can also be applied to the online denitrif3 calculation...
			% 
			% Usage:
			%	[n2oao1_cons,n2osiden_cons,n2osoden_cons,n2oatm_cons] = n2ocons_offline(vars_in,denitrif3,mask);	
			%
			% Inputs:
			%	vars_in  = Structure containing decomposed N2O tracers 
			%	           (i.e. vars_in.N2O_AO1, vars_in.N2O_SODEN, vars_in.N2O_SIDEN, vars_in.N2O_ATM) 
			%	params   = NitrOMZ parameter scheme (i.e. params.kao, params.nitrif_par_lim)
			%	mask     = Land-mask for 3D data (usually in obj.grid.mask_rho)
			%
			% Outputs:
			%	n2oao1_cons   = Reduction of N2O from nitrification
			%	n2osiden_cons = Reduction of N2O from boundaries
			%	n2osoden_cons = Reduction of N2O from denitrification
			%	n2oatm_cons   = Reductino of N2O from atmosphere
			% --------------------------------------------------------------------

			% Very small number
			epsN = 10^(-30);

			% Load 3D filtered fields
			n2oao1   = vars_in.N2O_AO1.*mask;
			n2osiden = vars_in.N2O_SIDEN.*mask; 
			n2osoden = vars_in.N2O_SODEN.*mask; 
			n2oatm   = vars_in.N2O_ATM.*mask;

			% Filter bad tracer fields	
			n2oao1(n2oao1<0) = 0;
			n2osiden(n2osiden<0) = 0;
			n2osoden(n2osoden<0) = 0;
			n2oatm(n2oatm<0) = 0;

			% Get total n2o
			n2otot = n2oao1 + n2osiden + n2osoden + n2oatm + 4*epsN;

			% Offline calculation
			n2oao1_cons   = ((n2oao1+epsN)./n2otot).*denitrif3;
			n2osiden_cons = ((n2osiden+epsN)./n2otot).*denitrif3;
			n2osoden_cons = ((n2osoden+epsN)./n2otot).*denitrif3;
			n2oatm_cons   = ((n2oatm+epsN)./n2otot).*denitrif3;
		end % end static method n2ocons_offline

        %--------------------------------------------------------------------------------
		function [PAR] = par_offline(vars_in,gridfile,params,mask)
			% --------------------------------------------------------------------
			% Calculates PAR offline using PARinc and TOT_CHL
			% NOTE:
			%	This calculation is based off of the ETH ROMS
			%	Update this for UCLA ROMS!
			%
            % Usage:
            %	[PAR] = par_offline(vars_in,gridfile,params,mask) 
            %
            % Inputs:
            %	vars_in  = Structure containing TOT_CHL, PARinc, and zeta (i.e. vars_in.PARinc) 
			%	gridfile = Full path and filename of ROMS gridfile
            %	params   = zlevs4 parameters for simulation (i.e. params.theta_s) 
            %	mask     = Land-mask for 3D data (usually in obj.grid.mask_rho)
			%
			% Outputs:
			%	PAR = Structure containing PAR.in, PAR.out, and PAR.lay (i.e. averaged over layer) 
			% --------------------------------------------------------------------
			% Load 3D filtered fields
			totchl = vars_in.TOT_CHL.*mask;
			
			% Load 2D filtered fields
			parinc = vars_in.PARinc.*mask;

			% Load zeta (for dz calc)
			zeta = vars_in.zeta.*mask;
			
			% Load bathymetry
			h = ncread(gridfile,'h').*mask;

			% Get z_w and dzb
			zw = zlevs4(h,zeta,params.theta_s,params.theta_b,params.hc,params.nz,'w','new2012');
			zw = permute(zw,[2 3 1]);
			dz = diff(zw,1,3); 

			% Initialize to NaN (not 0, which would allow for rates)
			PAR.out = nan(size(dz));
			PAR.in  = nan(size(dz));
			PAR.lay = nan(size(dz));

			% Get attenuation coefficient, apply ETH Zurich(?) fixes
			% E3SM
			% https://github.com/E3SM-Project/Ocean-BGC/blob/master/BGC_mod.F90
			PAR.kPARdz = max(totchl,0.02.*ones(size(totchl))).*mask;
			PAR.kPARdz(PAR.kPARdz<0.13224)  = 0.0919.*((PAR.kPARdz(PAR.kPARdz<0.13224).^0.3536));
			PAR.kPARdz(PAR.kPARdz>=0.13224) = 0.1131.*((PAR.kPARdz(PAR.kPARdz>=0.13224).^0.4562)); 
			PAR.kPARdz = PAR.kPARdz.*dz;
			
			% Get surface attenuation
			par.kpardz = squeeze(PAR.kPARdz(:,:,end));	
			par.parin = parinc;											    % top of cell is just incoming PAR
			par.parout = parinc.*exp(-par.kpardz);						    % attenuate over top layer
			PAR.in(:,:,end)  = par.parin;								    % incoming
			PAR.out(:,:,end) = par.parout;								    % attenuate
			PAR.lay(:,:,end) = par.parin.*(1-exp(-par.kpardz))./par.kpardz; % average over layer
			for z = size(dz,3)-1:-1:1
				par.kpardz = squeeze(PAR.kPARdz(:,:,z));				    % layer attenuation
				par.parin  = squeeze(PAR.out(:,:,z+1));					    % in == out of above cell
				par.parout = par.parin.*exp(-par.kpardz);					% attentuate to get new out
				par.parlay = par.parin.*(1-exp(-par.kpardz))./par.kpardz;	% average over layer
				PAR.in(:,:,z)  = par.parin.*mask;							% update PAR.in
				PAR.out(:,:,z) = par.parout.*mask;							% update PAR.out
				PAR.lay(:,:,z) = par.parlay.*mask;							% update PAR.lay
			end
		end % end static method par_offline
	end % end static methods declarations
	%----------------------------------------------------------------------------------------
end % end classdef
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
function [params] = getParams 
	% ----------------------
	% Get N-cycle parameters
    filedir = ['/data/project1/demccoy/iNitrOMZ/analysis/run/bgc_params.mat'];
    tmp = load(filedir);
    params = tmp.params;
end % end method getParams

%--------------------------------------------------------------------------------
function [constants] = getConstants
	% ----------------------
	% Save constants
    constants.mmolNs_TgNy  = [14*86400*365.25]/(1e15); % mmolN/s to TgN/yr
    constants.mmolN2s_TgNy = [28*86400*365.25]/(1e15); % mmolN2/s to TgN/yr
    constants.s_d          = 86400;                    % seconds to days
    constants.mmol_umol    = 1000;                     % mmol to umol
    constants.umol_nmol    = 1000;                     % umol to nmol
    constants.mmol_nmol    = 1000*1000;                % mmol to nmol
    constants.nmol_umol    = 1/1000;                   % nmol to umol
    constants.umol_mmol    = 1/1000;                   % umol to mmol
    constants.nmol_mmol    = (1/1000)*(1/1000);        % nmol to mmol
end % end method getConstants 

%--------------------------------------------------------------------------------
function [paths] = getDiagPaths;	
	% -------------------
	% Initialize diagnostic plots: set metadata, save paths, get axis limits 
	% temperature
	paths.diag.temp.file  = {'/data/project3/data/woa18/temperature/0p25/temp_woa18_clim.nc'};
	paths.diag.temp.type  = {'nc'};
	paths.diag.temp.var   = {'temp'};
	paths.diag.temp.zvar  = {'depth'};
	paths.diag.temp.dim   = {'xyzt'};
	paths.diag.temp.lon   = {'lon'};
	paths.diag.temp.lat   = {'lat'};
	paths.diag.temp.name  = {'WOA-18 Temperature'};
	paths.diag.temp.units = {'$^oC$'};

	% salinity
	paths.diag.salt.file  = {'/data/project3/data/woa18/salinity/0p25/salt_woa18_clim.nc'};
	paths.diag.salt.type  = {'nc'};
	paths.diag.salt.var   = {'salt'};
	paths.diag.salt.zvar  = {'depth'};
	paths.diag.salt.dim   = {'xyzt'};
	paths.diag.salt.lon   = {'lon'};
	paths.diag.salt.lat   = {'lat'};
	paths.diag.salt.name  = {'WOA-18 Salinity'};
	paths.diag.salt.units = {'$PSU$'};

	% density
	paths.diag.sigma.file  = {'/data/project3/data/woa18/density/0p25/sigma_woa18_clim.nc'};
	paths.diag.sigma.type  = {'nc'};
	paths.diag.sigma.var   = {'sigma'};
	paths.diag.sigma.zvar  = {'depth'};
	paths.diag.sigma.dim   = {'xyzt'};
	paths.diag.sigma.lon   = {'lon'};
	paths.diag.sigma.lat   = {'lat'};
	paths.diag.sigma.name  = {'WOA-18 Density'};
	paths.diag.sigma.units = {'$kg$ $m^{-3}$'};

	% sea surface height
	paths.diag.SSH.file  = {'/data/project1/demccoy/ROMS/validation/AVISO/monthly_AVISO.mat'};
	paths.diag.SSH.type  = {'mat'};
	paths.diag.SSH.var   = {'adt_month_av'};
	paths.diag.SSH.zvar  = {[]};
	paths.diag.SSH.dim   = {'xyt'};
	paths.diag.SSH.lon   = {'lon_aviso'};
	paths.diag.SSH.lat   = {'lat_aviso'};
	paths.diag.SSH.name  = {'AVISO Absolute Dynamic Topography'};
	paths.diag.SSH.units = {'$m$'};

	% zonal wind stress
	paths.diag.zws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
	paths.diag.zws.type  = {'mat'};
	paths.diag.zws.var   = {'u'};
	paths.diag.zws.zvar  = {[]};
	paths.diag.zws.dim   = {'xyt'};
	paths.diag.zws.lon   = {'lon'};
	paths.diag.zws.lat   = {'lat'};
	paths.diag.zws.name  = {'SCOW-2010 Zonal Wind Stress'};
	paths.diag.zws.units = {'$N$ $m^{-2}$'};

	% meridional wind stress
	paths.diag.mws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
	paths.diag.mws.type  = {'mat'};
	paths.diag.mws.var   = {'v'};
	paths.diag.mws.zvar  = {[]};
	paths.diag.mws.dim   = {'xyt'};
	paths.diag.mws.lon   = {'lon'};
	paths.diag.mws.lat   = {'lat'};
	paths.diag.mws.name  = {'SCOW-2010 Meridional Wind Stress'};
	paths.diag.mws.units = {'$N$ $m^{-2}$'};
	
	% wind stress
	paths.diag.ws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
	paths.diag.ws.type  = {'mat'};
	paths.diag.ws.var   = {'ws'};
	paths.diag.ws.zvar  = {[]};
	paths.diag.ws.dim   = {'xyt'};
	paths.diag.ws.lon   = {'lon'};
	paths.diag.ws.lat   = {'lat'};
	paths.diag.ws.name  = {'SCOW-2010 Wind Stress'};
	paths.diag.ws.units = {'$N$ $m^{-2}$'};

	% wind stress curl
	paths.diag.wsc.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
	paths.diag.wsc.type  = {'mat'};
	paths.diag.wsc.var   = {'wsc'};
	paths.diag.wsc.zvar  = {[]};
	paths.diag.wsc.dim   = {'xyt'};
	paths.diag.wsc.lon   = {'lon'};
	paths.diag.wsc.lat   = {'lat'};
	paths.diag.wsc.name  = {'SCOW-2010 Wind Stress Curl'};
	paths.diag.wsc.units = {'$N$ $m^{-2}$'};

	% u velocity
	paths.diag.u.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
	paths.diag.u.type  = {'mat'};
	paths.diag.u.var   = {'U'};
	paths.diag.u.zvar  = {'dep'};
	paths.diag.u.dim   = {'xyzt'};
	paths.diag.u.lon   = {'lon'};
	paths.diag.u.lat   = {'lat'};
	paths.diag.u.name  = {'GODAS U-Velocity'};
	paths.diag.u.units = {'$m$ $s^{-1}'}; 

	% v velocity
	paths.diag.v.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
	paths.diag.v.type  = {'mat'};
	paths.diag.v.var   = {'V'};
	paths.diag.v.zvar  = {'dep'};
	paths.diag.v.dim   = {'xyzt'};
	paths.diag.v.lon   = {'lon'};
	paths.diag.v.lat   = {'lat'};
	paths.diag.v.name  = {'GODAS V-Velocity'};
	paths.diag.v.units = {'$m$ $s^{-1}'}; 

	% oxygen (only 1.00 available)
	paths.diag.O2.file  = {'/data/project3/data/woa18/oxygen/1p0/o2_woa18_clim.nc'};
	paths.diag.O2.type  = {'nc'};
	paths.diag.O2.var   = {'o2'};
	paths.diag.O2.zvar  = {'depth'};
	paths.diag.O2.dim   = {'xyzt'};
	paths.diag.O2.lon   = {'lon'};
	paths.diag.O2.lat   = {'lat'};
	paths.diag.O2.name  = {'WOA-18 Oxygen'};
	paths.diag.O2.units = {'$mmol$ $m^{-3}$'};

	% nitrate (only 1.00 available)
	paths.diag.NO3.file  = {'/data/project3/data/woa18/nitrate/1p0/no3_woa18_clim.nc'};
	paths.diag.NO3.type  = {'nc'};
	paths.diag.NO3.var   = {'no3'};
	paths.diag.NO3.zvar  = {'depth'};
	paths.diag.NO3.dim   = {'xyzt'};
	paths.diag.NO3.lon   = {'lon'};
	paths.diag.NO3.lat   = {'lat'};
	paths.diag.NO3.name  = {'WOA-18 Nitrate'};
	paths.diag.NO3.units = {'$mmol$ $m^{-3}$'};
	
	% Copy for NOx (WOA-18 NO3 is actually NO3+NO2)
	paths.diag.NOX = paths.diag.NO3;
	paths.diag.NOX.name = {'WOA-18 $NO_x$'};

	% phosphate (only 1.00 available)
	paths.diag.PO4.file  = {'/data/project3/data/woa18/phosphate/1p0/po4_woa18_clim.nc'};
	paths.diag.PO4.type  = {'nc'};
	paths.diag.PO4.var   = {'po4'};
	paths.diag.PO4.zvar  = {'depth'};
	paths.diag.PO4.dim   = {'xyzt'};
	paths.diag.PO4.lon   = {'lon'};
	paths.diag.PO4.lat   = {'lat'};
	paths.diag.PO4.name  = {'WOA-18 Phosphate'};
	paths.diag.PO4.units = {'$mmol$ $m^{-3}$'};

	% N* (only 1.00 available)
	paths.diag.nstar.file  = {'/data/project3/data/woa18/nstar/1p0/nstar_woa18_clim.nc'};
	paths.diag.nstar.type  = {'nc'};
	paths.diag.nstar.var   = {'nstar'};
	paths.diag.nstar.zvar  = {'depth'};
	paths.diag.nstar.dim   = {'xyzt'};
	paths.diag.nstar.lon   = {'lon'};
	paths.diag.nstar.lat   = {'lat'};
	paths.diag.nstar.name  = {'WOA-18 N*'};
	paths.diag.nstar.units = {'$mmol$ $m^{-3}$'};
				 
	% silicate (only 1.00 available)
	paths.diag.SiO3.file  = {'/data/project3/data/woa18/silicate/1p0/si_woa18_clim.nc'};
	paths.diag.SiO3.type  = {'nc'};
	paths.diag.SiO3.var   = {'si'};
	paths.diag.SiO3.zvar  = {'depth'};
	paths.diag.SiO3.dim   = {'xyzt'};
	paths.diag.SiO3.lon   = {'lon'};
	paths.diag.SiO3.lat   = {'lat'};
	paths.diag.SiO3.name  = {'WOA-18 Silicate'};
	paths.diag.SiO3.units = {'$mmol$ $m^{-3}$'};

	% N2O (only 1.00 available, reconstructed)
	%paths.diag.N2O.file   = {'/data/project1/dclements/N_Interior/n2opredRF_02-21-2021.nc'};
	%paths.diag.N2O.file   = {'/data/project2/yangsi/analysis/n2oInterior/processed/n2opredRF_05-27-2020.nc'};
	%paths.diag.N2O.type   = {'nc'};
	paths.diag.N2O.file   = {'/data/project1/demccoy/ROMS/validation/n2o/n2o_NN_format.mat'};
	paths.diag.N2O.type   = {'mat'};
	paths.diag.N2O.var    = {'n2o'};
	paths.diag.N2O.zvar   = {'depth'};
	paths.diag.N2O.dim    = {'xyzt'};
	paths.diag.N2O.lon    = {'lon'};
	paths.diag.N2O.lat    = {'lat'};
	paths.diag.N2O.name   = {'Clements et al. Nitrous Oxide'};
	paths.diag.N2O.name   = {'ML Nitrous Oxide Estimate'};
	paths.diag.N2O.units  = {'$mmol$ $m^{-3}$'};
	paths.diag.N2O.factor = {0.001}; %nmol to mmol

	% NO2 (only 1.00 available, reconstructed)
	%paths.diag.NO2.file  = {'/data/project1/demccoy/ROMS/validation/no2/no2_NN_format.mat'};
	%paths.diag.NO2.type  = {'mat'};
	paths.diag.NO2.file  = {'/data/project2/yangsi/analysis/no2Interior/processed/no2predRF_05-27-2020.nc'};
	paths.diag.NO2.type  = {'nc'};
	paths.diag.NO2.var   = {'no2'};
	paths.diag.NO2.zvar  = {'depth'};
	paths.diag.NO2.dim   = {'xyzt'};
	paths.diag.NO2.lon   = {'lon'};
	paths.diag.NO2.lat   = {'lat'};
	paths.diag.NO2.name  = {'Clements et al. Nitrite'};
	paths.diag.NO2.name  = {'ML Nitrite Estimate'};
	paths.diag.NO2.units = {'$mmol$ $NO_2 m^{-3}$'};

	% Chl 
	paths.diag.SFC_CHL.file  = {'/data/project1/data/MODIS-Aqua/backup/res0p25/chl_clim_0p25.nc'};
	paths.diag.SFC_CHL.type  = {'nc'};
	paths.diag.SFC_CHL.dim   = {'xyt'};
	paths.diag.SFC_CHL.var   = {'chl'};
	paths.diag.SFC_CHL.zvar  = {[]};
	paths.diag.SFC_CHL.lon   = {'lon'};
	paths.diag.SFC_CHL.lat   = {'lat'};
	paths.diag.SFC_CHL.name  = {'MODIS-Aqua Chlorophyll'};
	paths.diag.SFC_CHL.units = {'$mg$ $C$ $m^{-3}$'};
		
	% NPP products (VGPM, CBPM, CAFE)	
	paths.diag.NPP.file  = {'/data/project2/yangsi/analysis/NPPcode/std-VGPM/std-VGPMnpp_MODIS_clim2002-2018.nc',...
								'/data/project2/yangsi/analysis/NPPcode/CbPM2/CbPM2npp_MODIS_clim2002-2018.nc',...
								'/data/project1/data/MODIS-Aqua/CAFE_NPP/climatology/nppclim_CAFE_MODIS_9km.nc'}; 
	paths.diag.NPP.type  = {'nc','nc','nc'};
	paths.diag.NPP.var   = {'npp','npp','npp'};
	paths.diag.NPP.zvar  = {[],[],[]};
	paths.diag.NPP.dim   = {'xyt','xyt','xyt'};
	paths.diag.NPP.lon   = {'lon','lon','lon'};
	paths.diag.NPP.lat   = {'lat','lat','lat'};
	paths.diag.NPP.name  = {'NPP-VGPM','NPP-CBPM','NPP-CAFE'};
	paths.diag.NPP.units = {'$mg$ $C$ $m^{-2}$ $d^{-1}$','$mg$ $C$ $m^{-2}$ $d^{-1}$','$mg$ $C$ $m^{-2}$ $d^{-1}$'};

	% MLD products (WOCE/NODC/ARGO or ARGO only);
	paths.diag.MLD.file  = {'/data/project1/demccoy/ROMS/validation/MLD/Argo_mixedlayers_monthlyclim_12112019.nc',...
								'/data/project1/demccoy/ROMS/validation/MLD/mld_DR003_c1m_reg2.0.nc'};
	paths.diag.MLD.type  = {'nc','nc'};
	paths.diag.MLD.var   = {'mld_da_mean','mld'};
	paths.diag.MLD.zvar  = {[],[]};
	paths.diag.MLD.dim   = {'txy','xyt'};
	paths.diag.MLD.lon   = {'lon','lon'};
	paths.diag.MLD.lat   = {'lat','lat'};
	paths.diag.MLD.name  = {'Argo Mixed Layer Depth','IFREMER Mixed Layer Depth'}; 
	paths.diag.MLD.units = {'$m$','$m$'};

	% POC_FLUX_IN
	paths.diag.POC_FLUX_IN.file   = {'/data/project1/demccoy/ROMS/validation/POC_FLUX_IN/clements_100m_flux.mat'};
	paths.diag.POC_FLUX_IN.file   = {'/data/project1/demccoy/ROMS/validation/POC_FLUX_IN/clements_100m_flux.mat'};
	paths.diag.POC_FLUX_IN.type   = {'mat'};
	paths.diag.POC_FLUX_IN.var    = {'flux_mean'};
	paths.diag.POC_FLUX_IN.zvar   = {[]};
	paths.diag.POC_FLUX_IN.dim    = {'xyt'};
	paths.diag.POC_FLUX_IN.lon    = {'lon'};
	paths.diag.POC_FLUX_IN.lat    = {'lat'};
	paths.diag.POC_FLUX_IN.name   = {'Clements et al. (2022) 100m POC Flux'};
	paths.diag.POC_FLUX_IN.units  = {'$mmol$ $C$ $m^{-2}$ $s^{-1}$'};
	paths.diag.POC_FLUX_IN.factor = {(1/12.01*86400)}; % mgC/m2/d to mmolC/m2/s 

end % end method initPlots

