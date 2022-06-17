
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
			A.lonbounds = [floor(min(obj.grid.lon_rho(:))) ceil(max(obj.grid.lon_rho(:)))];
			A.latbounds = [floor(min(obj.grid.lat_rho(:))) ceil(max(obj.grid.lat_rho(:)))];
			A.tfont     = 18;
			A.cbfont    = 14;
			A           = parse_pv_pairs(A,varargin); % parse method arguments to A
			
			% Process inputs
			if isempty(A.time)
				A.time = 1:1:obj.grid.nt;
			end

			% First generate uniform color bars for each axis	
			% Go through each variables
			for i = 1:length(vars)
				
				% Go through each time-record (or dont, based on input)
				for t = 1:length(A.time)
					tt = A.time(t);	
					
					% Gather terms
					dcdt = obj.data.budget.(vars{i}).intdcdt(:,:,tt);
					adv  = obj.data.budget.(vars{i}).intadx(:,:,tt) + ...
						   obj.data.budget.(vars{i}).intady(:,:,tt) + ...
						   obj.data.budget.(vars{i}).intadz(:,:,tt);
					dfz  = obj.data.budget.(vars{i}).intdfz(:,:,tt);
					sms  = obj.data.budget.(vars{i}).intsms(:,:,tt);
					fg   = obj.data.budget.(vars{i}).intfg(:,:,tt);
					sed  = obj.data.budget.(vars{i}).intsed(:,:,tt);
					net  = obj.data.budget.(vars{i}).intnet(:,:,tt);
					
					% Blank land
					dcdt = dcdt .* obj.grid.mask_rho;
					adv  = adv  .* obj.grid.mask_rho;
					dfz  = dfz  .* obj.grid.mask_rho;
					sms  = sms  .* obj.grid.mask_rho;
					fg   = fg   .* obj.grid.mask_rho;
					sed  = sed  .* obj.grid.mask_rho;
					net  = net  .* obj.grid.mask_rho;

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
					dcdt = obj.data.budget.(vars{i}).intdcdt.(:,:,tt);
					adv  = obj.data.budget.(vars{i}).intadx.(:,:,tt) + ...
						   obj.data.budget.(vars{i}).intady.(:,:,tt) + ...
						   obj.data.budget.(vars{i}).intadz.(:,:,tt);
					dfz  = obj.data.budget.(vars{i}).intdfz.(:,:,tt);
					sms  = obj.data.budget.(vars{i}).intsms.(:,:,tt);
					fg   = obj.data.budget.(vars{i}).intfg.(:,:,tt);
					sed  = obj.data.budget.(vars{i}).intsed.(:,:,tt);
					net  = obj.data.budget.(vars{i}).intnet.(:,:,tt);

					% Blank land
					dcdt = dcdt .* obj.grid.mask_rho;
					adv  = adv  .* obj.grid.mask_rho;
					dfz  = dfz  .* obj.grid.mask_rho;
					sms  = sms  .* obj.grid.mask_rho;
					fg   = fg   .* obj.grid.mask_rho;
					sed  = sed  .* obj.grid.mask_rho;
					net  = net  .* obj.grid.mask_rho;
					
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
						[ax] = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,...
									'lonbounds',A.lonbounds,'latbounds',A.latbounds);
						clevs        	    = clevs(1):(diff(clevs)/31):clevs(2);
						dat(dat<clevs(1))   = clevs(1);
						dat(dat>clevs(end)) = clevs(end);
						[tmp hc]            = m_contourf(obj.grid.lon_rho,obj.grid.lat_rho,dat,clevs,'LineStyle','none');
						% Plot region box
						m_plot(obj.grid.lon_rho(1,:),obj.grid.lat_rho(1,:),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(:,1),obj.grid.lat_rho(:,1),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'--k','linewidth',2);
						m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'--k','linewidth',2);
						cb = colorbar;
						title({titstr,[obj.info.time_string_idv{tt},' ',num2str(obj.info.runYear)]},...
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
		function plot1D(obj,vars,file,varargin)
			% ------------------
			% - Plot the relevant rates versus depth at a lat,lon point
			%
			% Usage:
			% - [fig] = plot1D(obj,vars,varargin);
			%
			% Inputs:
			% - vars = variable(s) to plot
			%
			% Inputs (varargin):
			% - time - time record to plot (default--> all)
			% - lon  - longitude point (can be a vector)
			% - lat  - latitude point (can be a vector)
			% - zlim - limit in depth space (i.e. [1200] for 0 - 1200)
			%
			% - Example:
			%   plot1D(obj,'vars',{'NO2','N2O'},'time',1:12,'lon',200,'lat',0)
			% ------------------

			disp('UPDATE FOR VARIABLE NT');
			return

			% Defaults for optional  arguments
			A.time  = [];
			A.lon   = [];
			A.lat   = [];
			A.zlim  = [];
			A       = parse_pv_pairs(A,varargin); % parse method arguments to A

			% Process inputs
			if isempty(A.time)
				A.time = 1:obj.grid.nt;
			end 	
			if isempty(A.lon) | isempty(A.lat)
				disp('No lat/lon point selected, killing');
				return
			end
			
			% Load data
			for i = 1:length(vars)
				try
					obj = loadData(obj,vars(i),'grid','raw');	
				catch
					obj = computeVar(obj,vars(i),'type','z_avg');
				end
			end

			% Get index of nearest lat/lon grid cell
			lon_idx = [];
			lat_idx = [];
			for i = 1:length(A.lon)
				all_dist   = distance(A.lat(i),A.lon(i),obj.grid.lat_rho,obj.grid.lon_rho);
				[lon_idx(i),lat_idx(i)]  = find(all_dist == min(all_dist(:)));
			end

			% Get colormix
			clrs = colormap(cmocean('phase'));
			idx  = floor(linspace(1,length(clrs),length(A.time)+1));
			clrs = clrs(idx,:); clrs(end,:) = [];

			% Process time
			if obj.grid.nt == 12
				tstr = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
			elseif obj.grid.nt == 12 & length(A.time) < 12
				ttstr = [];
				tstr = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
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
						tmpdata = squeeze(obj.data.(A.vars{j}).data(lon_idx(i),lat_idx(i),:,A.time(k)));

						% Gather depth
						depth = -squeeze(obj.grid.z_r(lon_idx(i),lat_idx(i),:,A.time(k)));

						% Make plot(s)
						if k == 1
							fig = piofigs('lfig',1.5);
							set(gca,'YDir','Reverse');
							xlabel(obj.data.(A.vars{j}).units);
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
						if length(A.time) < obj.grid.nt
							fname = [A.vars{j},'_vs_z_Lon_',num2str(A.lon(i)),'_Lat_',num2str(A.lat(i)),'_',tstr{A.time(k)}];
							print('-djpeg',[obj.paths.plots.roms.profile,fname])
							close all
						end
					end
					
					% Add averages if many months plotted
					if length(A.time) == obj.grid.nt
						% Gather data
						tmpdata = squeeze(obj.data.(A.vars{j}).data(lon_idx(i),lat_idx(i),:,:));
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
						tmp.(vars{i}) = nan([obj.grid.nx obj.grid.ny length(obj.grid.z_avg_dep) obj.grid.nt length(yr_range)]);
					else
						tmp.(vars{i}) = nan([obj.grid.nx obj.grid.ny length(A.lvl) obj.grid.nt length(yr_range)]);
					end
				else % 2D variable
					tmp.(vars{i}) = nan([obj.grid.nx obj.grid.ny obj.grid.nt]);
				end
				for j = 1:length(yr_range)
					% Get filename
					current_file = ['avg_',num2str(yr_range(j))];
					disp(['...',num2str(yr_range(j))]);
					if ~ismember(vars{i},obj.info.var2d)
						try
							obj = loadData(obj,vars(i),'type','z_avg','file',current_file);
						catch
							obj = changeInputs(obj,'runYear',yr_range(j));
							obj = computeVar(obj,vars(i),'type','z_avg');
						end
						if isempty(A.lvl)
							tmp.(vars{i})(:,:,:,:,j) = obj.data.(A.type).(vars{i}).data; 
						else
							ind = [];
							for k = 1:length(A.lvl)
								ind(k) = find(obj.grid.z_avg_dep == A.lvl(k));
							end
							tmp.(vars{i})(:,:,:,:,j) = obj.data.(A.type).(vars{i}).data(:,:,ind,:);   
						end
						obj.data.(A.type).(vars{i}).data = [];
					else
						obj                         = loadData(obj,vars(i),'grid','raw','file',current_file);
						tmp.(vars{i})(:,:,:,j)      = obj.data.(A.type).(vars{i}).data; 
						obj.data.(A.type).(vars{i}).data = [];
					end
				end
				if ~ismember(vars{i},obj.info.var2d)
					obj.data.(A.type).(vars{i}).data = nanmean(tmp.(vars{i}),5);
				else
					obj.data.(A.type).(vars{i}).data = nanmean(tmp.(vars{i}),4);
				end
			end
		end % end method getAvgData	

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
			romsu  = ncread(obj.paths.z_avg,'u',[obj.info.(A.type).ind4D(1,:)],[obj.info.(A.type).ind4D(2,:)]);
											   
			romsln = obj.grid.lon_rho;
			romslt = obj.grid.lat_rho;
			romsdp = obj.grid.woa1p0.depth;
			
			% Take slice of data for each month
			dmsn = obj.grid.ny;
			nl   = obj.info.(A.type).time(file);
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
			obj.data.u.slice  = cat(3,gridu{:});
			obj.data.u.name  = obj.info.name3d{idx};
			obj.data.u.units  = obj.info.unit3d{idx};
			obj.diag.u.slice  = datu./100;
			obj.diag.u.name   = diagname;
			obj.diag.u.units  = obj.info.unit3d{idx};
			obj.slice.deg   = latu;
			obj.slice.depth = depu;
			if nanmean(depu) < 0
				obj.slice.depth = -depu;
			end
			obj.slice.coord = 'longitude';
			obj.slice.sect  = sect_avg;
		end % end method equatorUcmp
