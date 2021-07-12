% Preprocessing Steps
% Replace strings inbetween $...$
% 1) create version folder under simulation folder 
%
%    > cd /data/project1/demccoy/ROMS/$r_tag$/
%    > mkdir $vers$
%
% 2) create subdirectories (clim, grid, zcoord)
%
%    > cd $vers$
%    > mkdir clim grid zcoord
%
% 3) link monthly data from parent directory 
%
%    > cd clim
%    > ln -s $data_directory$/avg_Y$year$*.nc .
%
% 4) concatenate monthly files into yearly file
%
%    > ncrcat avg_Y$year$*.nc avg_$year$.nc
%
% 5) link grid file into grid directory
%
%    > cd ../grid
%    > ln -s $grid_dir$/$grid_fname$ .
%
% 6) link zcoord.mat file
%
%    > cd ../zcoord
%    > ln -s $zcoord_dir$/$zcoord_fname$ .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear figures
close all

% - romsMaster directory (update for new users!!!!)
RMdir = '/data/project1/demccoy/ROMS/';	

% - initiate simslist
simslist = []; simscnt = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 NYF4
% meta info
simslist{simscnt}    = ['peru_NYF4']; simscnt = simscnt + 1;
peru_NYF4.r_tag      = ['peru_chile_0p1']; % simulation
peru_NYF4.vers       = ['dccoy_NYF4']; % version
peru_NYF4.ext        = ['2009']; 
peru_NYF4.simPath    = ['/data/project2/model_output/',peru_NYF4.r_tag,'/',[peru_NYF4.r_tag,'_',peru_NYF4.vers,'/']];
peru_NYF4.avgfile    = ['avg_',peru_NYF4.ext,'.nc']; % file to diagnose
peru_NYF4.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_NYF4.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_NYF4.year       = [2009];

% region info (indices)
peru_NYF4.rlati  = [ ];%[ 321 462];
peru_NYF4.rloni  = [ ];%[  11 341];
peru_NYF4.rdepi  = [ ];%[-750 inf];
peru_NYF4.rcoast = [ ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 CVKV2
% meta info
simslist{simscnt}    = ['peru_CKV2']; simscnt = simscnt + 1;
peru_CKV2.r_tag      = ['peru_chile_0p1']; % simulation
peru_CKV2.vers       = ['dccoy_CKV2']; % version
peru_CKV2.ext        = ['2009'];
peru_CKV2.simPath    = ['/data/project2/model_output/',peru_CKV2.r_tag,'/',[peru_CKV2.r_tag,'_',peru_CKV2.vers,'/']];
peru_CKV2.avgfile    = ['avg_',peru_CKV2.ext,'.nc']; % file to diagnose
peru_CKV2.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_CKV2.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_CKV2.year       = [2009];

% region info (indices)
peru_CKV2.rlati  = []; %[ 321 462];
peru_CKV2.rloni  = []; %[  11 341];
peru_CKV2.rdepi  = []; %[-750 inf];
peru_CKV2.rcoast = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p VKV1
% meta info
simslist{simscnt}    = ['peru_VKV1']; simscnt = simscnt + 1;
peru_VKV1.r_tag      = ['peru_chile_0p1']; % simulation
peru_VKV1.vers       = ['dccoy_VKV1']; % version
peru_VKV1.ext        = ['2009'];
peru_VKV1.simPath    = ['/data/project2/model_output/',peru_VKV1.r_tag,'/',[peru_VKV1.r_tag,'_',peru_VKV1.vers,'/']];
peru_VKV1.avgfile    = ['avg_',peru_VKV1.ext,'.nc']; % file to diagnose
peru_VKV1.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_VKV1.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_VKV1.year       = [2009];

% region info (indices)
peru_VKV1.rlati  = []; %[ 321 462];
peru_VKV1.rloni  = []; %[  11 341];
peru_VKV1.rdepi  = []; %[-750 inf];
peru_VKV1.rcoast = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV3
% meta info
simslist{simscnt}    = ['peru_VKV3']; simscnt = simscnt + 1;
peru_VKV3.r_tag      = ['peru_chile_0p1']; % simulation
peru_VKV3.vers       = ['dccoy_VKV3']; % version
peru_VKV3.ext        = ['2009'];
peru_VKV3.simPath    = ['/data/project2/model_output/',peru_VKV3.r_tag,'/',[peru_VKV3.r_tag,'_',peru_VKV3.vers,'/']];
peru_VKV3.avgfile    = ['avg_',peru_VKV3.ext,'.nc']; % file to diagnose
peru_VKV3.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_VKV3.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_VKV3.year       = [2009];

% region info (indices)
peru_VKV3.rlati  = []; %[ 321 462];
peru_VKV3.rloni  = []; %[  11 341];
peru_VKV3.rdepi  = []; %[-750 inf];
peru_VKV3.rcoast = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4
% meta info
simslist{simscnt}    = ['peru_VKV4']; simscnt = simscnt + 1;
peru_VKV4.r_tag      = ['peru_chile_0p1']; % simulation
peru_VKV4.vers       = ['dccoy_VKV4']; % version
peru_VKV4.ext        = ['2009'];
peru_VKV4.simPath    = ['/data/project2/model_output/',peru_VKV4.r_tag,'/',[peru_VKV4.r_tag,'_',peru_VKV4.vers,'/']];
peru_VKV4.avgfile    = ['avg_',peru_VKV4.ext,'.nc']; % file to diagnose
peru_VKV4.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_VKV4.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_VKV4.year       = [2009];

% region info (indices)
peru_VKV4.rlati  = []; %[ 321 462];
peru_VKV4.rloni  = []; %[  11 341];
peru_VKV4.rdepi  = []; %[-750 inf];
peru_VKV4.rcoast = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE1 (lower Ji's, higher KNo, lower KO2Den2)
% meta info
simslist{simscnt}          = ['peru_VKV4_TUNE1']; simscnt = simscnt + 1;
peru_VKV4_TUNE1.r_tag      = ['peru_chile_0p1']; % simulation
peru_VKV4_TUNE1.vers       = ['dccoy_VKV4_tune1']; % version
peru_VKV4_TUNE1.ext        = ['2009'];
peru_VKV4_TUNE1.simPath    = ['/data/project2/model_output/',peru_VKV4_TUNE1.r_tag,'/',[peru_VKV4_TUNE1.r_tag,'_',peru_VKV4_TUNE1.vers,'/']];
peru_VKV4_TUNE1.avgfile    = ['avg_',peru_VKV4_TUNE1.ext,'.nc']; % file to diagnose
peru_VKV4_TUNE1.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_VKV4_TUNE1.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_VKV4_TUNE1.year       = [2009];

% region info (indices)
peru_VKV4_TUNE1.rlati  = []; %[ 321 462];
peru_VKV4_TUNE1.rloni  = []; %[  11 341];
peru_VKV4_TUNE1.rdepi  = []; %[-750 inf];
peru_VKV4_TUNE2.rcoast = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE2 (lower Ji's, higher KNo, lower KO2Den2, double KAo)
% meta info
simslist{simscnt}          = ['peru_VKV4_TUNE2']; simscnt = simscnt + 1;
peru_VKV4_TUNE2.r_tag      = ['peru_chile_0p1']; % simulation
peru_VKV4_TUNE2.vers       = ['dccoy_VKV4_tune2']; % version
peru_VKV4_TUNE2.ext        = ['2010']; % year of budget
peru_VKV4_TUNE2.simPath    = ['/data/project2/model_output/',peru_VKV4_TUNE2.r_tag,'/',[peru_VKV4_TUNE2.r_tag,'_',peru_VKV4_TUNE2.vers,'/']];
peru_VKV4_TUNE2.avgfile    = ['avg_',peru_VKV4_TUNE2.ext,'.nc']; % file to diagnose
peru_VKV4_TUNE2.gfile      = ['/data/project1/yangsi/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];
peru_VKV4_TUNE2.zcfile     = ['/data/project2/yangsi/analysis/ROMS_coords/peru_chile_0p1/peru_chile_0p1_zcoords.mat'];
peru_VKV4_TUNE2.year       = [2010];

% region info (indices)
peru_VKV4_TUNE2.rlati  = [ 321 452];%[ 321 462];
peru_VKV4_TUNE2.rloni  = [ 21  281];%[  11 341];
peru_VKV4_TUNE2.rdepi  = [-750 inf];
peru_VKV4_TUNE2.rcoast = [100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([RMdir,'current_sims.mat'],simslist{:})
clear(simslist{:});
clear simslist simscnt RMdir

%{
% pacmed_0p25 2006 (more rates, less tracers)
% meta info
r_tag      = ['pacmed_0p25']; % simulation
vers       = ['25TRY']; % version
ext        = ['2006']; % from 2010
simPath    = ['/data/project1/demccoy/ROMS/',r_tag,'/',vers,'/']; % path to data and folders
avgfile    = ['avg_',ext,'.nc']; % file to diagnose
gfile      = [simPath,'clim/',r_tag,'_grd.nc']; % grid file
zcfile     = [];

% region info (indices)
rlati = [];
rloni = [];
rdepi = [];

% pacmed_0p25 2010 (more tracers, less rates)
% meta info
r_tag      = ['pacmed_0p25']; % simulation
vers       = ['25TRY']; % version
ext        = ['2010']; % from 2010
simPath    = ['/data/project1/demccoy/ROMS/',r_tag,'/',vers,'/']; % path to data and folders
avgfile    = ['avg_',ext,'.nc']; % file to diagnose
gfile      = [simPath,'clim/',r_tag,'_grd.nc']; % grid file
zcfile     = [];

% region info (indices)
rlati = [];
rloni = [];
rdepi = [];

%}



