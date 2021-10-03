% Preprocessing Steps
% Replace strings inbetween $...$
% 1) create version folder under simulation folder 
%
%    > cd /data/project2/model_output/peru_chile_0p1/$r_tag$/
%    > mkdir $vers$
%
% 2) concatenate monthly files into yearly file
%
%    > ncrcat avg_Y$year$*.nc avg_$year$.nc
%
% 3) Update info below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear figures
close all

% - romsMaster directory (update for new users!!!!)
RMdir = '/data/project1/demccoy/ROMS/';	

% - initiate simslist
simslist = []; simscnt = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERU 10km PRODUCTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 NYF4
% meta info
simslist{simscnt} = ['peru_NYF4']; simscnt = simscnt + 1;
tmp.r_tag	      = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_NYF4']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = [ ];%[ 321 462];
tmp.rloni  = [ ];%[  11 341];
tmp.rdepi  = [ ];%[-750 inf];
tmp.rcoast = [ ];

% rename and save
peru_NYF4 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 CVKV2
% meta info
simslist{simscnt} = ['peru_CKV2']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_CKV2']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = []; %[ 321 462];
tmp.rloni  = []; %[  11 341];
tmp.rdepi  = []; %[-750 inf];
tmp.rcoast = [];

% rename and save
peru_CKV2 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p VKV1
% meta info
simslist{simscnt} = ['peru_VKV1']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV1']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = []; %[ 321 462];
tmp.rloni  = []; %[  11 341];
tmp.rdepi  = []; %[-750 inf];
tmp.rcoast = [];

% rename and save
peru_VKV1 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV3
% meta info
simslist{simscnt} = ['peru_VKV3']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV3']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = []; %[ 321 462];
tmp.rloni  = []; %[  11 341];
tmp.rdepi  = []; %[-750 inf];
tmp.rcoast = [];

% rename and save
peru_VKV3 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4
% meta info
simslist{simscnt} = ['peru_VKV4']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV4']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = []; %[ 321 462];
tmp.rloni  = []; %[  11 341];
tmp.rdepi  = []; %[-750 inf];
tmp.rcoast = [];

% rename and save
peru_VKV4 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE1 (lower Ji's, higher KNo, lower KO2Den2)
% meta info
simslist{simscnt} = ['peru_VKV4_TUNE1']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV4_tune1']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = []; %[ 321 462];
tmp.rloni  = []; %[  11 341];
tmp.rdepi  = []; %[-750 inf];
tmp.rcoast = []; 

% rename and save
peru_VKV4_TUNE1 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE2 (lower Ji's, higher KNo, lower KO2Den2, double KAo)
% meta info
simslist{simscnt} = ['peru_VKV4_TUNE2']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV4_tune2']; % version
tmp.year          = [2019];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = [301  461];%[ 321 452];%[ 321 462];
tmp.rloni  = [ 31  341];%[ 21  281];%[  11 341];
tmp.rdepi  = [-750 inf];%[-750 inf];
tmp.rcoast = [20];        %[100];

% rename and save
peru_VKV4_TUNE2 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE2 Ji test(double Ji params from VKv4t2)
% meta info
simslist{simscnt} = ['peru_VKV4_TUNE2_JI_TEST']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV4_tune2_Ji_test']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = [301  461];%[ 321 452];%[ 321 462];
tmp.rloni  = [ 31  341];%[ 21  281];%[  11 341];
tmp.rdepi  = [-750 inf];%[-750 inf];
tmp.rcoast = [20];        %[100];

% rename and save
peru_VKV4_TUNE2_JI_TEST = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peru_chile_0p1 VKV4_TUNE2 Ji test(double Ji params from VKv4t2)
% meta info
simslist{simscnt} = ['peru_VKV4_TUNE2_O2_TEST']; simscnt = simscnt + 1;
tmp.r_tag         = ['peru_chile_0p1']; % simulation
tmp.vers          = ['dccoy_VKV4_tune2_O2_test']; % version
tmp.year          = [2009];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% region info (indices)
tmp.rlati  = [301  461];%[ 321 452];%[ 321 462];
tmp.rloni  = [ 31  341];%[ 21  281];%[  11 341];
tmp.rdepi  = [-750 inf];%[-750 inf];
tmp.rcoast = [20];        %[100];

% rename and save
peru_VKV4_TUNE2_O2_TEST = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PACIFIC 25km PRODUCTS  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pacmed_0p25 NYF4 (Simon's original parameter set)
% meta info
simslist{simscnt} = ['pacmed_NYF4']; simscnt = simscnt + 1;
tmp.r_tag         = ['pacmed_0p25']; % simulation
tmp.vers          = ['dccoy_NYF4']; % version
tmp.year          = [2004];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath		  = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile		  = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/pacmed_0p25/grid/pacmed_0p25_grd_corrected.nc'];

% region info (indices)
tmp.rlati  = [];
tmp.rloni  = [];
tmp.rdepi  = [];
tmp.rcoast = [];

% rename and save
pacmed_NYF4 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pacmed_0p25 VKV4_TUNE2 (lower Ji's, higher KNo, lower KO2Den2, double KAo)
% meta info
simslist{simscnt} = ['pacmed_VKV4_TUNE2']; simscnt = simscnt + 1;
tmp.r_tag         = ['pacmed_0p25']; % simulation
tmp.vers          = ['dccoy_VKV4_tune2']; % version
tmp.year          = [2004];
tmp.ext           = [num2str(tmp.year)]; 
tmp.simPath       = ['/data/project2/model_output/',tmp.r_tag,'/',[tmp.r_tag,'_',tmp.vers,'/']];
tmp.avgfile       = ['avg_',tmp.ext,'.nc']; % file to diagnose
tmp.gfile         = ['/data/project2/model_output/pacmed_0p25/grid/pacmed_0p25_grd_corrected.nc'];

% region info (indices)
tmp.rlati  = [];
tmp.rloni  = [];
tmp.rdepi  = [];
tmp.rcoast = [];

% rename and save
pacmed_VKV4_TUNE2 = tmp; clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([RMdir,'current_sims.mat'],simslist{:})
clear(simslist{:});
clear simslist simscnt RMdir
