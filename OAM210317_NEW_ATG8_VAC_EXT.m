

% This version has problems with V_loc and W5_loc. Most of the locations
% were missing.
parpool('local',12)

type='.tif';
path_h='/archive/bioinformatics/Danuser_lab/Quiescence/shared/experiments/';
%path_h='Z:\shared\experiments\';
exp_name='OAM_201207_ATG8_cdc14'; 
final_time=300;
positions=[0:3 7:9 12 15 17:22];% 
disp(positions)
% BF TFP GFP/mNG mKok mRuby mNeptune CYOFP
% =========================================================================
vac_modifier=1;
STD_Thr=2;

%============flags for markers - used in allocations==========================================
%I cannot change these flags, because it gives an error inside the parfor
%loop.
%============load segmentation=============================================
for pos_here=positions
    % pos_here=0
    ptmp_h=pos_here;
     path_seg='//archive/bioinformatics/Danuser_lab/Quiescence/shared/segmented_extracted_data/';
   %path_seg='Z:\shared\segmented_extracted_data\';
    load([path_seg exp_name '/' exp_name '_pos' num2str(ptmp_h) '_' num2str(final_time)])
    
    x_size =520;   %1040;
    y_size =692;   %1388;
    
    appr_vol                      =zeros(no_obj,numbM);
    cell_area                     =zeros(no_obj,numbM);
    %---------------------------------------------------------------------
    %Allocation of features extracted from all fluorophores regardless of the
    %protein tagged
        mean_int_per_area_BF        =zeros(no_obj,numbM);
        Total_int_per_area_BF         =zeros(no_obj,numbM);
     
        all_obj.med_G                =zeros(1,numbM);
        Total_int_NG                 =zeros(no_obj,numbM);
        mean_int_NG                  =zeros(no_obj,numbM);

    %Allocation of protein specific features
    %Vacuolar focci===========================================================================
        Non_Vac_mean_int_ATG8                           =zeros(no_obj,numbM);
        Non_Vac_STD_int_ATG8                            =zeros(no_obj,numbM); 
        Non_Vac_Total_ATG8                              =zeros(no_obj,numbM);   
        Non_Vac_focci                                   =zeros(no_obj,numbM);
        Non_Vac_focci_area                              =zeros(no_obj,numbM);
        Non_Vac_focci_cellarea_ratio                    =zeros(no_obj,numbM);
        Non_Vac_focci_Total_int                         =zeros(no_obj,numbM);
        Non_Vac_focci_mean_int                          =zeros(no_obj,numbM);
        Non_Vac_focci_median_int                        =zeros(no_obj,numbM);
        Non_Vac_focci_std_int                           =zeros(no_obj,numbM);
    
    
Vac_mean_int_ATG8                               =zeros(no_obj,numbM);
Vac_STD_int_ATG8                                =zeros(no_obj,numbM); 
Vac_Total_ATG8                                  =zeros(no_obj,numbM);   
Vac_focci                                       =zeros(no_obj,numbM);
Vac_focci_area                                  =zeros(no_obj,numbM);
Vac_focci_cellarea_ratio                        =zeros(no_obj,numbM);
Vac_focci_Total_int                             =zeros(no_obj,numbM);
Vac_focci_mean_int                              =zeros(no_obj,numbM);
Vac_focci_median_int                            =zeros(no_obj,numbM);
Vac_focci_std_int                               =zeros(no_obj,numbM);
   
Vac_Total_BF                     =zeros(no_obj,numbM); % total GFP concentration
Vac_mean_int_BF                  =zeros(no_obj,numbM); 
Vac_no_obj                       =zeros(no_obj,numbM);
NON_Vac_mean_int_BF              =zeros(no_obj,numbM);
NON_Vac_Total_BF                 =zeros(no_obj,numbM); % total GFP concentration
   

    tmp_sizes=zeros(1,no_obj); %allocation
    Lcells1=all_obj.cells(:,:,numbM);
    parfor i5=1:no_obj
        tmp_sizes(i5)=tmp_sizes(i5)+sum(sum(Lcells1==i5));
    end
    max_allowed_cell_size=max_size_vs_largest_cell*max(tmp_sizes);
    s=round(sqrt(max_allowed_cell_size))+50; %this will determine the size of the matrix we will put put_vac etc into
    
    for c_time=numbM:-1:1
       %  c_time=7;
        disp(c_time)
        Lcells=all_obj.cells(:,:,c_time); %a broadcast variable
        
        %Variables to be used in construction of V_loc E6_loc nucl_mC W5_loc nucl_K
        V_loc_tmp=zeros(s,s,no_obj);
        location_cell=zeros(no_obj,4);
 
        %Find Median Fluorescence
        image_number=sprintf('%09d',c_time);
 
            IBF=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_BrightField_000.tif']);
            IBF=double(IBF);
%             imagesc(IBF)
            IBF=medfilt2(IBF,'symmetric');

            IG=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_505 mNG_000.tif']); 
            IG=double(IG);
            IG=medfilt2(IG,'symmetric');
            backgr_G=(((IG)+1).*(~Lcells));
            backgr_G=sort(backgr_G(:));
            [tmpV,posH]=max(backgr_G>0);
            backgr_G=median(backgr_G(max([1 posH-1]):end))-1;
            IG=IG-backgr_G;
            all_obj.med_G(1,c_time)  =backgr_G;
       
        parfor cell_no=1:no_obj
                %  cell_no=202
             if c_time>=cell_exists(cell_no,2)
                
                ccell=(Lcells==cell_no); %a temporary variable
                appr_vol(cell_no,c_time)  = appr_vol(cell_no,c_time)+(sum(ccell(:))).^(1.5); % 2D area ^(3/2) %reduction variable
                % figure; imagesc(Lcells)
                %Get Cell Location
                cell_margin=1;
                [x_cn,y_cn]=get_wind_coord1_EZGI(ccell,cell_margin);
                location_cell(cell_no,:)=location_cell(cell_no,:)+[y_cn(1) y_cn(end) x_cn(1) x_cn(end)];
                
                %Vacuole first
                
                ccell2=double(ccell(y_cn,x_cn)).*double(IBF(y_cn,x_cn));
                 
               Total_int_per_area_BF (cell_no,c_time) = Total_int_per_area_BF (cell_no,c_time)+sum(ccell2(ccell2(:)>0));
                mean_int_per_area_BF (cell_no,c_time) =mean_int_per_area_BF (cell_no,c_time)+mean(ccell2(ccell2(:)>0));
                                     
                % figure;imagesc(ccell2)
                 
                %put inside a function later

                   put_vac=(ccell2>=(vac_modifier.*mean(ccell2(ccell2~=0)))); % manual threshold by vac_mean_modifier
                   % figure;imagesc(put_vac)
                if isempty(put_vac)==0 && sum(put_vac(:))~=0
                       
                   %%%%% PUT VAC OR PUTVAC-2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
%                        put_vac2=retain_largest_obj1(put_vac);
                         put_vac1= bwdist(~put_vac); % linear transformation
                         % figure;imagesc(put_vac1)
                         put_vac2=put_vac1; % bait
                         put_vac2(put_vac2<=1)=0; % destroy spagetti
                         % figure;imagesc(put_vac2)
                       
                         put_vac2b=medfilt2(put_vac2); % get rid of noise pixels
                         % figure;imagesc(put_vac2b)
                        
                         put_vac3 = imgaussfilt(put_vac2b); % smooth the area for vacuole
% %                        % figure;imagesc(put_vac3)
                         se = offsetstrel('ball',5,5); % structuring element
                         put_vac4=imdilate(put_vac3,se); % blurred area to include associated foci 
                         % figure;imagesc(put_vac4)
                         med=median(put_vac4(:)); 
                         put_vac5= put_vac4-med; % make mask with bacjground 0 
                         % figure;imagesc(put_vac5)
% ==========================================================================================
                         put_vac5b=bwlabel(put_vac5); %  label what remains as vaoule
                         % figure;imagesc(put_vac5b)
                         Vac_no_obj(cell_no,c_time) = max( put_vac5b(:)); 
                         put_vac6=put_vac5b;
                         put_vac6(put_vac6>0)=1; % binary mask to multiply other channels
                           % figure;imagesc(put_vac6)
%===========================================================================================
                         % vacuole from brightfield/phase
                         vacuole=double(put_vac6).*double(ccell2); % get vouale pixels in phase contrast figure;imagesc(vacuole)
                         NO_vacuole=double(~put_vac6).*double(ccell2) ; % get non vacuole pixels in contrast  figure;imagesc(NO_vacuole)
                         
                         % vacoule in BF
                         Vac_mean_int_BF(cell_no,c_time)   = Vac_mean_int_BF(cell_no,c_time) +mean(vacuole(vacuole>0));
                         Vac_Total_BF(cell_no,c_time)      = Vac_Total_BF(cell_no,c_time)+sum(vacuole(:));
                          
                         NON_Vac_mean_int_BF(cell_no,c_time)   = NON_Vac_mean_int_BF(cell_no,c_time) +mean(NO_vacuole(NO_vacuole>0));
                         NON_Vac_Total_BF(cell_no,c_time)      = NON_Vac_Total_BF(cell_no,c_time)+sum(NO_vacuole(:));

 %==============================================================================
                
            % NeonGreen

             ccell3=double(ccell(y_cn,x_cn)).*double(IG(y_cn,x_cn)); % obtain vacuole area from fluorescent channel
             % figure;imagesc(ccell3)
             % totals
             Total_int_NG(cell_no,c_time)=sum(ccell3(ccell3>0));
             mean_int_NG(cell_no,c_time)=mean(ccell3(ccell3>0));
             cell_area(cell_no,c_time)=numel(ccell3(ccell3>0));
            
              % caulate vacuolar associated foci   
              vacuole_ATG8=double(put_vac6).*double(ccell3); % get vacuole in NG
              % figure;imagesc(vacuole_ATG8) 
              % general vacuolar features
              Vac_mean_int_ATG8(cell_no,c_time) = Vac_mean_int_ATG8(cell_no,c_time) + mean(vacuole_ATG8(vacuole_ATG8>0));
              Vac_STD_int_ATG8 (cell_no,c_time) = Vac_STD_int_ATG8(cell_no,c_time) + std(vacuole_ATG8(vacuole_ATG8>0));
              Vac_Total_ATG8(cell_no,c_time)    = Vac_Total_ATG8(cell_no,c_time)+sum(vacuole_ATG8(:));
              % Get the number of foci and signals
              focci_inside=find(vacuole_ATG8<= (mean(vacuole_ATG8(vacuole_ATG8>0))+ STD_Thr*std(vacuole_ATG8(vacuole_ATG8>0))));% find the foci by intensity
              vacuole_ATG8_2=vacuole_ATG8; % bait
              vacuole_ATG8_2(focci_inside)=0;  % left only the foci figure;imagesc(vacuole_ATG8_2);
            
              vacuole_ATG8_3=bwlabel(vacuole_ATG8_2); % ; labelled objects
              % figure;imagesc(vacuole_ATG8_3);
              Vac_focci(cell_no,c_time) =max(vacuole_ATG8_3(:)); % get number of objects in the vacuole
              
              Vac_focci_area(cell_no,c_time)=numel(vacuole_ATG8_3(vacuole_ATG8_3>0)); % count all pixels id'ed as vacoular
              Vac_focci_cellarea_ratio(cell_no,c_time)= (Vac_focci_area(cell_no,c_time)/cell_area(cell_no, c_time))*100; % vacuolar area covered by pixels
              Vac_focci_Total_int(cell_no, c_time)=sum(vacuole_ATG8_2(vacuole_ATG8_2>0));
              Vac_focci_mean_int(cell_no, c_time)=nanmean(vacuole_ATG8_2(vacuole_ATG8_2>0));
              Vac_focci_median_int(cell_no, c_time)=nanmedian(vacuole_ATG8_2(vacuole_ATG8_2>0));
              Vac_focci_std_int(cell_no, c_time)=nanstd(vacuole_ATG8_2(vacuole_ATG8_2>0));
              
              % calculate vacuolar associated foci   
              Non_vacuole_ATG8=double(~put_vac6).*double(ccell3); % get NOnvacuole in NG
              % figure;imagesc( Non_vacuole_ATG8) 
              % general vacuolar features
              Non_Vac_mean_int_ATG8(cell_no,c_time) = Non_Vac_mean_int_ATG8(cell_no,c_time) + mean(Non_vacuole_ATG8(Non_vacuole_ATG8>0));
              Non_Vac_STD_int_ATG8 (cell_no,c_time) = Non_Vac_STD_int_ATG8(cell_no,c_time) + std(Non_vacuole_ATG8(Non_vacuole_ATG8>0));
              Non_Vac_Total_ATG8(cell_no,c_time)    = Non_Vac_Total_ATG8(cell_no,c_time)+sum(Non_vacuole_ATG8(:));
              % Get the number of foci and signals
              focci_outside=find(Non_vacuole_ATG8<= (mean(Non_vacuole_ATG8(Non_vacuole_ATG8>0))+ STD_Thr*std(Non_vacuole_ATG8(Non_vacuole_ATG8>0))));% find the foci by intensity
              Non_vacuole_ATG8_2=Non_vacuole_ATG8; % bait
              Non_vacuole_ATG8_2(focci_outside)=0;  % left only the foci figure;imagesc(vacuole_ATG8_2);
            
              Non_vacuole_ATG8_3=bwlabel(Non_vacuole_ATG8_2); % ; labelled objects
              % figure;imagesc(vacuole_ATG8_3);
              Non_Vac_focci(cell_no,c_time) =max(Non_vacuole_ATG8_3(:)); % get number of objects in the vacuole
              
              Non_Vac_focci_area(cell_no,c_time)=numel(Non_vacuole_ATG8_3(Non_vacuole_ATG8_3>0)); % count all pixels id'ed as vacoular
              Non_Vac_focci_cellarea_ratio(cell_no,c_time)= (Non_Vac_focci_area(cell_no,c_time)/cell_area(cell_no, c_time))*100; % vacuolar area covered by pixels
              Non_Vac_focci_Total_int(cell_no, c_time)=sum(Non_vacuole_ATG8_2(Non_vacuole_ATG8_2>0));
              Non_Vac_focci_mean_int(cell_no, c_time)=nanmean(Non_vacuole_ATG8_2(Non_vacuole_ATG8_2>0));
              Non_Vac_focci_median_int(cell_no, c_time)=nanmedian(Non_vacuole_ATG8_2(Non_vacuole_ATG8_2>0));
              Non_Vac_focci_std_int(cell_no, c_time)=nanstd(Non_vacuole_ATG8_2(Non_vacuole_ATG8_2>0));
              
               else
               end

            end %end if cell exists
        end %parfor
        

    end %time-loop
    
    %add all the fields to all_obj structure
    all_obj.appr_vol                      =appr_vol;
    
    all_obj.mean_int_per_area_BF           =mean_int_per_area_BF;
    all_obj.Total_int_per_area_BF          =Total_int_per_area_BF   ;
    all_obj.Vac_mean_int_BF                =Vac_mean_int_BF;
    all_obj.Vac_no_obj                     =Vac_no_obj;
    all_obj.NON_Vac_mean_int_BF            =NON_Vac_mean_int_BF;
    all_obj.NON_Vac_Total_BF               =NON_Vac_Total_BF;
    all_obj.Vac_conc_T                     =Vac_Total_BF;
 

              all_obj.cell_area       =cell_area;

              all_obj.Total_int_NG       =Total_int_NG;
              all_obj.mean_int_NG      =mean_int_NG;

     
              all_obj.Vac_mean_int_ATG8          =Vac_mean_int_ATG8;
              all_obj.Vac_STD_int_ATG8           =Vac_STD_int_ATG8; 
              all_obj.Vac_Total_ATG8             =Vac_Total_ATG8;   
              all_obj.Vac_focci                  =Vac_focci;
              all_obj.Vac_focci_area             =Vac_focci_area;
              all_obj.Vac_focci_cellarea_ratio   =Vac_focci_cellarea_ratio;
              all_obj.Vac_focci_Total_int        =Vac_focci_Total_int;
              all_obj.Vac_focci_mean_int         =Vac_focci_mean_int;
              all_obj.Vac_focci_median_int       =Vac_focci_median_int;
              all_obj.Vac_focci_std_int          =Vac_focci_std_int;

              
              all_obj.Non_Vac_mean_int_ATG8          =Non_Vac_mean_int_ATG8;
              all_obj.Non_Vac_STD_int_ATG8           =Non_Vac_STD_int_ATG8; 
              all_obj.Non_Vac_Total_ATG8             =Non_Vac_Total_ATG8;   
              all_obj.Non_Vac_focci                  =Non_Vac_focci;
              all_obj.Non_Vac_focci_area             =Non_Vac_focci_area;
              all_obj.Non_Vac_focci_cellarea_ratio   =Non_Vac_focci_cellarea_ratio;
              all_obj.Non_Vac_focci_Total_int        =Non_Vac_focci_Total_int;
              all_obj.Non_Vac_focci_mean_int         =Non_Vac_focci_mean_int;
              all_obj.Non_Vac_focci_median_int       =Non_Vac_focci_median_int;
              all_obj.Non_Vac_focci_std_int          =Non_Vac_focci_std_int;

    path_save='/archive/bioinformatics/Danuser_lab/Quiescence/shared/segmented_extracted_data/';
    name1=[exp_name '/' exp_name '_pos' num2str(ptmp_h) '_VacFoci_LTG'];
    save(fullfile(path_save,name1), 'all_obj','vac_modifier','STD_Thr'); 
end
delete(gcp('nocreate'))


