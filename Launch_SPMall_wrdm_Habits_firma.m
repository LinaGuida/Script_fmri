%onset image = when image is displayed
% onset write == endtrial_time, does not mean writing moment. this is normally very inmediate after onset image
% fixation_time (duration of fixation) = random duration between 0.4 - 1. to discount from onset
% image and create fixation onset
% probably, we need to enter trial duration in each onset: trial_dur = endtrial_time - onset image


clear; close all; clc;

% source_dir = '/media/lina/Data_Disk/Habits/Early_fmri_data/post/fmri';
% struct_dir = '/media/lina/Data_Disk/Habits/Early_fmri_data/post/structural';
% res_dir =    '/media/lina/Data_Disk/Habits/Early_fmri_data/post/fMRI_SPM_writing';
% subjs = dir([source_dir '/Early_0*']);
% cd(source_dir);
source_dir = '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_habits';
struct_dir = '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI/Structural_habits';
res_dir =    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing';
subjs = dir([source_dir '/Habits_057*']);
cd(source_dir);


for id_sub = 1:length(subjs)
    
    
%     addpath([matlab_tboxfileRef = fullfile(fmrifolderin,'vol1_0001.nii');
%         fileIn = fullfile(fmap_dir,'fmap_rads.nii.gz');
%         fileInRes = fullfile(fmap_dir,'fmap_rads.nii.gz');
%         sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
%         [status,result] = system(sentence)es_path 'spm12'])
%     rmpath([matlab_tboxes_path 'spm12'])
%     addpath([matlab_tboxes_path 'spm12'])
%     
    subjsname = subjs(id_sub).name;
    fprintf('Doing subject %s......\n',subjsname);
    subject_dir = fullfile(source_dir,subjsname);
    
    EPIfolders = dir([subject_dir '/EPI_ESCRITURA_*']);
    if length(EPIfolders)>1
        error('Comprueba carpetas de escritura en este sujeto');
    end
    res_dir_subj = fullfile(res_dir,subjsname);
    if ~exist(res_dir_subj,'dir'); mkdir(res_dir_subj); end
    
    %--------------------------------------------------------------------------
    % Prepare files and folders
    fMRIdir_all = [res_dir_subj '/fMRI_all_firma'];
    mkdir(fMRIdir_all);
    cont = 0; % Counter to store the volumes from all blocks in ascending order
    for i=1:length(EPIfolders)
        sub_tasks{i} = ['/fMRI_B' num2str(i)]; % Store each block folder
        fMRIdir = fullfile(res_dir_subj,sub_tasks{i});
        if exist(fMRIdir,'dir')
            system(['rm -Rf ' fMRIdir]);
            mkdir(fMRIdir);
        else
            mkdir(fMRIdir);
        end
        
        % Convert dicom to nifti
        runline = sprintf('dcm2nii -o %s %s',fMRIdir,[subject_dir '/' EPIfolders(i).name]);
        [status,res] = system(runline);
        
        % Split 4D volumes into 3Ds
        filein = dir([fMRIdir '/EPI_*']);
        if isempty(filein); filein = dir([fMRIdir '/20*']); end
        runline = sprintf('fslsplit %s %s -t',fullfile(fMRIdir,filein.name),[fMRIdir '/vol']);
        [status,res] = system(runline); delete(fullfile(fMRIdir,filein.name));
        system(['gunzip ' fMRIdir '/*nii.gz']);
        vols = dir([fMRIdir '/vol*']);
        
        %         Store number of volumes per block
        %         Nvols_x_bloq (id_sub,1) = length(vols);
        %         Store acquisition time per block related to the TR (2000ms)
        Tvols_x_bloq{id_sub}(i,1) = 2000*length(vols);
        %
        for id_v = 1:length(vols)
            cont = cont + 1;
            system(['mv ' fMRIdir '/vol' num2str(id_v-1,'%04d') '.nii ' fMRIdir_all '/vol' num2str(i) '_' num2str(cont,'%04d') '.nii']);
        end
        system(['rm -Rf ' fMRIdir]);
    end
    
    % Read triggers file
    triggers_trials_file = dir([subject_dir '/Writing_*/*_trial*.mat']);
    cd ([triggers_trials_file(1).folder]);
    %individual trials infos
    for q = 1:length(triggers_trials_file)
        trials(q,1) = load(triggers_trials_file(q).name);
    end
    %all trials onsets
    triggers_onsets_file = dir([subject_dir '/Writing_*/*alltrials_onsets.mat']);
    all_trials_onsets = load(triggers_onsets_file.name);
    
    triggers.nTrial = all_trials_onsets.all_trials.S.trial_n;
    triggers.onset_image = all_trials_onsets.all_trials.S.onset_image;
    triggers.onset_write = all_trials_onsets.all_trials.S.onset_write;
    triggers.onset_fixation = all_trials_onsets.all_trials.S.onset_fixation;
    triggers.ITI_time = all_trials_onsets.all_trials.S.fixation_time;
    triggers.trial_end = all_trials_onsets.all_trials.S.endtrial_time; % when we press keyboard at trial end
    
    cont_B = 0; % Counter to identify the beginning of each block within the triggers file
    blockvector = [];
    
    for id_t = 1:size(triggers.nTrial,2)
        sample = triggers.nTrial(1,id_t);
        if sample == 1
            cont_B = cont_B + 1;
            blockvector(id_t) = cont_B;
            trigger_init(cont_B) = triggers.onset_image(1,id_t); % condition Onset - Stores per block the starting point
            trigger_rowinit(cont_B) = id_t; % Stores per block the first row within the triggers file
            %             if cont_B > 1
            %                 trigger_finish(cont_B) = triggers.onset_image(1,id_t) + (triggers.trial_end(1,id_t) - triggers.onset_image(1,id_t));%
            %                 trigger_rowfinish(cont_B) = id_t;
            %             end
        else
            blockvector(id_t) = cont_B;
        end
    end
    
    trigger_finish(cont_B) =  triggers.trial_end(1,id_t);
    trigger_rowfinish(cont_B) = id_t;
    trigger_flag_reset = [0 trigger_init] < [0 trigger_finish];
    
    time_x_block = (trigger_finish - trigger_init); % Tiempo que dura cada bloque de acuerdo a los triggers
    excessTime_x_block = (Tvols_x_bloq{id_sub}')/1000 - time_x_block; % Tiempo por encima que tenemos en volumenes
    excessVols_x_block = floor(abs(excessTime_x_block/2)).*sign(excessTime_x_block); % El numero de volumenes positivos se tienen que eliminar del final
    
    % Hay que a??adir al punto de inicio de cada bloque el tiempo que le
    % queda al bloque anterior para terminar la adquisicion del ultimo
    % volumen
    
    % Delete excess volumes
    for i_Bloque = 1:length(EPIfolders)
        if excessVols_x_block(i_Bloque) > 0
            vols = dir([fMRIdir_all '/vol' num2str(i_Bloque) '*']);
            for i_del = 1:excessVols_x_block(i_Bloque)
                delete([fMRIdir_all '/' vols(end-i_del+1).name]);
            end
        end
    end
    
    d = all_trials_onsets.all_trials.S.onset_image - all_trials_onsets.all_trials.S.fixation_time;
    FixOnset = d';
    offset = FixOnset(1);
    FixOnset = (FixOnset - FixOnset(1));
    trial_duration = all_trials_onsets.all_trials.S.endtrial_time - all_trials_onsets.all_trials.S.onset_image;
    trial_duration  = trial_duration';
    
    
    % Trial Types: 1 spanish, 2 greek, 3 polish, 4 signatures
    condition = all_trials_onsets.all_trials.S.trial_cond;
    for id_t = 1:size(triggers.onset_image,2)
        %onset image
        spanish = [];greek = [];polish = [];signature = [];
        cont = 0;
        for id_t = 1:size(triggers.onset_image,2)
            if condition(1,id_t) == 1
                cont = cont + 1;
                spanish(cont) = id_t;
                Onset_spanish_i(cont,1) = (triggers.onset_image(1,id_t) - offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_image,2)
            if condition(1,id_t) == 2
                cont = cont + 1;
                greek(cont) = id_t;
                Onset_greek_i(cont,1) = (triggers.onset_image(1,id_t) - offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_image,2)
            if condition(1,id_t) == 3
                cont = cont + 1;
                polish(cont) = id_t;
                Onset_polish_i(cont,1) = (triggers.onset_image(1,id_t)- offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_image,2)
            if condition(1,id_t) == 4
                cont = cont + 1;
                signature(cont) = id_t;
                Onset_signature_i(cont,1) = (triggers.onset_image(1,id_t)- offset);
            end
        end
        
        % sobrarian estos onsets write....
        %onset write
        cont = 0; spanish_w = [];greek_w = [];polish_w = [];signature_w = [];
        for id_t = 1:size(triggers.onset_write,2)
            if condition(1,id_t) == 1
                cont = cont + 1;
                spanish_w(cont) = id_t;
                Onset_spanish_w(cont,1) = (triggers.onset_write(1,id_t)- offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_write,2)
            if condition(1,id_t) == 2
                cont = cont + 1;
                greek_w(cont) = id_t;
                Onset_greek_w(cont,1) = (triggers.onset_write(1,id_t)- offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_write,2)
            if condition(1,id_t) == 3
                cont = cont + 1;
                polish_w(cont) = id_t;
                Onset_polish_w(cont,1) = (triggers.onset_write(1,id_t)- offset);
            end
        end
        cont = 0;
        for id_t = 1:size(triggers.onset_write,2)
            if condition(1,id_t) == 4
                cont = cont + 1;
                signature_w(cont) = id_t;
                Onset_signature_w(cont,1) = (triggers.onset_write(1,id_t)- offset);
            end
        end
        
        % Trial Types: 1 spanish, 2 greek, 3 polish, 4 signatures
        trial_spanish = condition==1;trial_greek = condition==2;
        trial_polish = condition==3;trial_sign = condition==4;
        
        trial_duration_spa = trial_duration(trial_spanish==1);
        trial_duration_greek = trial_duration(trial_greek==1);
        trial_duration_polish = trial_duration(trial_polish==1);
        trial_duration_sign = trial_duration(trial_sign==1);
        
    end
    
    save(fullfile(fMRIdir_all,'Onsets_write.mat'),'FixOnset','Onset_spanish_i','Onset_greek_i','Onset_polish_i','Onset_signature_i','Onset_spanish_w','Onset_greek_w','Onset_polish_w','Onset_signature_w','trial_duration_spa','trial_duration_greek', 'trial_duration_polish', 'trial_duration_sign');
    cd(fMRIdir_all)
end