
% remove some movement regressors
% add time per trial
% check individual movement values
% bajar smooth a 6?

clear; close all; clc;
fs=filesep;
% 
% source_dir = '/media/lina/Data_Disk/Habits/Early_fmri_data/post/fmri';
% struct_dir = '/media/lina/Data_Disk/Habits/Early_fmri_data/post/structural';
% res_dir =    '/media/lina/Data_Disk/Habits/Early_fmri_data/post/fMRI_SPM_writing';
% subjs = dir([source_dir '/Early_0*']);
% cd(source_dir);
source_dir = '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_habits';
struct_dir = '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI/Structural_habits';
res_dir =    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing';
subjs = dir([source_dir '/Habits_05*']);
cd(source_dir);

if ~exist(res_dir,'dir'); mkdir(res_dir); end
fslpreline = 'export FSLDIR=/usr/local/fsl ; source $FSLDIR/etc/fslconf/fsl.sh ; $FSLDIR/bin/';

prefix_dcm2nii = 'vol';

% Fieldmap params (same for all subjects)
param.esp = 0.53;
param.asset = 0.5;

PATHS.spm12_path = '/home/lina/Documents/MATLAB/spm12';
rmpath(genpath(PATHS.spm12_path));
addpath(genpath(PATHS.spm12_path));

    for id_s = 1:length(subjs)        
        
        do_fmap = 1;
        subjname = subjs(id_s).name;
        subject_dir = fullfile(res_dir,subjname);
        cd(fullfile(res_dir,subjname));
        
        % Preproc T1
        t1folderin = fullfile(struct_dir,subjname);
        t1dcmfolder = dir([t1folderin, fs, '3D*']);
        T1dir = fullfile(t1folderin,'T1');
        if ~exist(T1dir,'dir')
            mkdir(T1dir)
            
            %dcm2nii
            t1folderin_dcm = fullfile(t1folderin,t1dcmfolder.name);
            runline_dcm2nii = sprintf('dcm2nii %s %s', t1folderin_dcm , t1folderin_dcm );
            system(runline_dcm2nii);
            %gunzip decompress
            decompress = dir(fullfile(t1folderin,t1dcmfolder.name, '20*.nii.gz'));
            decompress = decompress(1);
            decompress_file = fullfile(t1folderin,t1dcmfolder.name,decompress.name);
            T1file = sprintf('gunzip %s', decompress_file);
            system(T1file);
            
            %move to T1 folder
            t1_decompress1 = dir(fullfile(t1folderin,t1dcmfolder.name, '20*.nii*'));
            t1_decompress1 = t1_decompress1(1);
            t1_decompress2 = fullfile(t1folderin,t1dcmfolder.name,t1_decompress1.name);
            t1destiny = fullfile(T1dir,'T1.nii');
            runline = sprintf('cp %s %s', t1_decompress2, t1destiny);
            system(runline);
        end
         
        fmrifolderin = fullfile(res_dir,subjname,'fMRI_all_firma');
        fmrifolderinNoFix = fullfile(res_dir,subjname,'fMRI_all_StatNoScrubbing_NoFix');
        if ~exist(fmrifolderinNoFix,'dir'); mkdir(fmrifolderinNoFix); end
        
        % Run new_segment + DARTEL
        if ~exist(fullfile(T1dir,'rc3T1.nii'),'file')
            clear matlabbatch; spm_jobman('initcfg');
            matlabbatch{1}.spm.spatial.preproc.channel.vols = {[T1dir, '/T1.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[PATHS.spm12_path '/tpm/TPM.nii,6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
            spm('defaults', 'fMRI');
            spm_jobman('run', matlabbatch);
        end
        
        % Fieldmap preparation
        clc;
        fprintf('........................\n');
        fprintf('Fieldmap preparation.\n');
        
        subjdir = fullfile(source_dir,subjname);
        fmap_folders = dir([subjdir, fs, 'GRE_*']);
        fmap_dir = fullfile(subject_dir,'fmap');
        if ~exist(fmap_dir,'dir'); mkdir(fmap_dir); end
        
        if ~isempty(fmap_folders)
            
            % Convert to nifti magnitud image and fieldmap
            fmap_dir_dcm = [subjdir, fs, fmap_folders(end).name];
            runline = ['dcm2nii -o ' fmap_dir ' ' fmap_dir_dcm];
            [status,res] = system(runline);
            
            prefixfmap = '20';
            if isempty(dir([fmap_dir, '/20*.ni*']))
                prefixfmap = 'GRE';
            end
            [~,~] = system(['gzip ' fmap_dir, fs, prefixfmap '*.nii']);
            [status,res] = system(['mv ' fmap_dir, fs, prefixfmap '*nii.gz ' fmap_dir, fs, 'fmap.nii.gz']);
            
            % depende de fecha, magnitud se recoge antes del fmap. en
            % tal caso cambias end por 1 abajo
            mag_dir_dcm = [subjdir, fs, fmap_folders(end).name];
            runline = ['dcm2nii -o ' fmap_dir ' ' mag_dir_dcm];
            [status,res] = system(runline);
            [~,~] = system(['gzip ' fmap_dir, fs, prefixfmap '*nii']);
            nfmapfiles = dir([fmap_dir, fs, prefixfmap '*nii.gz']);
            if length(nfmapfiles) > 1
                delete([fmap_dir, fs, nfmapfiles(1).name]);
            end
            pause(2)
            [status,res] = system(['mv ' fmap_dir, fs, prefixfmap '*nii.gz ' fmap_dir, fs, 'mag.nii.gz']);
            runline = sprintf('fslmaths %s -Tmean %s',[fmap_dir, fs, 'mag.nii.gz'],[fmap_dir, fs, 'mag.nii.gz']);
            [status,res] = system(runline);
            
            % Bet Magnitude Image
            fileIn = fullfile(fmap_dir,'mag.nii.gz');
            fileOut = fullfile(fmap_dir,'mag_brain.nii.gz');
            sentence = sprintf('bet %s %s -f 0.3 -m',fileIn,fileOut);
            [status,result] = system(sentence);
            % Brain mask fieldmap
            fileIn = fullfile(fmap_dir,'fmap.nii.gz');
            fileMas = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            fileOut = fullfile(fmap_dir,'fmap.nii.gz');
            sentence = sprintf('fslmaths %s -mas %s %s',fileIn,fileMas,fileOut);
            [status,result] = system(sentence);
            
            %----------------------------------------------------------------------
            % Prepare fieldmap
            fileIn = fullfile(fmap_dir,'fmap.nii.gz');
            fileMagRes = fullfile(fmap_dir,'mag_brain.nii.gz');
            fileOut = fullfile(fmap_dir,'fmap_rads.nii.gz');
            sentence = sprintf('fsl_prepare_fieldmap SIEMENS %s %s %s %f',fileIn,fileMagRes,fileOut,2.46);
            [status,result] = system(sentence);
            if status; disp('Fallo fsl_prepare_fieldmap'); error(result); end
            pause(3);
            
            % Reslice fmap into fMRI dimensions
            fileRef = fullfile(fmrifolderin,'vol1_0001.nii');
            fileIn = fullfile(fmap_dir,'fmap_rads.nii.gz');
            fileInRes = fullfile(fmap_dir,'fmap_rads.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            pause(3);
            fileIn = fullfile(fmap_dir,'mag_brain.nii.gz');
            fileInRes = fullfile(fmap_dir,'mag_brain.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            pause(3);
            fileIn = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            fileInRes = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            sentence = sprintf('WarpImageMultiTransform 3 %s %s -R %s --tightest-bounding-box --reslice-by-header',fileIn,fileInRes,fileRef);
            [status,result] = system(sentence);
            
            % Delete residual files
            file1 = fullfile(fmap_dir,'mag.nii.gz');
            delete(file1);
        else
            do_fmap = 0;
        end
        
        %--------------------------------------------------------------------------
        % Prepare files and folders in fMRI data
        fprintf('........................\n');
        fprintf('Preparing paths, and running slice timing and realignment.\n');
        clear matlabbatch; spm_jobman('initcfg');
        epivol = dir([fmrifolderin, fs, 'vol*']);
        
        paths_vols = [];
        for id_epi = 1:length(epivol)
            paths_vols{id_epi,1} = [fmrifolderin, fs, epivol(id_epi).name ',1'];
        end
        param.NSlices = 35;
        param.TR = 2000;
        
        matlabbatch{1}.spm.temporal.st.scans = {paths_vols}';
        matlabbatch{1}.spm.temporal.st.nslices = param.NSlices;
        matlabbatch{1}.spm.temporal.st.tr = param.TR/1000;
        matlabbatch{1}.spm.temporal.st.ta = (param.TR/1000) - ((param.TR/1000)/(param.NSlices));
        matlabbatch{1}.spm.temporal.st.so = [1:2:param.NSlices 2:2:param.NSlices];
        matlabbatch{1}.spm.temporal.st.refslice = 1;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear paths_vols;
        clear matlabbatch;
        spm_jobman('initcfg');
        
        volumes1 = dir([fmrifolderin, fs, 'avol1*']);
        for mv = 1:length(volumes1)
            paths_vols1{mv,1} = [fmrifolderin, fs, volumes1(mv).name ',1'];
        end
        
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {paths_vols1}';
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear paths_vols1;
        
        %------------------------------------------------------------------
        % Fieldmap correction for the EPI. Nota, dcm2nii cuando convierte mis
        % imagenes a nifti hace que empiecen por 2017, en tu caso no es asi.
        % Asegurate de como es, y cambia el variable de abajo por el prefijo
        % que le pongan
        fprintf('........................\n');
        fprintf('Running fieldmap correction.\n');
        subjdir = fullfile(source_dir,subjname);
        fmap_dir = fullfile(subject_dir,'fmap');
        
        if do_fmap
            %------------------------------------------------------------------
            % Unwarp EPI data
            fileOut = fullfile(fmap_dir,'merged4Depi.nii.gz');
            sentence = sprintf('%sfslmerge -t %s %s%sravol*',fslpreline,fileOut,fmrifolderin,'/');
            [status,result] = system(sentence);
            pause(3);
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            filefmap = fullfile(fmap_dir,'fmap_rads.nii.gz');
            fileout = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileMask = fullfile(fmap_dir,'mag_brain_mask.nii.gz');
            
            % Dwell time asked here is the time between even and odd echoes.
            % Echo spacing = 1 / [(0019,1028) * (0051,100b component #1)]
            % Where (0019,1028) is the bandwidth per pixel phase encode
            % and (0051,100b component #1) is the acquisition matrix
            % This formula works regardless of the grappa factor, and should
            % match the "ESpacing" in the sequence tab. If we however select
            % the Esp from the sequence tab, we must keep in consideration that
            % if we use GRAPPA we must divide the Esp by the GRAPPA factor.
            runline = sprintf('%sfugue -i %s --dwell=%f --loadfmap=%s -u %s --unwarpdir=y- --mask=%s -s 0.5',fslpreline,fileIn,(param.esp*param.asset)/1000,filefmap,fileout,fileMask);
            [status,result] = system(runline);
            delete(fileIn);
            pause(2);
            %------------------------------------------------------------------
            % Extract mean dewarped volume as a reference
            fileIn = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileOut = fullfile(fmap_dir,'exf_brain.nii.gz');
            sentence = sprintf('%sfslmaths %s -Tmean %s',fslpreline,fileIn,fileOut);
            [status,result] = system(sentence);
            pause(2);
            %------------------------------------------------------
            % Get Mean fMRI, brain extract, then convert to nifti fieldmap,
            % brain extract, coregister to fMRI and reslice, prepare fieldmap
            % with fsl and apply fieldmap correction
            fileIn = fullfile(fmap_dir,'exf_brain.nii.gz');
            fileOut = fullfile(fmrifolderin,'EPI_mean.nii.gz');
            system(['mv ' fileIn ' ' fileOut]);
            system(['gunzip ' fileOut]);
            pause(3);
            
            fileIn = fullfile(fmap_dir,'merged4Depi_dewarped.nii.gz');
            fileOut = fullfile(fmrifolderin,'dravol');
            sentence = sprintf('fslsplit %s %s/dravol -t',fileIn,fmrifolderin);
            [status,result] = system(sentence);
            system(['gunzip ' fmrifolderin '/dravol*']);
            pause(3);
            delete(fileIn);
        else
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            sentence = sprintf('fslmerge -t %s %s/ravol*',fileIn,fmrifolderin);
            [status,result] = system(sentence);
            
            fileOut = fullfile(fmap_dir,'exf_brain.nii.gz');
            sentence = sprintf('fslmaths %s -Tmean %s',fileIn,fileOut);
            [status,result] = system(sentence);
            
            fileIn = fullfile(fmap_dir,'exf_brain.nii.gz');
            fileOut = fullfile(fmrifolderin,'EPI_mean.nii.gz');
            system(['mv ' fileIn ' ' fileOut]);
            system(['gunzip ' fileOut]);
            
            fileIn = fullfile(fmap_dir,'merged4Depi.nii.gz');
            fileOut = fullfile(fmrifolderin,'dravol');
            sentence = sprintf('fslsplit %s %s/dravol -t',fileIn,fmrifolderin);
            [status,result] = system(sentence);
            system(['gunzip ' fmrifolderin '/dravol*']);
            pause(5);
            delete(fileIn);
            
        end
        %------------------------------------------------------------------
        % Estimate coregister to T1 high resolution image, not reslice
        Ref_file = fullfile(T1dir,'T1.nii');
        Source_file = fullfile(fmrifolderin,'EPI_mean.nii');
        Def_file = fullfile(T1dir,'y_T1.nii');
        
        volumes = dir([fmrifolderin, fs, 'dravol*']);
        %  clear paths_vols;
          paths_vols = [];
        for mv = 1:length(volumes)
            paths_vols{mv,1} = [fmrifolderin, fs, volumes(mv).name ',1'];
        end
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[Ref_file ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[Source_file ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = paths_vols;
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        
        % normalization
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {Def_file};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = paths_vols;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
            NaN NaN NaN];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        spm('defaults', 'FMRI');
        
        spm_jobman('run', matlabbatch);
        
        fileIn = fullfile(fmrifolderin, fs,'ravol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin, fs,'dravol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin, fs,'avol*.nii'); system(['rm ' fileIn]);
        fileIn = fullfile(fmrifolderin, fs,'wdravol*.nii'); system(['rm ' fileIn]);
        
        clear paths_vols; volumes = dir([fmrifolderin, fs, 'swdravol*']);
        for mv = 1:length(volumes)
            paths_vols{mv,1} = [fmrifolderin, fs, volumes(mv).name ',1'];
        end
        realign_files = dir([fmrifolderin, fs, 'rp*txt']);
        
        % ADD ALSO DERIVATIVES; AND QUADRATIC FORMS
        for mv = 1:length(realign_files)
            tmp = importdata([fmrifolderin, fs, realign_files(mv).name]);
            if mv == 1
                tmp2 = [zeros(1,6);diff(tmp)];
                R = [tmp tmp2 tmp.^2];
            else
                tmp2 = [zeros(1,6);diff(tmp)];
                R(end+1:end+size(tmp,1),:) = [tmp tmp2 tmp.^2];
            end
        end
        
%         d_circ = 2*pi*50 ; % Head is a sphere of radius equal to 50mm
%         Angle2Dist = d_circ * R(:,4:6) ./ pi;
%         diffR = [zeros(1,6);diff([R(:,1:3) Angle2Dist])];
%         FD = sum(abs(diffR),2);
%         scrubbing = FD < 0.9;
%         scrubbing = scrubbing & [scrubbing(2:end);1] & [1;scrubbing(1:end-1)] & [1;1;scrubbing(1:end-2)];
%         scrubbing = ~scrubbing;
%         scrubbing_idx = find(scrubbing);
%         scrubbing_mat = zeros(size(scrubbing,1),sum(scrubbing));
%         for i_scrub = 1:size(scrubbing_mat,2)
%             scrubbing_mat(scrubbing_idx(i_scrub),i_scrub) = 1;
%         end
%         R = [R scrubbing_mat];
        cd(fmrifolderin);
        rp_file = load(fullfile(fmrifolderin, fs,'rp_avol1_0001.txt'));
        figure (1)
        plot(rp_file(:,1:3))
        saveas(gcf,'movement_trans_fig1.jpg')
        figure (2)
        plot(rp_file(:,4:6))
        saveas(gcf,'movement_rot_fig2.jpg')
        
        save([fmrifolderin, fs, 'motionres.mat'],'R');
        
        load(fullfile(fmrifolderin, fs,'Onsets_write.mat'));
        
        clear matlabbatch; spm_jobman('initcfg');
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fmrifolderinNoFix};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = paths_vols;
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Onset_spanish_i';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = Onset_spanish_i;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Onset_greek_i';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = Onset_greek_i;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Onset_polish_i';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = Onset_polish_i;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Onset_signature_i';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = Onset_signature_i;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
        %
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[fmrifolderin, fs, 'motionres.mat']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'spanish';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1	0	0	0	0	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'greek';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0	0	1	0	0	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'polish';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0	0	0	0	1	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'signature';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0	0	0	0	0	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'signature > spanish';
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1	0	0	0	0	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'spanish > signature';
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1	0	0	0	0	0	-1	0];
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'spanish > polish';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [1	0	0	0	-1	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'spanish > greek';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1	0	-1	0	0	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'signature > polish';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0	0	0	0	-1	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'signature > greek';
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0	0	-1	0	0	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'habit > GDS';
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [1	0	-1	0	-1	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'GDS > habit';
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [-1	0	1	0	1	0	-1	0];
        matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'habit';
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [1	0	0	0	0	0	1	0];
        matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'GDS';
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [0	0	1	0	1	0	0	0];
        matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'greek > signature';
        matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [0	0	1	0	0	0	-1	0];
        matlabbatch{3}.spm.stats.con.delete = 0;
                
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);        
%         fileIn = fullfile(fmrifolderin, fs,'Res*.nii'); system(['rm ' fileIn]);
        
        % Compute FD, DVARS, and motion representative values to reject outlier
        % subjects
%         R(:,4:6) = R(:,4:6)*180/pi;
%         Motion(id_s).MaxDisp = max(abs(R),[],1);
%         Motion(id_s).MaxRelDisp = max(abs(diff(R,1,1)),[],1);
%         Motion(id_s).MeanRMS = rms(R,1);
%         Motion(id_s).id = subjname;
%         clear R path_vo* Ons* volumes*
%         save([fmrifolderin, fs, 'All_motionparam_write.mat'],'Motion');
        clear matlabbatch;
    end
% end