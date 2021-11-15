directory = '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/TTest_writing_withoutsubjectssinactivation/Ttest_0013_habit';
matlabbatch=[];

matlabbatch{1}.spm.stats.factorial_design.dir = {directory};
if ~(exist(directory))
    mdcom=['!mkdir ' directory];
    eval(mdcom);
end
%%
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {
     '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_006/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'    
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_007/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
%      '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_010/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1' 
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_025/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_029/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'    
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_035/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
%mov    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_036/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'     
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_038/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_041/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_042/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_043/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_044/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_046/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'    
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_047/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_048/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_049/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_050/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_051/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
%                                                           
                                                           };
%%
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {
%     '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_001/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'    
%     '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_003/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_004/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'     
 %   '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_005/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
   %mov '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_008/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
   % '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_009/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_011/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
   % '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_012/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_013/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_014/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_015/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
   % '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_016/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_017/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_018/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_019/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'    
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_020/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_022/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_023/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
  %  '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_024/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_026/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_028/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_029/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    %'/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_031/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_032/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
   % '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_033/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
 %mov   '/media/lina/Data_Disk/Habits/habit_fmri_data/fMRI_SPM_writing/Habits_034/fMRI_all_StatNoScrubbing_NoFix/con_0013.nii,1'
    

                                                           };
%%
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PD > HC (habit)';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'HC > PD (habit)';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'HC (Habit)';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'PD (Habit)';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 0];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm('defaults','FMRI');
spm_jobman('interactive', matlabbatch); %serial
