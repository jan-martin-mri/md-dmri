function m_bs = jm_bs2m( bootstrap_path, bootstrap_idx )

% Collect all m vectors contained in the directory 'bootstrap_path'
% indicated by 'bootstrap_idx' and concatenate them into a 5D array with
% dimensions x x y x z x n x bs.
% Example: m_bs = jm_bs2m(fn, [1 3 5 6]) collects bootstraps number 1, 3,
% 5, and 6.

m_bs = [];
for bs = bootstrap_idx
    fn = fullfile(bootstrap_path, num2str(bs), 'mfs.mat');
    disp(['Importing ' fn]);
    mfs = mdm_mfs_load(fn);
    m = mfs.m;
    m_bs = cat(5, m_bs, m);
end

end

