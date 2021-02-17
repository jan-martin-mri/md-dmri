elstat_dir = fullfile(pwd, 'Elstat_multiopt');
N = 79;

for n = 3:3
    txt_file = fullfile(elstat_dir, strcat('Grad_dirs_multiopt_', num2str(n), '.txt'));
    matlab_file = fullfile(pwd, strcat('uvec_elstat_', num2str(n), '.m')); 
    
    fid_txt = fopen(txt_file, 'r');
    u_txt = cell([n, 1]);
    count = 1;
    while ~feof(fid_txt)
        u_txt{count} = fgetl(fid_txt);
        count = count + 1;
%         i = regexp(line, '\w', 'match')
    end
    
    
%     while ~feof(fid_initial)
%         line = [fgetl(fid_initial) newline];   % read line
%         for i = 1:length(str_initial_list)
%             if contains(line, str_initial_list{i})
%                 line = strrep(line, str_initial_list{i}, str_final_list{i});
%             end
%         end
%         fprintf(fid_matlab, '%s', line);
%     end
    
    
%     fid_matlab = fopen(text_file_final, 'w');
    fclose(fid_txt);
%     fclose(fid_matlab);

end

