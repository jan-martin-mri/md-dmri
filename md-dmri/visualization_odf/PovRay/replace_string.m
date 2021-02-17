function replace_string(text_file_initial, text_file_final, str_initial_list, str_final_list)

fid_initial = fopen(text_file_initial, 'r');
fid_final = fopen(text_file_final, 'w');

while ~feof(fid_initial)
  line = [fgetl(fid_initial) newline];   % read line
  for i = 1:length(str_initial_list)
      if contains(line, str_initial_list{i})
          line = strrep(line, str_initial_list{i}, str_final_list{i});
      end
  end
  fprintf(fid_final, '%s', line);
end

fclose(fid_initial);
fclose(fid_final);

end