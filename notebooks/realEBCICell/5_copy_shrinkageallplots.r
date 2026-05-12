# done on server bc windows limits path file length....
# copy over the all.pdf and all.png pics over from previously saved EBCI/shrinkage/ folder
# this is to more easily copy just these plots to overleaf, without uploading too many large things
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset

dir.create(sprintf('%s/final/', plot_path))

# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== select only some prev saved files =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////
# copy over the all.pdf and all.png pics over from previously saved EBCI/shrinkage/ folder
# this is to more easily copy just these plots to overleaf, without uploading too many large things
assertthat::assert_that(dir.exists(sprintf('%s/final/', plot_path)))
dir.create(sprintf('%s/final/selAllPlots/', plot_path))




files_to_copy = list.files(path = sprintf('%s/shrinkage/', plot_path) , 
                           pattern = 'all', recursive=TRUE, full.names = TRUE)
dir_to_create = sub(pattern = sprintf('%s/shrinkage/', plot_path), replacement = '', x = files_to_copy) # just the folder names
dir_to_create = sub(pattern = 'all.((png)|(pdf))', replacement = '', x = dir_to_create) # just the folder names

for(i in 1:length(files_to_copy)) {
  cur_file   = files_to_copy[i]
  cur_folder = dir_to_create[i]
  print(sprintf('%s/final/selAllPlots/%s', plot_path, cur_folder))
  # dir.create(sprintf('%s/final/selAllPlots/%s', plot_path, cur_folder), recursive=TRUE)
  
  
  
  
  file.copy(from = cur_file, 
            to = sprintf('%s/final/selAllPlots/%s%s', 
                         plot_path,                                                 # original plot path
                         gsub(pattern = '/', replacement = '_', x = cur_folder), # no sub directories, just name w _ between
                         substr(cur_file, nchar(cur_file)-6, nchar(cur_file))),     # preserve pdf or png
            overwrite = TRUE)
  
  
}





