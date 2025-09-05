srcDir  = '\\home.ansatt.ntnu.no\emiliajm\Documents\MATLAB\huc_slow_osc';     % where scatter_bin*.fig live
outDir  = '\\home.ansatt.ntnu.no\emiliajm\Documents\MATLAB\huc_slow_osc';
pattern = 'scatter_*.fig';            % matches your list
dpi     = 300;                        % 150–300 is a good range
bg      = 'none';                     % or 'white' if you prefer

batch_flatten_scatter_to_svg(srcDir, outDir, pattern, dpi, bg);
