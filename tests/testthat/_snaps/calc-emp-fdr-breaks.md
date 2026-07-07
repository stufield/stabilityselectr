# `calc_emp_fdr_breaks()` function generates `breaks` entry

    Code
      print(empFDR$breaks, n = Inf)
    Output
      # A tibble: 5 x 4
        FDR_breaks MeanFPs n_selected piThresh
             <dbl>   <dbl>      <int>    <dbl>
      1        0.5     4.8          7      0.9
      2        1       4.8          7      0.9
      3        2       4.8          7      0.9
      4        3       4.8          7      0.9
      5        5      19.2         20      0.8

