# `calc_emp_fdr_breaks()` function generates `breaks` entry

    Code
      print(empFDR$breaks, n = Inf)
    Output
      # A tibble: 5 x 4
        FDR_breaks MeanFPs n_selected piThresh
             <dbl>   <dbl>      <int>    <dbl>
      1        0.5     5.8          7      0.9
      2        1       5.8          7      0.9
      3        2       5.8          7      0.9
      4        3       5.8          7      0.9
      5        5       5.8          7      0.9

