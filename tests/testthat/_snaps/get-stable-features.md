# `get_stable_features()` returns correct values at thresh = 0.55

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 20 x 3
         feature MaxSelectProb FDRbound
         <chr>           <dbl>    <dbl>
       1 feat_t          0.95   0.02500
       2 feat_d          0.92   0.05000
       3 feat_m          0.92   0.07500
       4 feat_s          0.92   0.10000
       5 feat_g          0.9    0.1250 
       6 feat_c          0.895  0.1500 
       7 feat_j          0.885  0.1750 
       8 feat_a          0.875  0.2000 
       9 feat_l          0.875  0.2250 
      10 feat_b          0.865  0.2500 
      11 feat_i          0.865  0.2750 
      12 feat_r          0.865  0.3000 
      13 feat_o          0.86   0.3250 
      14 feat_q          0.855  0.3500 
      15 feat_f          0.835  0.3750 
      16 feat_h          0.825  0.4000 
      17 feat_k          0.825  0.4250 
      18 feat_e          0.82   0.4500 
      19 feat_n          0.81   0.4750 
      20 feat_p          0.81   0.5000 

# `get_stable_features()` returns correct values at thresh <= 0.5

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 20 x 3
         feature MaxSelectProb FDRbound
         <chr>           <dbl>    <dbl>
       1 feat_t          0.95        NA
       2 feat_d          0.92        NA
       3 feat_m          0.92        NA
       4 feat_s          0.92        NA
       5 feat_g          0.9         NA
       6 feat_c          0.895       NA
       7 feat_j          0.885       NA
       8 feat_a          0.875       NA
       9 feat_l          0.875       NA
      10 feat_b          0.865       NA
      11 feat_i          0.865       NA
      12 feat_r          0.865       NA
      13 feat_o          0.86        NA
      14 feat_q          0.855       NA
      15 feat_f          0.835       NA
      16 feat_h          0.825       NA
      17 feat_k          0.825       NA
      18 feat_e          0.82        NA
      19 feat_n          0.81        NA
      20 feat_p          0.81        NA

# `get_stable_features` returns correct values at thresh = 0.90

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 5 x 3
        feature MaxSelectProb FDRbound
        <chr>           <dbl>    <dbl>
      1 feat_t           0.95 0.003125
      2 feat_d           0.92 0.00625 
      3 feat_m           0.92 0.009375
      4 feat_s           0.92 0.0125  
      5 feat_g           0.9  0.01562 

# `get_stable_features()` with permutation adds `EmpFDR` column

    Code
      print(s_feat_perm, n = Inf)
    Output
      # A tibble: 10 x 4
         feature MaxSelectProb FDRbound EmpFDR
         <chr>           <dbl>    <dbl>  <dbl>
       1 feat_t          0.95  0.003289   0.14
       2 feat_m          0.925 0.006579   0.34
       3 feat_r          0.92  0.009868   0.35
       4 feat_d          0.91  0.01316    0.45
       5 feat_s          0.91  0.01645    0.45
       6 feat_a          0.905 0.01974    0.54
       7 feat_g          0.905 0.02303    0.54
       8 feat_c          0.89  0.02632    0.63
       9 feat_f          0.88  0.02961    0.74
      10 feat_l          0.88  0.03289    0.74

