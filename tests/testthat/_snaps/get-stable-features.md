# `get_stable_features()` returns correct values at thresh = 0.55

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 20 x 3
         feature MaxSelectProb FDR_bound
         <chr>           <dbl>     <dbl>
       1 feat_t          0.93    0.02500
       2 feat_d          0.91    0.05000
       3 feat_q          0.905   0.07500
       4 feat_a          0.9     0.10000
       5 feat_g          0.895   0.1250 
       6 feat_j          0.89    0.1500 
       7 feat_s          0.89    0.1750 
       8 feat_r          0.885   0.2000 
       9 feat_m          0.88    0.2250 
      10 feat_n          0.87    0.2500 
      11 feat_b          0.865   0.2750 
      12 feat_e          0.865   0.3000 
      13 feat_c          0.85    0.3250 
      14 feat_f          0.85    0.3500 
      15 feat_h          0.845   0.3750 
      16 feat_l          0.845   0.4000 
      17 feat_o          0.825   0.4250 
      18 feat_p          0.815   0.4500 
      19 feat_k          0.81    0.4750 
      20 feat_i          0.8     0.5000 

# `get_stable_features()` returns correct values at thresh <= 0.5

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 20 x 3
         feature MaxSelectProb FDR_bound
         <chr>           <dbl>     <dbl>
       1 feat_t          0.93         NA
       2 feat_d          0.91         NA
       3 feat_q          0.905        NA
       4 feat_a          0.9          NA
       5 feat_g          0.895        NA
       6 feat_j          0.89         NA
       7 feat_s          0.89         NA
       8 feat_r          0.885        NA
       9 feat_m          0.88         NA
      10 feat_n          0.87         NA
      11 feat_b          0.865        NA
      12 feat_e          0.865        NA
      13 feat_c          0.85         NA
      14 feat_f          0.85         NA
      15 feat_h          0.845        NA
      16 feat_l          0.845        NA
      17 feat_o          0.825        NA
      18 feat_p          0.815        NA
      19 feat_k          0.81         NA
      20 feat_i          0.8          NA

# `get_stable_features` returns correct values at thresh = 0.90

    Code
      print(s_feat, n = Inf)
    Output
      # A tibble: 4 x 3
        feature MaxSelectProb FDR_bound
        <chr>           <dbl>     <dbl>
      1 feat_t          0.93   0.003125
      2 feat_d          0.91   0.00625 
      3 feat_q          0.905  0.009375
      4 feat_a          0.9    0.0125  

# `get_stable_features()` with permutation adds `EmpFDR` column

    Code
      print(s_feat_perm, n = Inf)
    Output
      # A tibble: 10 x 4
         feature MaxSelectProb FDR_bound EmpFDR
         <chr>           <dbl>     <dbl>  <dbl>
       1 feat_t          0.935  0.003289   0.28
       2 feat_d          0.915  0.006579   0.38
       3 feat_a          0.905  0.009868   0.44
       4 feat_j          0.905  0.01316    0.44
       5 feat_q          0.905  0.01645    0.44
       6 feat_g          0.9    0.01974    0.48
       7 feat_r          0.9    0.02303    0.48
       8 feat_s          0.89   0.02632    0.62
       9 feat_e          0.885  0.02961    0.71
      10 feat_m          0.885  0.03289    0.71

# `get_stable_features()` accepts a vector of thresholds

    Code
      print(thresh_feat_tbl)
    Output
      $thresh_0.6
      # A tibble: 20 x 3
         feature MaxSelectProb FDR_bound
         <chr>           <dbl>     <dbl>
       1 feat_t          0.93     0.0125
       2 feat_d          0.91     0.025 
       3 feat_q          0.905    0.0375
       4 feat_a          0.9      0.05  
       5 feat_g          0.895    0.0625
       6 feat_j          0.89     0.075 
       7 feat_s          0.89     0.0875
       8 feat_r          0.885    0.1   
       9 feat_m          0.88     0.1125
      10 feat_n          0.87     0.125 
      11 feat_b          0.865    0.1375
      12 feat_e          0.865    0.15  
      13 feat_c          0.85     0.1625
      14 feat_f          0.85     0.175 
      15 feat_h          0.845    0.1875
      16 feat_l          0.845    0.2   
      17 feat_o          0.825    0.2125
      18 feat_p          0.815    0.225 
      19 feat_k          0.81     0.2375
      20 feat_i          0.8      0.25  
      
      $thresh_0.7
      # A tibble: 20 x 3
         feature MaxSelectProb FDR_bound
         <chr>           <dbl>     <dbl>
       1 feat_t          0.93    0.00625
       2 feat_d          0.91    0.0125 
       3 feat_q          0.905   0.01875
       4 feat_a          0.9     0.025  
       5 feat_g          0.895   0.03125
       6 feat_j          0.89    0.0375 
       7 feat_s          0.89    0.04375
       8 feat_r          0.885   0.05   
       9 feat_m          0.88    0.05625
      10 feat_n          0.87    0.0625 
      11 feat_b          0.865   0.06875
      12 feat_e          0.865   0.075  
      13 feat_c          0.85    0.08125
      14 feat_f          0.85    0.0875 
      15 feat_h          0.845   0.09375
      16 feat_l          0.845   0.1    
      17 feat_o          0.825   0.1063 
      18 feat_p          0.815   0.1125 
      19 feat_k          0.81    0.1188 
      20 feat_i          0.8     0.125  
      
      $thresh_0.8
      # A tibble: 20 x 3
         feature MaxSelectProb FDR_bound
         <chr>           <dbl>     <dbl>
       1 feat_t          0.93   0.004167
       2 feat_d          0.91   0.008333
       3 feat_q          0.905  0.0125  
       4 feat_a          0.9    0.01667 
       5 feat_g          0.895  0.02083 
       6 feat_j          0.89   0.025   
       7 feat_s          0.89   0.02917 
       8 feat_r          0.885  0.03333 
       9 feat_m          0.88   0.0375  
      10 feat_n          0.87   0.04167 
      11 feat_b          0.865  0.04583 
      12 feat_e          0.865  0.05    
      13 feat_c          0.85   0.05417 
      14 feat_f          0.85   0.05833 
      15 feat_h          0.845  0.0625  
      16 feat_l          0.845  0.06667 
      17 feat_o          0.825  0.07083 
      18 feat_p          0.815  0.075   
      19 feat_k          0.81   0.07917 
      20 feat_i          0.8    0.08333 
      
      $thresh_0.9
      # A tibble: 4 x 3
        feature MaxSelectProb FDR_bound
        <chr>           <dbl>     <dbl>
      1 feat_t          0.93   0.003125
      2 feat_d          0.91   0.00625 
      3 feat_q          0.905  0.009375
      4 feat_a          0.9    0.0125  
      

# `get_stable_features()` trips warning when no stable features found

    Code
      print(thresh_feat_tbl)
    Output
      $thresh_0.92
      # A tibble: 1 x 3
        feature MaxSelectProb FDR_bound
        <chr>           <dbl>     <dbl>
      1 feat_t           0.93  0.002976
      
      $thresh_0.97
      # A tibble: 0 x 3
      # i 3 variables: feature <chr>, MaxSelectProb <dbl>, FDR_bound <dbl>
      

