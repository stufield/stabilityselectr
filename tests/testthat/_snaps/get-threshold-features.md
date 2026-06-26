# `get_threshold_features()` returns the correct values

    Code
      print(thresh_feat_tbl)
    Output
      $thresh_0.6
      # A tibble: 20 x 3
         feature MaxSelectProb FDRbound
         <chr>           <dbl>    <dbl>
       1 feat_t          0.95    0.0125
       2 feat_d          0.92    0.025 
       3 feat_m          0.92    0.0375
       4 feat_s          0.92    0.05  
       5 feat_g          0.9     0.0625
       6 feat_c          0.895   0.075 
       7 feat_j          0.885   0.0875
       8 feat_a          0.875   0.1   
       9 feat_l          0.875   0.1125
      10 feat_b          0.865   0.125 
      11 feat_i          0.865   0.1375
      12 feat_r          0.865   0.15  
      13 feat_o          0.86    0.1625
      14 feat_q          0.855   0.175 
      15 feat_f          0.835   0.1875
      16 feat_h          0.825   0.2   
      17 feat_k          0.825   0.2125
      18 feat_e          0.82    0.225 
      19 feat_n          0.81    0.2375
      20 feat_p          0.81    0.25  
      
      $thresh_0.7
      # A tibble: 20 x 3
         feature MaxSelectProb FDRbound
         <chr>           <dbl>    <dbl>
       1 feat_t          0.95   0.00625
       2 feat_d          0.92   0.0125 
       3 feat_m          0.92   0.01875
       4 feat_s          0.92   0.025  
       5 feat_g          0.9    0.03125
       6 feat_c          0.895  0.0375 
       7 feat_j          0.885  0.04375
       8 feat_a          0.875  0.05   
       9 feat_l          0.875  0.05625
      10 feat_b          0.865  0.0625 
      11 feat_i          0.865  0.06875
      12 feat_r          0.865  0.075  
      13 feat_o          0.86   0.08125
      14 feat_q          0.855  0.0875 
      15 feat_f          0.835  0.09375
      16 feat_h          0.825  0.1    
      17 feat_k          0.825  0.1063 
      18 feat_e          0.82   0.1125 
      19 feat_n          0.81   0.1188 
      20 feat_p          0.81   0.125  
      
      $thresh_0.8
      # A tibble: 20 x 3
         feature MaxSelectProb FDRbound
         <chr>           <dbl>    <dbl>
       1 feat_t          0.95  0.004167
       2 feat_d          0.92  0.008333
       3 feat_m          0.92  0.0125  
       4 feat_s          0.92  0.01667 
       5 feat_g          0.9   0.02083 
       6 feat_c          0.895 0.025   
       7 feat_j          0.885 0.02917 
       8 feat_a          0.875 0.03333 
       9 feat_l          0.875 0.0375  
      10 feat_b          0.865 0.04167 
      11 feat_i          0.865 0.04583 
      12 feat_r          0.865 0.05    
      13 feat_o          0.86  0.05417 
      14 feat_q          0.855 0.05833 
      15 feat_f          0.835 0.0625  
      16 feat_h          0.825 0.06667 
      17 feat_k          0.825 0.07083 
      18 feat_e          0.82  0.075   
      19 feat_n          0.81  0.07917 
      20 feat_p          0.81  0.08333 
      
      $thresh_0.9
      # A tibble: 5 x 3
        feature MaxSelectProb FDRbound
        <chr>           <dbl>    <dbl>
      1 feat_t           0.95 0.003125
      2 feat_d           0.92 0.00625 
      3 feat_m           0.92 0.009375
      4 feat_s           0.92 0.0125  
      5 feat_g           0.9  0.01562 
      

# `get_threshold_features()` trips warning when no stable features found

    Code
      print(thresh_feat_tbl)
    Output
      $thresh_0.92
      # A tibble: 4 x 3
        feature MaxSelectProb FDRbound
        <chr>           <dbl>    <dbl>
      1 feat_t           0.95 0.002976
      2 feat_d           0.92 0.005952
      3 feat_m           0.92 0.008929
      4 feat_s           0.92 0.01190 
      
      $thresh_0.97
      # A tibble: 0 x 3
      # i 3 variables: feature <chr>, MaxSelectProb <dbl>, FDRbound <dbl>
      

