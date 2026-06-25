# `stability_selection()` generates correct object

    Code
      lapply(ss, class)
    Output
      $stabpath_matrix
      [1] "matrix" "array" 
      
      $lambda
      [1] "numeric"
      
      $alpha
      [1] "numeric"
      
      $Pw
      [1] "numeric"
      
      $kernel
      [1] "character"
      
      $n_iter
      [1] "numeric"
      
      $standardize
      [1] "logical"
      
      $impute_outliers
      [1] "logical"
      
      $lambda_min_ratio
      [1] "numeric"
      
      $perm_data
      [1] "logical"
      
      $permpath_list
      [1] "NULL"
      
      $perm_lambda
      [1] "numeric"
      
      $permpath_max
      [1] "tbl_df"     "tbl"        "data.frame"
      
      $beta
      [1] "dgCMatrix"
      attr(,"package")
      [1] "Matrix"
      
      $r_seed
      [1] "numeric"
      

# `stability_selection()` generates the correct values

    Code
      print(tibble::enframe(rowSums(ss$stabpath_matrix), "feat", "rowsum"), n = Inf)
    Output
      # A tibble: 20 x 2
         feat    rowsum
         <chr>    <dbl>
       1 feat_a 64.96  
       2 feat_b 52.73  
       3 feat_c 60.875 
       4 feat_d 77.0550
       5 feat_e 54.44  
       6 feat_f 59.56  
       7 feat_g 53.285 
       8 feat_h 53.49  
       9 feat_i 49.225 
      10 feat_j 75.985 
      11 feat_k 55.3   
      12 feat_l 63.265 
      13 feat_m 86.2400
      14 feat_n 46.3   
      15 feat_o 47.75  
      16 feat_p 47.83  
      17 feat_q 56.51  
      18 feat_r 65.8   
      19 feat_s 66.645 
      20 feat_t 79.1750

# the S3 `summary` method generates correct output

    Code
      summary(ss, thresh = 0.49, warn = TRUE)
    Condition
      Warning:
      FDR upper bound not defined for `thresh <= 0.5`
    Output
      # A tibble: 20 x 4
         feature MaxSelectProb      AUC FDRbound
         <chr>           <dbl>    <dbl>    <dbl>
       1 feat_d          0.945 0.422319       NA
       2 feat_t          0.935 0.415582       NA
       3 feat_a          0.915 0.342093       NA
       4 feat_j          0.915 0.430052       NA
       5 feat_s          0.905 0.343445       NA
       6 feat_m          0.9   0.410396       NA
       7 feat_l          0.89  0.326203       NA
       8 feat_f          0.885 0.299708       NA
       9 feat_q          0.88  0.290240       NA
      10 feat_g          0.87  0.269070       NA
      11 feat_e          0.86  0.283133       NA
      12 feat_n          0.86  0.226053       NA
      13 feat_r          0.86  0.352070       NA
      14 feat_c          0.855 0.326408       NA
      15 feat_k          0.85  0.292292       NA
      16 feat_b          0.845 0.275154       NA
      17 feat_o          0.84  0.243122       NA
      18 feat_h          0.835 0.274262       NA
      19 feat_p          0.83  0.236150       NA
      20 feat_i          0.785 0.241943       NA

---

    Code
      summary(ss, thresh = 0.9)
    Output
      # A tibble: 6 x 4
        feature MaxSelectProb      AUC FDRbound
        <chr>           <dbl>    <dbl>    <dbl>
      1 feat_d          0.945 0.422319 0.003125
      2 feat_t          0.935 0.415582 0.00625 
      3 feat_a          0.915 0.342093 0.009375
      4 feat_j          0.915 0.430052 0.0125  
      5 feat_s          0.905 0.343445 0.015625
      6 feat_m          0.9   0.410396 0.01875 

# `stab_sel` S3 print method returns expected known output

    == Stability Selection (Kernel: l1-logistic) ===================================
    * Weakness (alpha)            0.8
    * Weakness Probability (Pw)   0.5
    * Number of Iterations        100
    * Standardized                'Yes'
    * Imputed Outliers            'No'
    * Lambda Max                  0.144
    * Lambda Min Ratio            0.1
    * Permuted Data               'No'
    * Random Seed                 101
    ================================================================================

# `stability_selection()` generates expected values for the Cox kernel

    Code
      print(tibble::enframe(rowSums(ss_cox$stabpath_matrix), "feat", "rowsum"), n = Inf)
    Output
      # A tibble: 40 x 2
         feat        rowsum
         <chr>        <dbl>
       1 seq.2802.68  0.1  
       2 seq.9251.29  0.38 
       3 seq.1942.70  0.34 
       4 seq.5751.80  0.59 
       5 seq.9608.12  0.06 
       6 seq.3459.49  3.1  
       7 seq.3865.56  0.015
       8 seq.3363.21  0.155
       9 seq.4487.88  0.135
      10 seq.5994.84  0    
      11 seq.9011.72  0.24 
      12 seq.2902.23  6.11 
      13 seq.2260.48  3.565
      14 seq.4936.96  1.58 
      15 seq.2277.95  0.965
      16 seq.2953.31  0.02 
      17 seq.3032.11  2.94 
      18 seq.4330.4   0    
      19 seq.4914.10  0.055
      20 seq.3896.5   1.77 
      21 seq.5002.7   0.125
      22 seq.3476.4   0    
      23 seq.1130.49  0.025
      24 seq.6356.60  0.005
      25 seq.4579.40  0.16 
      26 seq.8344.24  0.02 
      27 seq.8441.53  0.29 
      28 seq.9360.55  0.015
      29 seq.7841.8   0.015
      30 seq.8142.63  0    
      31 seq.4461.56  0.01 
      32 seq.9297.97  3.83 
      33 seq.9396.38  0    
      34 seq.3300.26  0.05 
      35 seq.2772.14  0.095
      36 seq.6615.18  0    
      37 seq.8797.98  0.29 
      38 seq.9879.88  0    
      39 seq.8993.16  0.005
      40 seq.9373.82  0    

