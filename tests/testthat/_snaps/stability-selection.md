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
         feat   rowsum
         <chr>   <dbl>
       1 feat_a  60.94
       2 feat_b  51.04
       3 feat_c  61.62
       4 feat_d  76.32
       5 feat_e  50.31
       6 feat_f  52.08
       7 feat_g  58.62
       8 feat_h  51.34
       9 feat_i  51.09
      10 feat_j  75.84
      11 feat_k  54.74
      12 feat_l  52.37
      13 feat_m  84.61
      14 feat_n  48.41
      15 feat_o  49.42
      16 feat_p  46.21
      17 feat_q  54.37
      18 feat_r  66.18
      19 feat_s  64.82
      20 feat_t  80.20

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
       1 feat_t          0.95  0.413037       NA
       2 feat_d          0.92  0.394679       NA
       3 feat_m          0.92  0.374027       NA
       4 feat_s          0.92  0.328683       NA
       5 feat_g          0.9   0.303614       NA
       6 feat_c          0.895 0.328239       NA
       7 feat_j          0.885 0.413677       NA
       8 feat_a          0.875 0.323781       NA
       9 feat_l          0.875 0.269080       NA
      10 feat_b          0.865 0.245967       NA
      11 feat_i          0.865 0.249322       NA
      12 feat_r          0.865 0.347174       NA
      13 feat_o          0.86  0.240215       NA
      14 feat_q          0.855 0.273082       NA
      15 feat_f          0.835 0.251229       NA
      16 feat_h          0.825 0.261538       NA
      17 feat_k          0.825 0.284152       NA
      18 feat_e          0.82  0.254872       NA
      19 feat_n          0.81  0.247171       NA
      20 feat_p          0.81  0.219970       NA

---

    Code
      summary(ss, thresh = 0.9)
    Output
      # A tibble: 5 x 4
        feature MaxSelectProb      AUC FDRbound
        <chr>           <dbl>    <dbl>    <dbl>
      1 feat_t           0.95 0.413037 0.003125
      2 feat_d           0.92 0.394679 0.00625 
      3 feat_m           0.92 0.374027 0.009375
      4 feat_s           0.92 0.328683 0.0125  
      5 feat_g           0.9  0.303614 0.015625

# `stab_sel` S3 print method returns expected known output

    == Stability Selection (Kernel: binomial) ======================================
    * Weakness (alpha)            0.8
    * Weakness Probability (Pw)   0.5
    * Number of Iterations        100
    * Standardized                'Yes'
    * Imputed Outliers            'No'
    * Lambda Max                  0.1471
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
       1 seq.2802.68  0.01 
       2 seq.9251.29  0.145
       3 seq.1942.70  0.08 
       4 seq.5751.80  0.505
       5 seq.9608.12  0    
       6 seq.3459.49  3.525
       7 seq.3865.56  0.225
       8 seq.3363.21  0    
       9 seq.4487.88  0    
      10 seq.5994.84  0    
      11 seq.9011.72  0.42 
      12 seq.2902.23  5.91 
      13 seq.2260.48  3.23 
      14 seq.4936.96  0.88 
      15 seq.2277.95  0.835
      16 seq.2953.31  0.045
      17 seq.3032.11  2.53 
      18 seq.4330.4   0.02 
      19 seq.4914.10  0.105
      20 seq.3896.5   1.28 
      21 seq.5002.7   0.22 
      22 seq.3476.4   0    
      23 seq.1130.49  0.085
      24 seq.6356.60  0    
      25 seq.4579.40  0.005
      26 seq.8344.24  0.065
      27 seq.8441.53  0.055
      28 seq.9360.55  0    
      29 seq.7841.8   0.005
      30 seq.8142.63  0.145
      31 seq.4461.56  0.045
      32 seq.9297.97  3.44 
      33 seq.9396.38  0    
      34 seq.3300.26  0.16 
      35 seq.2772.14  0    
      36 seq.6615.18  0    
      37 seq.8797.98  0.345
      38 seq.9879.88  0    
      39 seq.8993.16  0    
      40 seq.9373.82  0.06 

