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
       1 feat_a  61.50
       2 feat_b  55.08
       3 feat_c  60.34
       4 feat_d  74.63
       5 feat_e  56.81
       6 feat_f  53.00
       7 feat_g  54.04
       8 feat_h  56.48
       9 feat_i  52.46
      10 feat_j  76.23
      11 feat_k  49.62
      12 feat_l  52.37
      13 feat_m  84.85
      14 feat_n  52.56
      15 feat_o  49.20
      16 feat_p  46.06
      17 feat_q  59.02
      18 feat_r  66.58
      19 feat_s  62.3 
      20 feat_t  79.24

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
       1 feat_t          0.93  0.436532       NA
       2 feat_d          0.91  0.408626       NA
       3 feat_q          0.905 0.256633       NA
       4 feat_a          0.9   0.286384       NA
       5 feat_g          0.895 0.222145       NA
       6 feat_j          0.89  0.440141       NA
       7 feat_s          0.89  0.291720       NA
       8 feat_r          0.885 0.342005       NA
       9 feat_m          0.88  0.533218       NA
      10 feat_n          0.87  0.220431       NA
      11 feat_b          0.865 0.245937       NA
      12 feat_e          0.865 0.253644       NA
      13 feat_c          0.85  0.287507       NA
      14 feat_f          0.85  0.224109       NA
      15 feat_h          0.845 0.243903       NA
      16 feat_l          0.845 0.221616       NA
      17 feat_o          0.825 0.198936       NA
      18 feat_p          0.815 0.186961       NA
      19 feat_k          0.81  0.216318       NA
      20 feat_i          0.8   0.235128       NA

---

    Code
      summary(ss, thresh = 0.9)
    Output
      # A tibble: 4 x 4
        feature MaxSelectProb      AUC FDRbound
        <chr>           <dbl>    <dbl>    <dbl>
      1 feat_t          0.93  0.436532 0.003125
      2 feat_d          0.91  0.408626 0.00625 
      3 feat_q          0.905 0.256633 0.009375
      4 feat_a          0.9   0.286384 0.0125  

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
       1 seq.2802.68  0    
       2 seq.9251.29  0.26 
       3 seq.1942.70  0.305
       4 seq.5751.80  0.335
       5 seq.9608.12  0    
       6 seq.3459.49  3.13 
       7 seq.3865.56  0    
       8 seq.3363.21  0    
       9 seq.4487.88  0.13 
      10 seq.5994.84  0.015
      11 seq.9011.72  0.62 
      12 seq.2902.23  6.225
      13 seq.2260.48  2.875
      14 seq.4936.96  0.69 
      15 seq.2277.95  0.915
      16 seq.2953.31  0    
      17 seq.3032.11  2.695
      18 seq.4330.4   0.025
      19 seq.4914.10  0.04 
      20 seq.3896.5   1.345
      21 seq.5002.7   0.05 
      22 seq.3476.4   0    
      23 seq.1130.49  0.065
      24 seq.6356.60  0    
      25 seq.4579.40  0    
      26 seq.8344.24  0    
      27 seq.8441.53  0.395
      28 seq.9360.55  0.08 
      29 seq.7841.8   0    
      30 seq.8142.63  0.03 
      31 seq.4461.56  0.075
      32 seq.9297.97  3.31 
      33 seq.9396.38  0.075
      34 seq.3300.26  0.08 
      35 seq.2772.14  0.135
      36 seq.6615.18  0    
      37 seq.8797.98  0.35 
      38 seq.9879.88  0    
      39 seq.8993.16  0.03 
      40 seq.9373.82  0.185

