# the object returned by `progeny_cluster()` is correct

    Code
      lapply(pclust, class)
    Output
      $scores
      [1] "matrix" "array" 
      
      $mean_scores
      [1] "numeric"
      
      $ci95_scores
      [1] "matrix" "array" 
      
      $random_scores
      [1] "matrix" "array" 
      
      $mean_random_scores
      [1] "numeric"
      
      $D_max
      [1] "numeric"
      
      $D_gap
      [1] "numeric"
      
      $clust_iter
      [1] "integer"
      
      $repeats
      [1] "integer"
      
      $n_iter
      [1] "integer"
      
      $size
      [1] "numeric"
      
      $call
      [1] "call"
      

# the object elements are the numerically expected results

    Code
      unclass(pclust)
    Output
      $scores
                k=2       k=3       k=4      k=5      k=6
      [1,] 4.600000 28.800000 17.524528 13.13455 11.90964
      [2,] 4.600000 16.280000  9.823729 15.13651 14.52632
      [3,] 4.877419  8.203636 17.293458 12.87273 12.96774
      [4,] 3.050000 34.480000 13.600000 12.76460 14.17241
      [5,] 2.462069 43.200000 13.554545 12.82949 10.25397
      
      $mean_scores
            k=2       k=3       k=4       k=5       k=6 
       3.917898 26.192727 14.359252 13.347575 12.766016 
      
      $ci95_scores
                 k=2       k=3      k=4      k=5      k=6
      2.5%  2.520862  9.011273 10.19681 12.77109 10.41954
      97.5% 4.849677 42.328000 17.50142 14.93631 14.49093
      
      $random_scores
                 k=2       k=3       k=4       k=5      k=6
      [1,]  3.013043  4.360000 10.800000  9.028571 22.43878
      [2,]  3.306667  5.369231  8.903226  9.887640 25.20670
      [3,]  2.568000  6.266667  8.364179  6.708934 16.32558
      [4,]  7.200000 10.074725 12.562238 16.080000 26.93023
      [5,] 11.742857  3.896774 12.891176 15.760428 11.96481
      
      $mean_random_scores
            k=2       k=3       k=4       k=5       k=6 
       5.566113  5.993479 10.704164 11.493115 20.573221 
      
      $D_max
            k=2       k=3       k=4       k=5       k=6 
      -1.648216 20.199248  3.655088  1.854460 -7.807205 
      
      $D_gap
              k=2         k=3         k=4         k=5         k=6 
      -34.1083048  34.1083048 -10.8217982  -0.4301176   0.4301176 
      
      $clust_iter
      [1] 2 3 4 5 6
      
      $repeats
      [1] 5
      
      $n_iter
      [1] 10
      
      $size
      [1] 6
      
      $call
      progeny_cluster(data = progeny_data, clust_iter = 2:6L, repeats = 5L, 
          r_seed = 101, n_iter = 10L, size = 6)
      

# `pclust` S3 print method returns expected known output

    == Progeny Clustering ==========================================================
    * Call              'progeny_cluster(data = progeny_data, clust_iter = 2:6L, repeats = 5L, r_seed = 101, n_iter = 10L, size = 6)'
    * Progeny Size      6
    * K iterations      '2 > 6'
    * No. iterations    10
    * No. repeats       5
    
    -- Mean & CI95 Stability Scores ------------------------------------------------
           k=2 *k=3*  k=4  k=5  k=6
    2.5%  2.52  9.01 10.2 12.8 10.4
    mean  3.92 26.19 14.4 13.3 12.8
    97.5% 4.85 42.33 17.5 14.9 14.5
    
    -- Maximum Distance Scores -----------------------------------------------------
      k=2 *k=3*   k=4   k=5   k=6 
    -1.65 20.20  3.66  1.85 -7.81 
    
    -- Gap Distance Scores ---------------------------------------------------------
       k=2  *k=3*    k=4    k=5    k=6 
    -34.11  34.11 -10.82  -0.43   0.43 
    ================================================================================

