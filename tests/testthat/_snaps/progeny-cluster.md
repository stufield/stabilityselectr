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
                k=2      k=3       k=4       k=5      k=6
      [1,] 3.800000 14.47273 15.096774 13.952941 12.65831
      [2,] 3.156522  9.60000 10.279518 16.745455 15.37500
      [3,] 4.354286 15.69836 18.819802 15.365079 12.01818
      [4,] 3.234783 20.50000  8.557068  9.724907 14.90253
      [5,] 3.257143 34.32000 10.732919 13.430986 17.92500
      
      $mean_scores
            k=2       k=3       k=4       k=5       k=6 
       3.560547 18.918218 12.697216 13.843874 14.575803 
      
      $ci95_scores
                 k=2      k=3       k=4      k=5      k=6
      2.5%  3.164348 10.08727  8.729313 10.09551 12.08219
      97.5% 4.298857 32.93800 18.447499 16.60742 17.67000
      
      $random_scores
                k=2      k=3       k=4      k=5      k=6
      [1,] 2.987500 4.600000 10.264286 13.94286 17.25000
      [2,] 3.320930 5.213245  8.036364 10.86512 16.12928
      [3,] 2.828571 4.418182 13.552941 12.49258 18.19917
      [4,] 3.120000 8.022222  6.346988 18.34074 26.45455
      [5,] 2.469231 5.105732 11.023602 18.83774 23.72021
      
      $mean_random_scores
            k=2       k=3       k=4       k=5       k=6 
       2.945246  5.471876  9.844836 14.895805 20.350640 
      
      $D_max
             k=2        k=3        k=4        k=5        k=6 
       0.6153001 13.4463413  2.8523801 -1.0519317 -5.7748369 
      
      $D_gap
              k=2         k=3         k=4         k=5         k=6 
      -21.5786723  21.5786723  -7.3676586   0.4147277  -0.4147277 
      
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
           k=2 *k=3*   k=4  k=5  k=6
    2.5%  3.16  10.1  8.73 10.1 12.1
    mean  3.56  18.9 12.70 13.8 14.6
    97.5% 4.30  32.9 18.45 16.6 17.7
    
    -- Maximum Distance Scores -----------------------------------------------------
       k=2  *k=3*    k=4    k=5    k=6 
     0.615 13.446  2.852 -1.052 -5.775 
    
    -- Gap Distance Scores ---------------------------------------------------------
        k=2   *k=3*     k=4     k=5     k=6 
    -21.579  21.579  -7.368   0.415  -0.415 
    ================================================================================

