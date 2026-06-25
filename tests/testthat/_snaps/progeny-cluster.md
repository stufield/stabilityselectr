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
      [1,] 3.293023 58.000000 13.146667 20.87671 13.16933
      [2,] 2.723077 20.032653  5.393594 15.78610 14.77419
      [3,] 4.300000 24.526829 15.858621 18.63354 11.32584
      [4,] 3.868421  9.909677  7.152212 10.56371 11.85498
      [5,] 2.862500 16.160000  9.844068 13.01918 13.83221
      
      $mean_scores
            k=2       k=3       k=4       k=5       k=6 
       3.409404 25.725832 10.279032 15.775847 12.991313 
      
      $ci95_scores
                 k=2      k=3       k=4      k=5      k=6
      2.5%  2.737019 10.53471  5.569456 10.80925 11.37876
      97.5% 4.256842 54.65268 15.587425 20.65240 14.68000
      
      $random_scores
                k=2      k=3      k=4       k=5      k=6
      [1,] 2.284211 7.260504 5.787170 12.357895 19.81081
      [2,] 1.837500 3.837563 8.205076 11.084337 16.62791
      [3,] 2.080645 6.194118 6.801770  9.722182 12.61146
      [4,] 2.166102 3.210526 7.454717 12.164103 23.16495
      [5,] 4.235294 3.857592 9.515676 19.509677 13.15385
      
      $mean_random_scores
            k=2       k=3       k=4       k=5       k=6 
       2.520750  4.872061  7.552882 12.967639 17.073795 
      
      $D_max
             k=2        k=3        k=4        k=5        k=6 
       0.8886539 20.8537713  2.7261507  2.8082079 -4.0824825 
      
      $D_gap
             k=2        k=3        k=4        k=5        k=6 
      -37.763227  37.763227 -20.943614   8.281348  -8.281348 
      
      $clust_iter
      [1] 2 3 4 5 6
      
      $repeats
      [1] 5
      
      $n_iter
      [1] 10
      
      $size
      [1] 6
      
      $call
      progeny_cluster(data = progeny_data, clust_iter = 2:6L, reps = 5L, 
          n_iter = 10L, size = 6)
      

# `pclust` S3 print method returns expected known output

    == Progeny Clustering ==========================================================
    * Call              'progeny_cluster(data = progeny_data, clust_iter = 2:6L, reps = 5L, n_iter = 10L, size = 6)'
    * Progeny Size      6
    * K iterations      '2 > 6'
    * No. iterations    10
    * No. repeats       5
    
    -- Mean & CI95 Stability Scores ------------------------------------------------
           k=2 *k=3*   k=4  k=5  k=6
    2.5%  2.74  10.5  5.57 10.8 11.4
    mean  3.41  25.7 10.28 15.8 13.0
    97.5% 4.26  54.7 15.59 20.7 14.7
    
    -- Maximum Distance Scores -----------------------------------------------------
       k=2  *k=3*    k=4    k=5    k=6 
     0.889 20.854  2.726  2.808 -4.082 
    
    -- Gap Distance Scores ---------------------------------------------------------
       k=2  *k=3*    k=4    k=5    k=6 
    -37.76  37.76 -20.94   8.28  -8.28 
    ================================================================================

