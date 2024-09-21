# `pclust` S3 print method returns expected known output

    == Progeny Cluster Object ======================================================
       Call                     progenyCluster(data = progeny_data, clust.iter = 2:6, reps = 5, iter = 10, size = 6)
       Progeny Size             6
       No. of Iterations        10
       K Iterations             2 3 4 5 6
    
    == Mean & CI95 Stability Scores ================================================
           k=2 k=3*   k=4  k=5  k=6
    2.5%  2.25 15.2  6.17 13.5 11.9
          3.22 28.1 15.87 18.5 13.1
    97.5% 4.23 54.7 25.14 27.4 14.6
    
    == Maximum Distance Scores =====================================================
      k=2  k=3*   k=4   k=5   k=6 
    -3.38 20.90  3.97  4.00  1.72 
    
    == Gap Distance Scores =========================================================
       k=2   k=3*    k=4    k=5    k=6 
    -37.08  37.08 -14.86   8.04  -8.04 
    ================================================================================

