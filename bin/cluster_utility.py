"""
cluster_utility.py
7/28/2013- cgates

Wraps the DBSCAN clustering implementation in a more abstract interface to 
enable modular clustering substitution and mocking.
"""
import time
import numpy as np
from scipy.spatial import distance
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler


class DbscanClusterUtility():
    """
    default min samples is set at number of dimensions + 1 
            (i.e. 3 = 2 + 1)
    default epsilon is based on fact that we are not scaling the matrix
        of (gap_start, gap_width) data and also the idea that true split 
        reads should be piling up literally on top of one another.
    Reasonable clustering behavior was confirmed through visual 
        inspection of several representative transcripts."""
    
    class ShuntLogger():
        def log(self):
            pass
    
    def __init__(self, epsilon=3, min_samples=3, logger=None):
        #Explicitily seeding ensures deterministic clusters; without seeding
        #   boundary points may flop between clusters and cluster labels
        #   invariably shift from run to run (innocuous, but very annoying 
        #   when comparing results or otherwise testing). 
        random_seed = 42
        random_state = np.random.RandomState(seed = random_seed)
        self._dbscan = DBSCAN(epsilon, min_samples, random_state = random_state)
        self._logger = ShuntLogger() if logger == None else logger 
        self._logger.log(
            "DBSCAN clustering with epsilon [{0}] and min_samples [{1}] "
                "using random seed [{2}]"
                .format(epsilon, min_samples, random_seed))

    @staticmethod
    def _chromosome_partition(sorted_gaps):
        """assumes gaps are sorted by chromosome"""
        gaps=[]
        # pylint: disable=line-too-long
        current_chromosome = "" if len(sorted_gaps) == 0 else sorted_gaps[0].chromosome
        for gap in sorted_gaps:
          if gap.chromosome != current_chromosome:
            yield (current_chromosome, gaps)
            current_chromosome = gap.chromosome
            gaps = [gap]
          else:
            gaps.append(gap)
        yield (current_chromosome, gaps)

    def _cluster_gaps(self, chromosome_gaps):
        """Identifies clusters in specified gaps, and assigns cluster
        to each gap, returns count of clusters identified in this set. 
        Valid clusters are ints starting with 0. The cluster value -1 
        represents a "noise gap", i.e. a gap which was not placed into a 
        cluster.
        Note that because this method does not use an explicit distance 
        matrix, it will tolerate edge cases where chromsomes have only one gap
        or only several "identical" gaps. (Gaps in these chromosomes would be
        assigned clusters of -1.)
        """
        matrix=[]
        for gap in chromosome_gaps:
            matrix.append(
                [float(gap.gap_start), float(gap.gap_width())])
        clusters = self._dbscan.fit_predict(np.matrix(matrix))
        
        for gap, cluster in zip(chromosome_gaps, clusters):
            gap.cluster = int(cluster)
            
        # Number of clusters in labels, ignoring noise if present.
        cluster_count = len(set(clusters)) - (1 if -1 in clusters else 0)        
        return cluster_count
        
    def assign_clusters(self, sorted_gaps):
        """assumes gaps are sorted by chromosome"""
        chromosome_count = 0
        # pylint: disable=line-too-long
        total_chromosome_count = len(set([gap.chromosome for gap in sorted_gaps]))
        # pylint: disable=line-too-long
        for (chromosome,gaps) in DbscanClusterUtility._chromosome_partition(sorted_gaps):
            chromosome_count += 1
            if len(gaps) > 100000:
                self._logger.log(
                    "WARNING skipping clustering for {0} gaps in "
                    "chromosome {1} ({2}/{3})". \
                    format(len(gaps), chromosome, chromosome_count, 
                        total_chromosome_count))
                continue
            self._logger.log(
                "clustering {0} gaps in chromosome {1} ({2}/{3})". \
                format(len(gaps), chromosome, chromosome_count, 
                    total_chromosome_count))
            start_time = time.clock()
            cluster_count = self._cluster_gaps(gaps)
            self._logger.log(
                "found {0} clusters for {1} gaps in "
                "chromosome {2} in {3:.0f} seconds". \
                format(cluster_count, len(gaps), chromosome, 
                    time.clock()-start_time))
