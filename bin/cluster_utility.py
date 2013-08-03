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
    A utility for clustering gaps.
    This class wraps DBSCAN implementation from scikit-learn.
     
    default min samples is set at number of dimensions + 1 
            (i.e. 3 = 2 + 1)
    default epsilon is based on fact that we are not scaling the matrix
        of (gap_start, gap_width) data and also the idea that true split 
        reads should be piling up literally on top of one another.
    Reasonable clustering behavior was confirmed through visual 
        inspection of several representative transcripts."""
    
    class ShuntLogger():
        def log(self, message):
            pass
    
    def __init__(self, epsilon=3, min_samples=3, 
            deduplication_threshold=10000, logger=ShuntLogger()):
        #Explicitily seeding ensures deterministic clusters; without seeding
        #   boundary points may flop between clusters and cluster labels
        #   invariably shift from run to run (innocuous, but very annoying 
        #   when comparing results and testing). 
        random_seed = 42
        random_state = np.random.RandomState(seed = random_seed)
        self._min_samples = min_samples
        self._dbscan = DBSCAN(epsilon, min_samples, random_state=random_state)
        self._deduplication_threshold = deduplication_threshold
        self._logger = logger 
        self._logger.log( 
            "DBSCAN clustering with epsilon [{0}] and min_samples [{1}] " 
                "using random seed [{2}]" 
                .format(epsilon, min_samples, random_seed))

    @staticmethod
    def _gaps_in_a_chromosome(gaps_sorted_by_chromosome):
        gaps=[]
        # pylint: disable=line-too-long
        current_chromosome = "" if len(gaps_sorted_by_chromosome) == 0 else gaps_sorted_by_chromosome[0].chromosome
        for gap in gaps_sorted_by_chromosome:
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

    def _deduplicate_and_cluster_gaps(self, gaps):
        """Conditionally deduplicates gaps before processing. 
        This optimization is reasonable because for split read gap data, MANY
        gaps appear at the same coordinate (gap start, gap_width) so there is
        huge "duplication" of the clustered samples, and duplication above the
        minimum number of samples for a cluster doesn't materially impact the
        clustering results. So instead of processing all the gaps, we reduce it
        to a representative subset, cluster those and then assign clusters to
        their excluded brethren."""

        if len(gaps) < self._deduplication_threshold:
            return self._cluster_gaps(gaps)
    
        deduped_gap_dict={}
        for gap in gaps:
            gap_coordinate = (gap.gap_start, gap.gap_width())
            gaps_at_this_coordinate = deduped_gap_dict.setdefault(gap_coordinate,[])
            if len(gaps_at_this_coordinate) < self._min_samples:
                deduped_gap_dict[gap_coordinate].append(gap)

        deduped_gaps = []
        for coordinate, gaps_at_this_coordinate in deduped_gap_dict.iteritems():
            deduped_gaps.extend(gaps_at_this_coordinate)
            
        self._logger.log(
            "Deduplication reduced gap count from {0} to {1}". \
            format(len(gaps), len(deduped_gaps)))

        cluster_count = self._cluster_gaps(deduped_gaps)
            
        for gap in gaps:
            if gap.cluster == -1:
                gap_coordinate = (gap.gap_start, gap.gap_width())
                gap.cluster = deduped_gap_dict[gap_coordinate][0].cluster

        return cluster_count

        
    def assign_clusters(self, gaps_sorted_by_chromosome):
        chromosome_count = 0
        total_chromosome_count = \
            len(set([gap.chromosome for gap in gaps_sorted_by_chromosome]))
        for (chromosome,gaps) in DbscanClusterUtility. \
                _gaps_in_a_chromosome(gaps_sorted_by_chromosome):
            chromosome_count += 1
            self._logger.log(
                "clustering {0} gaps in chromosome {1} ({2}/{3})". \
                format(len(gaps), chromosome, chromosome_count, 
                    total_chromosome_count))
            start_time = time.clock()
            cluster_count = self._deduplicate_and_cluster_gaps(gaps)
            self._logger.log(
                "found {0} clusters for {1} gaps in "
                "chromosome {2} in {3:.0f} seconds". \
                format(cluster_count, len(gaps), chromosome, 
                    time.clock()-start_time))
