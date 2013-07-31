import unittest
from bin.cluster_utility import DbscanClusterUtility

class ClusterUtilityTestCase(unittest.TestCase):

    def test_cluster_gaps_simpleCluster(self):
        #ENSMUST00000000305 for Samples 21786, 21797
        #using eps=3, minSample=3, random_seed=42
        length = 5
        gaps = [MockGap(start, length) for start in [
            192,192,192,
            193,193,193,
            194,194,194,
            195,195,195,
            196,196,196
            ]]

        dbscan = DbscanClusterUtility(
            epsilon=3, min_samples=3, logger=MockLogger())
        actual_cluster_count = dbscan._cluster_gaps(gaps)
        
        self.assertEquals(1, actual_cluster_count)
        #plot_clusters("ENSMUST00000000305", gaps, dbscan._dbscan)

    def test_cluster_gaps_singletonCluster(self):
        length = 5
        gaps = [MockGap(start, length) for start in [42,42,42]]

        dbscan = DbscanClusterUtility(
            epsilon=3, min_samples=3, logger=MockLogger())
        actual_cluster_count = dbscan._cluster_gaps(gaps)
        
        self.assertEquals(1, actual_cluster_count)

    def test_cluster_gaps_noClusters(self):
        length = 5
        gaps = [MockGap(start, length) for start in [42,42,46]]

        dbscan = DbscanClusterUtility(
            epsilon=3, min_samples=3, logger=MockLogger())
        actual_cluster_count = dbscan._cluster_gaps(gaps)
        
        self.assertEquals(0, actual_cluster_count)

    def test_cluster_standardCluster(self):
        #the distinct start,lengths from 
        #ENSMUST00000063084 for Samples 21786, 21797
        gaps = [MockGap(start, length) for (start,length) in [
            (811,26),(812,26),(813,26),(814,26),(815,26),(816,26),(817,26),
                (818,26),(819,26),(820,26),(821,26),
            (1554,2),(1555,2),(1556,2),(1557,2),(1558,2),
            (773,146),(774,146),(775,146),(776,146),
            (918,67),(919,67),(920,67),(921,67),(922,67),
            (1324,2),(1325,2),(1326,2),(1327,2),(1328,2),
            ]]
        dbscan = DbscanClusterUtility(
            epsilon=3, min_samples=3, logger=MockLogger())
        actual_cluster_count = dbscan._cluster_gaps(gaps)
        
        self.assertEquals(5, actual_cluster_count)
        #plot_clusters("ENSMUST00000063084", gaps, dbscan._dbscan)


    def test_cluster_gaps_complexCluster(self):
        #ENSMUST00000012259 for Samples 21786, 21797
        #using eps=3, minSample=3, random_seed=42
        #Individual reads were re-ordered to group with expected clusters
        gaps = [MockGap(start, length) for (start,length) in [
            (176,1216), (177,1216), #noise
            (1516,27), (1516,24), (1516,24), (1516,27), (1517,27), (1517,24), 
                (1517,24), (1517,27), (1518,27), (1518,24), (1518,24),
                (1518,27), (1519,27), (1519,24), (1519,24), (1519,27),
                (1520,27), (1520,24), (1520,24), (1520,27), (1521,27),
                (1521,24), (1521,24), (1521,27), (1522,27), (1522,24),
                (1522,27), (1523,27), (1524,27), #cluster 0
            (1142,249), (1142,249), (1142,249), (1142,249), (1142,249), 
                (1142,249), (1142,249), (1142,249), (1142,249), (1142,249),
                (1142,249), (1143,249), (1143,249), (1143,249), (1143,249),
                (1143,249), (1143,249), (1143,249), (1143,249), (1143,249),
                (1144,249), (1144,249), (1144,249), (1144,249), (1144,249),
                (1144,249), (1144,249), (1144,249), (1144,249), (1145,249),
                (1145,249), (1145,249), (1145,249), (1145,249), (1145,249),
                (1145,249), (1145,249), (1145,249), (1146,249), (1146,249),
                (1146,249), (1146,249), (1146,249), (1146,249), (1146,249),
                (1146,249), (1146,249), (1146,249), (1147,249), (1147,249),
                (1147,249), (1147,249), (1147,249), (1147,249), (1147,249),
                (1147,249), (1147,249), (1147,249), (1148,249), (1148,249),
                (1148,249), (1148,249), (1148,249), (1148,249), (1148,249),
                (1148,249), (1148,249), #cluster 1
            (592,81), (593,81), (594,81), (595,81), (596,81), (597,81), 
                (598,81), (599,81),(600,81), (601,81), (602,81), (603,81),
                (604,81), (605,81), (606,81), (607,81),(608,81), #cluster 2
            (314,30), (314,30), (315,30), (315,30), (316,30), (316,30),
                (317,30), (317,30), (318,30), (318,30), #cluster 3
            (1780,4), (1781,4), (1782,4), (1783,4), (1784,4), (1785,4), 
                (1786,4), #cluster 4 
            (1388,43), (1389,43), (1390,43), (1391,43), (1392,43), (1393,43), (1394,43),
                (1395,43), (1396,43), (1397,43), (1398,43), #cluster 5
            (1142,129), (1142,129), (1143,129), (1143,129),(1145,129), 
                (1144,129), #cluster 6
            (345,213), (345,213), (345,213), (345,213),
                (345,213), (345,213), (345,213), (345,213), (345,213),
                (345,213), (345,213), (345,213), (345,213), (345,213),
                (346,213), (346,213), (346,213), (346,213), (346,213),
                (346,213), (346,213), (346,213), (346,213), (346,213),
                (346,213), (346,213), (346,213), (346,213), (346,213),
                (346,213), (347,213), (347,213), (347,213), (347,213),
                (347,213), (347,213), (347,213), (347,213), (347,213),
                (347,213), (347,213), (347,213), (347,213), (347,213),
                (347,213), #cluster 7
            (175,88), (175,88), (176,88), (176,88), (177,88), #cluster 8
            (3053,2), (3054,2), (3055,2), (3056,2), #cluster 9
            (2072,2), (2073,2), (2074,2), (2075,2), (2076,2), #cluster 10
            (583,4), (584,4), (585,4), #cluster 11
            (558,713), (559,713), (560,713), (561,713), #cluster 12 
            (2240,99), (2241,99), (2242,99), (2243,99), (2244,99), (2245,99), 
                (2246,99), #cluster 13
            (344,926), (345,926), (346,926), (347,926), (348,926), #cluster 14
            (1912,64), (1913,64),(1914,64), (1915,64),  #cluster 15
            (1089,2), (1090,2), (1091,2), (1092,2), (1093,2), #cluster 16 
            (363,2), (364,2), (365,2), (366,2), (367,2), #cluster 17
            ]]

        dbscan = DbscanClusterUtility(
            epsilon=3, min_samples=3, logger=MockLogger())
        actual_cluster_count = dbscan._cluster_gaps(gaps)
        
        self.assertEquals(18, actual_cluster_count)
        #plot_clusters("ENSMUST00000012259", gaps, dbscan._dbscan)


def plot_clusters(title, gaps, dbscan):
    """Renders the clusters to GUI. 
    Note this will not work in a headless environ but could be adapted 
    to use matplotlib (instead of pylab) to enable creation of graphic 
    files instead of a GUI."""
    import numpy as np
    import pylab as pl
    from itertools import cycle


    core_samples = dbscan.core_sample_indices_
    labels = dbscan.labels_
    
    pl.close('all')
    pl.figure(1)
    pl.clf()

    colors = cycle("aqua,blue,fuchsia,gray,green,lime,maroon,navy,olive,orange,purple,red,silver,teal,white,yellow".split(","))
    markershapes = cycle('ov^<>sp*hHDd')
    for k, col, markershape in zip(set(labels), colors, markershapes):
        if k == -1: #this point is noise
            col = 'black'
            markersize = 3
            markershape = 'x' 
        class_members = [index[0] for index in np.argwhere(labels == k)]
        cluster_core_samples = [index for index in core_samples
                                if labels[index] == k]
        for index in class_members:
            gap = gaps[index]
            if index in core_samples and k != -1:
                markersize = 10
            else:
                markersize = 10
            pl.plot(gap.gap_start, gap.gap_width(), 
                    marker=markershape, markerfacecolor='white',
                    markeredgecolor=col, markersize=markersize)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    pl.title("{0} ({1} clusters)".format(title,n_clusters))
    pl.xlabel('gap_start (base position)')
    pl.ylabel('gap_length (base count)')
    pl.show()


class MockGap():
    
    def __init__(self, gap_start, gap_width): 
        self.gap_start = gap_start
        self._gap_width = gap_width
        self.cluster = -1

    def gap_width(self):
        return self._gap_width

class MockLogger():
    def log(self, message):
        pass
