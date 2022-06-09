# GTWIG

Cluster whistle traces by the Graph-based Traces to Whistles Gathering algorithm (G-TWIG)

This is an implementation of the algorithm in:
D. Kipnis and R. Diamant, "Graph-Based Clustering of Dolphin Whistles." IEEE/ACM Transactions on Audio, Speech, and Language Processing 29 (2021): 2216-2227.

To fully run this algorithm you need an external inplementation of fractional Fourier transform, which is not included in this repository.
I used the function fracF.m that can be downloaded from:
http://www.ee.bilkent.edu.tr/~haldun/wileybook.html
(Matlab code for fast computation of the fractional Fourier transform)
Alternately, you can use other implementations, or disable the trace-continuity calculation by setting: AlgParams.calcContinuityFlag=0 (which is the default)

Expected results when running on demo:
For AlgParams.calcContinuityFlag=0, there should be clusters at 56s (multipath), 66.1s and 77.6s (harmonics).
For AlgParams.calcContinuityFlag=1, there should be clusters at 56s (multipath), 77.6s (harmonics), 84s (connectivity and harmonics) and 252.7s (connectivity) 
