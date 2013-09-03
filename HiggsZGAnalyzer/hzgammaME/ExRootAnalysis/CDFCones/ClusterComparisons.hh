#ifndef _CLUSTER_COMPARISONS_HH_
#define _CLUSTER_COMPARISONS_HH_

#include "Cluster.hh"

class ClusterFourVectorEtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.fourVector.Et() > c2.fourVector.Et();
  }
};

class ClusterCentroidEtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.centroid.Et > c2.centroid.Et;
  }
};

class ClusterPtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.fourVector.pt() > c2.fourVector.pt();
  }
};

#endif
