#ifndef _JETCLU_ALGORITHM_HH_
#define _JETCLU_ALGORITHM_HH_

#include "PhysicsTower.hh"
#include "Cluster.hh"
#include <vector>

class JetCluAlgorithm
{
 private:
  double _seedThreshold;
  double _coneRadius;
  int    _adjacencyCut;
  int    _maxIterations;
  int    _iratch;
  double _overlapThreshold;

 public:
  JetCluAlgorithm():
    _seedThreshold(1),
    _coneRadius(0.7),
    _adjacencyCut(2),
    _maxIterations(100),
    _iratch(1),
    _overlapThreshold(0.75)
  {}
  JetCluAlgorithm(double st, double cr, int ac, int mi, int ir, double ot):
    _seedThreshold(st),
    _coneRadius(cr),
    _adjacencyCut(ac),
    _maxIterations(mi),
    _iratch(ir),
    _overlapThreshold(ot)
  {}
  void makeSeedTowers(std::vector<PhysicsTower>& towers, std::vector<Cluster>& seedTowers);
  void buildPreClusters(std::vector<Cluster>& seedTowers, std::vector<PhysicsTower>& towers, std::vector<Cluster>& preClusters);
  void findStableCones(std::vector<Cluster>& preClusters, std::vector<PhysicsTower>& towers, std::vector<Cluster>& stableCones);
  void splitAndMerge(std::vector<Cluster>& stableCones, std::vector<Cluster>& jets);
  void run(std::vector<PhysicsTower>& towers, std::vector<Cluster>& jets);
};

#endif
