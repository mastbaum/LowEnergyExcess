#ifndef ERTOOL_ERALGOTAGEMULATEDDELETIONSCOSMIC_CXX
#define ERTOOL_ERALGOTAGEMULATEDDELETIONSCOSMIC_CXX

#include "ERAlgoTagEmulatedDeletionsCosmic.h"

namespace ertool {

  ERAlgoTagEmulatedDeletionsCosmic::ERAlgoTagEmulatedDeletionsCosmic(const std::string& name) : AlgoBase(name)
  {}

  void ERAlgoTagEmulatedDeletionsCosmic::Reset()
  {}

  void ERAlgoTagEmulatedDeletionsCosmic::AcceptPSet(const ::fcllite::PSet& cfg)
  {}

  void ERAlgoTagEmulatedDeletionsCosmic::ProcessBegin()
  {}

  bool ERAlgoTagEmulatedDeletionsCosmic::Reconstruct(const EventData &data, ParticleGraph& graph)
  {

    for (auto const& id : graph.GetParticleNodes()) {

      auto& part = graph.GetParticle(id);

      if (part.RecoType() == kTrack)
        if ( data.Track(part.RecoID()).front().at(1) < -999. )
          part.SetProcess(kCosmic);
      if (part.RecoType() == kShower)
        if ( data.Shower(part.RecoID()).Start().at(1) < -999. )
          part.SetProcess(kCosmic);

    }

    return true;

  }

  void ERAlgoTagEmulatedDeletionsCosmic::ProcessEnd(TFile * fout)
  {}

}

#endif
