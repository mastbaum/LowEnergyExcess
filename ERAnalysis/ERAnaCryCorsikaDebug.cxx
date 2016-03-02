#ifndef ERTOOL_ERANACRYCORSIKADEBUG_CXX
#define ERTOOL_ERANACRYCORSIKADEBUG_CXX

#include "ERAnaCryCorsikaDebug.h"

namespace ertool {

  ERAnaCryCorsikaDebug::ERAnaCryCorsikaDebug(const std::string& name) : AnaBase(name)
  {}

  void ERAnaCryCorsikaDebug::Reset()
  {}

  void ERAnaCryCorsikaDebug::AcceptPSet(const ::fcllite::PSet& cfg)
  {}

  void ERAnaCryCorsikaDebug::ProcessBegin()
  {
    _vactive  = ::geoalgo::AABox(0,
                                 -larutil::Geometry::GetME()->DetHalfHeight(),
                                 0,
                                 2 * larutil::Geometry::GetME()->DetHalfWidth(),
                                 larutil::Geometry::GetME()->DetHalfHeight(),
                                 larutil::Geometry::GetME()->DetLength());

    n_with_parent_track = 0;
    n_with_parent_shower = 0;
    n_with_parent_uncontained_track = 0;
    n_with_parent_contained_track = 0;
    n_showers_with_start_not_contained = 0;
    n_showers_with_no_vtx_activity = 0;
    n_showers_from_pizeros = 0;
    n_evts = 0;

    _vtx_sphere.Radius(5.);
  }

  bool ERAnaCryCorsikaDebug::Analyze(const EventData &data, const ParticleGraph &ps)
  {
    n_evts++;
    // Get MC particle set
    auto const& mc_graph = MCParticleGraph();
    // Get the MC data
    auto const& mc_data = MCEventData();

    // Count the # of showers in the detector that don't have parents in the detector
    for ( auto const & mc : mc_graph.GetParticleArray() ) {

      if (mc.RecoType() == kShower) {
        ertool::Shower mc_ertoolshower = mc_data.Shower(mc.RecoID());

        if (IsVertexActivity(mc, mc_graph, mc_data))
          n_showers_with_no_vtx_activity ++;

        if (!_vactive.Contain(mc_ertoolshower.Start()))
          n_showers_with_start_not_contained++;

        // std::cout << mc_ertoolshower._energy << std::endl;
        auto parent = mc_graph.GetParticle(mc.Parent());
        if (parent.PdgCode() == 111)
          n_showers_from_pizeros++;

        if (parent.RecoType() == kShower) {
          // std::cout << "parent is shower" << std::endl;
          // std::cout << mc_data.Shower(parent.RecoID())._energy << std::endl;
          n_with_parent_shower++;
        }

        if (parent.RecoType() == kTrack) {
          n_with_parent_track++;
          // std::cout << "parent is track" << std::endl;
          // std::cout << mc_data.Track(parent.RecoID())._energy << std::endl;
          bool parent_contained = false;
          for (auto const & pt : mc_data.Track(parent.RecoID())) {
            if (_vactive.Contain(pt)) {
              parent_contained = true;
              break;
            }
          }
          if (!parent_contained) {
            std::cout << "Parent of shower is a track not at all contained in tpc!" << std::endl;
            std::cout << "Note, is the shower start contained? " << _vactive.Contain(mc_ertoolshower.Start()) << std::endl;
            n_with_parent_uncontained_track++;
          }
          else
            n_with_parent_contained_track++;

        }
      }
    }

    return true;
  }

  void ERAnaCryCorsikaDebug::ProcessEnd(TFile* fout)
  {
    std::cout << "END!" << std::endl;
    std::cout << "per event, n_showers_with_start_not_contained = "
              << n_showers_with_start_not_contained / float(n_evts) << std::endl;
    std::cout << "per event, n_with_parent_contained_track = "
              << n_with_parent_contained_track / float(n_evts) << std::endl;
    std::cout << "per event, n_with_parent_uncontained_track = "
              << n_with_parent_uncontained_track / float(n_evts) << std::endl;
    std::cout << "per event, n_with_parent_shower = "
              << n_with_parent_shower / float(n_evts) << std::endl;
    std::cout << "per event, n_with_parent_track = "
              << n_with_parent_track / float(n_evts) << std::endl;
    std::cout << "per event, n_showers_with_no_vtx_activity = "
              << n_showers_with_no_vtx_activity / float(n_evts) << std::endl;
    std::cout << "per event, n_showers_from_pizeros = " 
              << n_showers_from_pizeros / float(n_evts) << std::endl;
  }

  bool ERAnaCryCorsikaDebug::IsVertexActivity(const Particle &shower, const ParticleGraph &ps, const EventData &data) {
    // Center the sphere on the shower start point
    _vtx_sphere.Center(shower.Vertex());

    // Loop over all tracks and showers in the event (excluding the shower's electron)
    // if any part of a track passes thru sphere, add that track's energy to vertex energy
    // if start point of a shower is in sphere, add that track's energy to vertex energy
    for (auto const& track : data.Track()) {
      // Only consider tracks that are longer than 0.3cm!
      if (track.Length() < 0.3) continue;

      for (auto const& pt : track)
        if (_vtx_sphere.Contain(pt))
          return true;
    }

    for (auto const& myshower : data.Shower()) {
      /// Don't include the "SingleE" energy in vertex energy calculation
      if (myshower.RecoID() == shower.RecoID())
        continue;
      if (_vtx_sphere.Contain(myshower.Start()))
        return true;
    }

    return false;
  }
}

#endif
