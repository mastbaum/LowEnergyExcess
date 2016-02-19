#ifndef LARLITE_ERSELSAVESINGLEEEVENTS_CXX
#define LARLITE_ERSELSAVESINGLEEEVENTS_CXX

#include "ERSelSaveSingleEEvents.h"
#include "DataFormat/mcshower.h"

namespace larlite {
  
  ERSelSaveSingleEEvents::ERSelSaveSingleEEvents( const ::ertool::io::StreamType_t in_strm,
					  const ::ertool::io::StreamType_t out_strm)
    : ERToolAnaBase(in_strm,out_strm)
  { 
    _name="ERSelSaveSingleEEvents"; 
  }

  bool ERSelSaveSingleEEvents::initialize() {

    return ERToolAnaBase::initialize();

  }
  
  bool ERSelSaveSingleEEvents::analyze(storage_manager* storage) {

    auto status = ERToolAnaBase::analyze(storage);
    if(!status) return status;

    bool mgr_status = _mgr.Process();
    if(!mgr_status) return mgr_status;

    // Manager is returning true, which happens when ERAna returns true
    // ERAnaLowEnergyExcess returns true when a nue is reconstructed
    // Get the particle graph, find the electron associated with the neutrino,
    // and store the electron's energy to double precision.
    // Then look for the mcshower in the event with that same energy
    // and set that mcshower's PDG to "1234567" to be saved in output
    auto const& graph = _mgr.ParticleGraph();
    auto const& particles = graph.GetParticleArray();
    auto const& data = _mgr.EventData();
    double mcshower_energy = -1.;

    // Loop over particles and find the nue's electron, save its energy
    for ( auto const & p : particles ) {
      if ( abs(p.PdgCode()) == 12 ) {
        // Loop over the neutrinos immediate children and find the electron
        for (auto const& d : p.Children()) {
          auto daught = graph.GetParticle(d);
          // This is the "ccsinglee" electron.
          if (daught.PdgCode() == 11) {
            mcshower_energy = data.Shower(daught.RecoID())._energy;
            break;
          }
        }
      }
    }

    // Loop over mcshowers and find the mcshower with the same energy, and alter its pdg
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    if (!ev_mcs) {
        print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, MCShower!"));
        return false;
    }
    if (!ev_mcs->size())
        return false;

    for (auto &mcs : *ev_mcs){
      if (mcs.DetProfile().E() == mcshower_energy){
        // std::cout<<"Found a mcshower with matching energy to ertool shower associated with nue!"<<std::endl;
        mcs.PdgCode(1234567);
      }
    }
    return mgr_status;
  }

  bool ERSelSaveSingleEEvents::finalize() {

    return ERToolAnaBase::finalize();

  }

}
#endif
