#ifndef LARLITE_MC_COSMIC_FILTER_CXX
#define LARLITE_MC_COSMIC_FILTER_CXX

#include "MC_cosmic_Filter.h"
#include "DataFormat/mctrack.h"

namespace larlite {

bool MC_cosmic_Filter::initialize() {

  _n_total_events = 0;
  _n_kept_events = 0;

  if(!_trigger_hack_tree){
    _trigger_hack_tree = new TTree("trigger_hack_tree","trigger_hack_tree");
    _trigger_hack_tree->Branch("hacked_trig_time",&hacked_trig_time,"hacked_trig_time/D");
  }

  return true;
}

bool MC_cosmic_Filter::analyze(storage_manager* storage) {

  hacked_trig_time = std::numeric_limits<double>::max();

  auto ev_mctruth = storage->get_data<event_mctruth>("generator");
  if (!ev_mctruth) {
    print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, mctruth!"));
    return false;
  }

  _n_total_events++;

  bool ret = false;

  for (auto const &mctruth : *ev_mctruth) {
  
    //Enforce Cosmic Origins
    if (mctruth.Origin() == simb::Origin_t::kCosmicRay) {ret = true;}
    //Corsika cosmics have origin == 0
    else if (mctruth.Origin() == simb::Origin_t::kUnknown) {
      print(larlite::msg::kWARNING, __FUNCTION__, Form("MCTruth origin unknown! Allowing it to pass cosmic filter."));
      ret = true;
    }
    else {ret = false;}

    if (ret == true) continue;

  }
  

  if (ret) _n_kept_events++;
  if(!ret) return false;

  // Loop over mctracks to hack a trigger since in-time cosmics are taking goddamn forever to generate
   auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
  if (!ev_mctrack) {
    print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
    return false;
  }

  // Make a vector of arrival time of all mctracks that pass thru the TPC
  std::vector<double> hacked_trig_times;
  hacked_trig_times.clear();
  for (auto const &mct : *ev_mctrack){
    //if track doesn't go through TPC fuck it
    if(!mct.size()) continue;
    hacked_trig_times.push_back(mct.front().T());
  }

  // Choose a random arrival time from the list
  int randomIndex = rand() % hacked_trig_times.size();
  hacked_trig_time = hacked_trig_times[randomIndex];
  _trigger_hack_tree->Fill();

  return true;

}

bool MC_cosmic_Filter::finalize() {

  std::cout << _n_total_events << " total events analyzed, " << _n_kept_events << " events passed MC_cosmic_Filter." << std::endl;

  if(_trigger_hack_tree && _fout){
    _fout->cd();
    _trigger_hack_tree->Write();
  }
  return true;
}

}
#endif
