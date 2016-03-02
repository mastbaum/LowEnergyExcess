#ifndef LARLITE_COSMICTRIGGERHACKER_CXX
#define LARLITE_COSMICTRIGGERHACKER_CXX

#include "CosmicTriggerHacker.h"

namespace larlite {

    bool CosmicTriggerHacker::initialize() {


        if (!_trigger_hack_tree) {
            _trigger_hack_tree = new TTree("trigger_hack_tree", "trigger_hack_tree");
            _trigger_hack_tree->Branch("hacked_trig_time", &hacked_trig_time, "hacked_trig_time/D");
        }

        return true;
    }

    bool CosmicTriggerHacker::analyze(storage_manager* storage) {

        hacked_trig_time = std::numeric_limits<double>::max();

        // Loop over mctracks to hack a trigger since in-time cosmics are taking goddamn forever to generate
        auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
        if (!ev_mctrack) {
            print(larlite::msg::kERROR, __FUNCTION__, Form("Did not find specified data product, track!"));
            return false;
        }

        // Make a vector of arrival time of all mctracks that pass thru the TPC
        std::vector<double> hacked_trig_times;
        hacked_trig_times.clear();
        for (auto const &mct : *ev_mctrack) {
            //if track doesn't go through TPC fuck it
            if (!mct.size()) continue;
            hacked_trig_times.push_back(mct.front().T());
        }

        // Choose a random arrival time from the list
        int randomIndex = rand() % hacked_trig_times.size();
        hacked_trig_time = hacked_trig_times[randomIndex];
        _trigger_hack_tree->Fill();

        return true;
    }

    bool CosmicTriggerHacker::finalize() {
        if (_trigger_hack_tree && _fout) {
            _fout->cd();
            _trigger_hack_tree->Write();
        }

        return true;
    }

}
#endif
