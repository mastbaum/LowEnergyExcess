/**
 * \file LEERW.h
 *
 * \ingroup LEEReweight
 *
 * \brief This class gets handed a neutrino interaction (mctruth object) and computes a weight.
 *  The weight takes into account the energy and angle of the outgoing electron in the interaction
 *  as well as cross section scalings and flux scalings from MiniBooNE to MicroBooNE.
 *  This weight will sculpt the neutrino energy spectrum to a shape that matches the MB scaled excess prediction.
 *  Knowing how many excess events are expected in MicroBooNE (from the scaling prediction), the weight is
 *  multiplied such that the final integral of the energy spectrum is what is expected.
 *  The idea is the user can generate nue interactions, apply some sort of filter to select one that
 *  pass MiniBooNE event selection cuts (one subevent, etc)
 *  then hand the neutrino interaction to this package to get a weight. Looping over events and
 *  histogramming them with this weight will produce a simulated low energy excess in MicroBooNE.
 *  (of course assuming the excess comes from an excess of nue interactions)
 *
 *  This takes a while to initialize. Its initialization loops over ALL of the events and stores the weight
 *  for each event. This is necessary because for absolute normalization, one needs to know the sum of all weights
 *  a priori.
 *
 * @author davidkaleko
 */

/** \addtogroup LEEReweight

    @{*/
#ifndef LEERW_H
#define LEERW_H

#include <iostream>
#include "Utilities/PlotReader.h"
#include "Base/larlite_base.h"
#include "DataFormat/mctruth.h"
#include "TH2.h"
#include "TGraph.h"

/**
   \class LEERW
   User defined class LEERW ... these comments are used to generate
   doxygen documentation!
 */
namespace lee {


	struct EventInfo_t {
		double electron_energy_MEV;
		double electron_uz;
		double nue_energy_GEV;
		bool   is_there_a_neutrino;
		size_t n_electrons = 0;
	};

	class LEERW : public larlite::larlite_base {

	public:

		/// Default constructor
		LEERW() {
			_source_filename = "";
			_events_generated_only_in_TPC = true;
		}

		/// Default destructor
		~LEERW() {}

		bool initialize();

		//This weight sculpts the energy/angle spectrum of the neutrino interaction to match that of miniboone
		//This weight (for now) is NOT NORMALIZED CORRECTLY. The user needs to do a separate event loop and sum up
		//all of the sculpting weights of all events, then scale the resulting histogram by n_events_analyzed/summed_weight
		double get_sculpting_weight(const larlite::mctruth* mytruth);

		//This weight is correctly normalized, and has to do with number of excess events (and some more sculpting) expected
		//to see in MicroBooNE after 6.6 POT running. The user can just use this weight for each event and everything will be fine.
		double get_normalized_weight(const larlite::mctruth* mytruth);

		double get_sculpting_weight(double electron_energy_MEV, double electron_uz);
		double get_normalized_weight(double neutrino_energy_GEV);

		void set_source_filename(std::string filename) { _source_filename = filename; }

		void set_n_generated_events(size_t david) { _n_generated_evts = david; }
		void set_generated_evis_uz_corr_name(std::string name) { _generated_evis_uz_corr_name = name; }
		
		/// If the generated neutrino sample was only generated in the TPC active volume,
		/// or if it was generated throughout the entire cryostat volume.
		void set_events_generated_only_in_TPC(bool kaleko) { _events_generated_only_in_TPC = kaleko; }

		void set_debug(bool david) { _debug = david; }

		/// Utility function to grab relevant stuff for reweighting from the mctruth object
		/// public for hacky reasons right now
		const EventInfo_t extract_event_info(const larlite::mctruth* mytruth);

		//microboone tonnage (total) in grams
		//from 1.3954 [g/cm^3] * pi * 190^2 * 1079 [cm^3]
		double UB_TONNAGE_GRAMS = 170756353.192;

		//miniboone tonnage (total) in grams
		//from 0.855 [g/cm^3] * (4/3) * pi * 610.6^3 [cm^3]
		double MINIBOONE_TONNAGE_GRAMS = 815313732.1;


	private:

		//Check if the LEERW package is properly initialized
		void check_is_initialized();

		//Utility to print event info to screen
		void print_evt_info(const EventInfo_t evt_info);

		TGraph _flux_ratio;
		TGraph _xsec_ratio;
		TH2D   _MB_evis_uz_corr;
		TH2D   _generated_evis_uz_corr;

		/// Number of generated nue events (for absolute normalization)
		double _n_generated_evts = 0;

		/// Name of input root file containing scaling flux ratio graphs, xsec ratio graphs, evis_uz_correlation histo
		std::string _source_filename;
		std::string _flux_ratio_name = "flux_ratio";
		std::string _xsec_ratio_name = "xsec_ratio";
		std::string _MB_evis_uz_corr_name = "hist_raw_uz_evis_smooth";
		std::string _generated_evis_uz_corr_name = "";//"initial_evis_uz_corr";

		double _pot_weight = 1.021671826625387; //6.6 / 6.46; //microboone POT over miniboone POT
		double _tonnage_weight = 0.20943637580123128;//UB_TONNAGE_GRAMS / MINIBOONE_TONNAGE_GRAMS; //roughly 0.2
		double _true_MB_excess_evts = 1212.114;
		bool _events_generated_only_in_TPC;

		bool _debug = false;

	};

} // end namespace lee
#endif
/** @} */ // end of doxygen group

