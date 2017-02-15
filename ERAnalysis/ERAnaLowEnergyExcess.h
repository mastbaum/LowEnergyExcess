/**
 * \file ERAnaLowEnergyExcess.h
 *
 * \ingroup LowEPlots
 *
 * \brief Class def header for a class ERAnaLowEnergyExcess
 *
 * @author kaleko
 */

/** \addtogroup LowEPlots

    @{*/

#ifndef ERTOOL_ERAnaLowEnergyExcess_H
#define ERTOOL_ERAnaLowEnergyExcess_H

#include "ERTool/Base/AnaBase.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include <string>
#include "DataFormat/storage_manager.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mceventweight.h"
// #include "ERToolBackend/ParticleID.h"
#include "../../LArLiteApp/fluxRW/fluxRW.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
#include "LEERW.h"
#include <cmath>
#include "GeoAlgo/GeoAlgo.h"
#include "ECCQECalculator.h"


namespace ertool {

    /**
       \class ERAnaLowEnergyExcess
       User custom Analysis class made by kazuhiro
     */
    class ERAnaLowEnergyExcess : public AnaBase {

    public:

        /// Default constructor
        ERAnaLowEnergyExcess(const std::string& name = "ERAnaLowEnergyExcess");

        /// Default destructor
        virtual ~ERAnaLowEnergyExcess() {}

        /// Reset function
        virtual void Reset() {}

        /// Function to accept fclite::PSet
        void AcceptPSet(const ::fcllite::PSet& cfg) {}

        /// Called @ before processing the first event sample
        void ProcessBegin();

        /// Function to evaluate input showers and determine a score
        bool Analyze(const EventData &data, const ParticleGraph &ps);

        /// Called after processing the last event sample
        void ProcessEnd(TFile* fout);

        /// setting result tree name for running the LowEExcess plotting code
        void SetTreeName(const std::string& name) { _treename = name; }

        void SetStorageManager(larlite::storage_manager* mgr) { _storage_manager = mgr; }

        //Set this to true if you're running over LEE sample (so it uses LEE reweighting package)
        void SetLEESampleMode(bool flag) { _LEESample_mode = flag; }
        // File with LEERW relevant histograms stored in it
        void SetLEEFilename(const std::string& name)     { _LEE_filename = name; }
        void SetLEENEvents(size_t n_evts_passing_filter) { _LEE_evts_passing_filter = n_evts_passing_filter; }
        void SetLEECorrHistName(const std::string& name) { _LEE_corrhist_name = name; }

    private:

        // Calc new E_nu^calo, with missing pT cut
        double EnuCaloMissingPt(const std::vector< ::ertool::NodeID_t >& Children, const ParticleGraph &graph);

        // Determine if the event is "simple" (1e, np, 0else)
        bool isInteractionSimple(const Particle &singleE, const ParticleGraph &ps, const EventData &data);

        /// Function to compute BITE relevant variables (in ttree) and fill them
        void FillBITEVariables(const Shower &singleE_shower, const Particle &p);

        /// Function to compute BNB flux RW weight, or LEE weight (if in LEE mode)
        double GetWeight(const ParticleGraph mc_graph);

        // keep movin'
	std::vector<double> GetMCWeights();

        /// Function to compute various neutrino energy definitions and fill them
        void FillRecoNuEnergies(const Particle &nue, const ParticleGraph &ps, const EventData &data);

        /// Function to sum energy for all tracks/showers coming within 5cm of neutrino vertex
        /// excluding the "singleE" energy
        void FillVertexEnergy(const Particle &nue, const ParticleGraph &ps, const EventData &data);

        // Result tree comparison for reconstructed events
        TTree* _result_tree;
        std::string _treename;

        double _e_nuReco;         /// Neutrino energy
        double _e_dep;            /// Neutrino energy
        double _weight;
        std::vector<double> _mcweight;
        int _parentPDG;           /// true PDG of parent of the electron (only for running on MC)
        int _ptype;               /// neutrino ptype to further break down nue slice in stacked histo
        int _mcPDG;               /// true PDG of "single electron" (probably 11 or 22)
        int _mcGeneration;        /// True generation of single electron (to characterize cosmics and other backgrounds)
        double _longestTrackLen;  /// longest track associated with the reconstructed neutrino
        double _x_vtx;            /// Neutrino vertex points (x,y,z separated)
        double _y_vtx;
        double _z_vtx;
        double _e_theta;          /// Electron's angle w.r.t/ z- axis
        double _e_phi;            /// Electron's phi angle
        double _e_Edep;           /// Electron's truth energy
        double _e_CCQE;           /// Electron's CCQE energy
        double _nu_p;             /// Neutrino reconstructed momentum magnitude
        double _nu_pt;            /// Component of nu momentum that is transverse (_nu_p*sin(_nu_theta))
        double _nu_theta;         /// Neutrino's reconstructed angle w.r.t. z- axis
        int _n_children;          /// Number of children associated with the neutrino interaction
        bool _is_simple;          /// Whether the interaction is 1e+np+0else (reconstructed)
        double _dedx;             /// dedx of "single electron" shower
        double _flash_time;       /// opflash associated with electron... flash time
        double _summed_flash_PE;  /// total reconstructed PE of the flash
        bool _maybe_pi0_MID;      /// whether the neutrino has a gamma tagged as one of its children
        int _n_ertool_showers;    /// total number of ertool::Showers in the event
        int _n_nues_in_evt;       /// # of nues reconstructed in the entire event (pi0 evts sometimes have two  )
        bool _has_muon_child;     /// If there is a muon associated with the reconstructed nue
        double _e_nuReco_better;    /// trying a better definition of energy
        double _vertex_energy; /// Summed energy of all things passing within 5cm of vertex, excluding the electron
        int _mc_origin;           /// mctruth/mctrack/mcshower Origin (==2 if from a cosmic)
        double _mc_time;          /// ertool::Shower._time for the single electron
        double _mc_nu_energy;     /// true neutrino energy if there is a neutrino
        double _trigger_hack_time; /// randomly selected cosmic track arrival time

        
        // prepare TTree with variables
        void PrepareTreeVariables();
        /// Function to re-set TTree variables
        void ResetTreeVariables();

        ::fluxRW _fluxRW;

        // ertool_helper::ParticleID singleE_particleID;
        ertool::Shower singleE_shower;

        bool _LEESample_mode = false;

        // Variables for B.I.T.E analysis
        double _dist_2wall_shr ;  /// Electron shower backwards distance 2 wall
        double _dist_2wall_vtx;   /// Vertex backwards distance 2 wall
        double _dist_2wall_longz_shr;
        double _perp_dist2wall_shr; ///e Shower's cloest perpendicular distance to TPC wall
        double _perp_dist2wall_vtx; ///Vertex's   cloest perpendicular distance to TPC wall
        ::geoalgo::AABox _vactive;
        ::geoalgo::AABox _vactive_longz;
        ::geoalgo::Sphere _vtx_sphere;

        std::string _LEE_filename = "";
        size_t _LEE_evts_passing_filter = 0;
        std::string _LEE_corrhist_name = "";

    protected:
        larlite::storage_manager* _storage_manager;

        ::lee::LEERW _rw;
        ::geoalgo::GeoAlgo _geoalg;
        ::lee::util::ECCQECalculator _eccqecalc;

    };
}
#endif

/** @} */ // end of doxygen group
