/**
 * \file ERAnaLowEnergyExcess.h
 *
 * \ingroup LowEPlots
 * 
 * \brief Class def header for a class ERAnaLowEnergyExcess
 *
 * @author jzennamo
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
#include "DataFormat/mctruth.h"
// #include "ERToolBackend/ParticleID.h"
#include "../../LArLiteApp/fluxRW/fluxRW.h"
#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
#include "LEERW.h"
#include <cmath>
 #include "GeoAlgo/GeoAlgo.h"


namespace ertool {

  /**
     \class ERAnaLowEnergyExcess
     User custom Analysis class made by kazuhiro
   */
  class ERAnaLowEnergyExcess : public AnaBase {
  
  public:

    /// Default constructor
    ERAnaLowEnergyExcess(const std::string& name="ERAnaLowEnergyExcess");

    /// Default destructor
    virtual ~ERAnaLowEnergyExcess(){}

    /// Reset function
    virtual void Reset() {}

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg){}

    /// Called @ before processing the first event sample
    void ProcessBegin() {};

    /// Function to evaluate input showers and determine a score
    bool Analyze(const EventData &data, const ParticleGraph &ps);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout);

    /// setting result tree name for running the LowEExcess plotting code
    void SetTreeName(const std::string& name){ _treename = name; }

    /// set the energy cut to be used when counting particles
    void SetECut(double c) { _eCut = c; }

    // geoalgo::AABox TPC;

    //Set this to true if you're running over LEE sample (so it uses LEE reweighting package)
    void SetLEESampleMode(bool flag) { _LEESample_mode = flag; }
    
  private:
    
    // Calc new E_nu^calo, with missing pT cut     
    double EnuCaloMissingPt(const std::vector< ::ertool::NodeID_t >& Children, const ParticleGraph &graph);

    // Result tree comparison for reconstructed events
    TTree* _result_tree;
    std::string _treename;

    float _eCut;

    double _e_nuReco;         /// Neutrino energy
    double _e_dep;            /// Neutrino energy
    double _weight;
    int _numEvts;
    bool _is_fiducial;
    int _parentPDG;           /// true PDG of parent of the electron (only for running on MC)
    int _mcPDG;               /// true PDG of "single electron" (probably 11 or 22)
    int _mcGeneration;        /// True generation of single electron (to characterize cosmics and other backgrounds)
    double _longestTrackLen;  /// longest track associated with the reconstructed neutrino
    double _x_vtx;            /// Neutrino vertex points (x,y,z separated)
    double _y_vtx;
    double _z_vtx;
    double _e_theta;          /// Electron's angle w.r.t/ z- axis
    double _e_phi;            /// Electron's phi angle
    double _nu_p;             /// Neutrino reconstructed momentum magnitude
    double _nu_pt;            /// Component of nu momentum that is transverse (_nu_p*sin(_nu_theta))
    double _nu_theta;         /// Neutrino's reconstructed angle w.r.t. z- axis


    // prepare TTree with variables
    void PrepareTreeVariables();
    /// Function to re-set TTree variables
    void ResetTreeVariables();

    ::fluxRW _fluxRW;

    // ertool_helper::ParticleID singleE_particleID;
    ertool::Shower singleE_shower;

    bool _LEESample_mode = false;

protected:

  ::lee::LEERW _rw;
  ::geoalgo::GeoAlgo _geoalg;

  };
}
#endif

/** @} */ // end of doxygen group 
