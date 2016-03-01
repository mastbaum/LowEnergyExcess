/**
 * \file ERAlgoTagEmulatedDeletionsCosmic.h
 *
 * \ingroup ERAnalysis
 * 
 * \brief Class def header for a class ERAlgoTagEmulatedDeletionsCosmic
 *
 * @author davidkaleko
 */

/** \addtogroup ERAnalysis

    @{*/

#ifndef ERTOOL_ERALGOTAGEMULATEDDELETIONSCOSMIC_H
#define ERTOOL_ERALGOTAGEMULATEDDELETIONSCOSMIC_H

#include "ERTool/Base/AlgoBase.h"

namespace ertool {

  /**
     \class ERAlgoTagEmulatedDeletionsCosmic
     User custom Algorithm class made by kazuhiro
   */
  class ERAlgoTagEmulatedDeletionsCosmic : public AlgoBase {
  
  public:

    /// Default constructor
    ERAlgoTagEmulatedDeletionsCosmic(const std::string& name="ERAlgoTagEmulatedDeletionsCosmic");

    /// Default destructor
    virtual ~ERAlgoTagEmulatedDeletionsCosmic(){};

    /// Reset function
    void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Reconstruct(const EventData &data, ParticleGraph& graph);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);

  };
}
#endif

/** @} */ // end of doxygen group 
