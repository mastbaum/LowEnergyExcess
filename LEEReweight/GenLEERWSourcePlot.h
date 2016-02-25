/**
 * \file GenLEERWSourcePlot.h
 *
 * \ingroup LEEReweight
 * 
 * \brief Class def header for a class GenLEERWSourcePlot
 *
 * @author davidkaleko
 */

/** \addtogroup LEEReweight

    @{*/

#ifndef LARLITE_GENLEERWSOURCEPLOT_H
#define LARLITE_GENLEERWSOURCEPLOT_H

#include "Analysis/ana_base.h"
#include "TH2.h"

namespace larlite {
  /**
     \class GenLEERWSourcePlot
     User custom analysis class made by SHELL_USER_NAME
   */
  class GenLEERWSourcePlot : public ana_base{
  
  public:

    /// Default constructor
    GenLEERWSourcePlot(){ _name="GenLEERWSourcePlot"; _fout=0; _plot_name = ""; _out_hist = 0;}

    /// Default destructor
    virtual ~GenLEERWSourcePlot(){}

    /** IMPLEMENT in GenLEERWSourcePlot.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GenLEERWSourcePlot.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GenLEERWSourcePlot.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setPlotName(const std::string plotname) { _plot_name = plotname; }

  protected:
    
    std::string _plot_name;
    TH2D* _out_hist;
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
