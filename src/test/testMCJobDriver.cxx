//____________________________________________________________________________
/*!

\program testMCJobDriver

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 22, 2005
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/XSecSplineList.h"

using std::string;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

void GetCommandLineArgs(int argc, char ** argv);

//command line options
bool   gOptBuildSplines; // spline building option
string gOptRootGeom;     // detector geometry ROOT file

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- Parse command line arguments
  GetCommandLineArgs(argc, argv);

  //-- Create the GENIE MC-job driver
  GMCJDriver mcj;

  //-- load and/or build splines if required
  XSecSplineList * xssl = 0;
  if(gOptBuildSplines) {
     xssl = XSecSplineList::Instance();
     // check whether there is a spline-list XML file to load
     string spllst_load_xmlfile =
             (gSystem->Getenv("GSPLOAD") ? gSystem->Getenv("GSPLOAD") : "");
     LOG("test", pINFO) << "$GSPLOAD env.var = " << spllst_load_xmlfile;

     if(spllst_load_xmlfile.size()>0) {
       LOG("test", pINFO) << "Loading cross section splines from an xml file";
       XmlParserStatus_t status = xssl->LoadFromXml(spllst_load_xmlfile);
       assert(status==kXmlOK);
     }
     // create any spline that is needed but is not loaded
     mcj.UseSplines();
  }

  //-- Specify a flux driver
/*
  LOG("Main", pINFO)  << "Creating [GFlukaAtmo3DFlux] flux driver";

  GFlukaAtmo3DFlux * flux = new GFlukaAtmo3DFlux;

  flux->SetNuMuFluxFile("/home/costas/tmp/downloads/sdave_numu07.dat");
  flux->SetNuMuBarFluxFile("/home/costas/tmp/downloads/sdave_anumu07.dat");
  flux->SetNuEFluxFile("/home/costas/tmp/downloads/sdave_nue07.dat");
  flux->SetNuEBarFluxFile("/home/costas/tmp/downloads/sdave_anue07.dat");
  flux->SetRadii(1000.,100.);

  flux->LoadFluxData();
*/

  LOG("Main", pINFO)  << "Creating [GCylindTH1Flux] flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  TF1 * f1 = new TF1("f1","1./x+2.",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu spectrum", 20,0.5,5);
  spectrum1->FillRandom("f1",10000);

  TVector3 direction(0,0,1);
  TVector3 beam_spot(0,0,-10);

  flux -> SetNuDirection      (direction);
  flux -> SetBeamSpot         (beam_spot);
  flux -> SetTransverseRadius (0.5);
  flux -> AddEnergySpectrum   (kPdgNuMu, spectrum1);

  //-- Specify the geometry analyzer

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptRootGeom);

  //-- Set the flux and the geometry analyzer to the GENIE MC driver

  LOG("Main", pINFO)
    << "Creating the GENIE MC Job Driver & specifying flux & geometry";

  GFluxI *        fluxb = dynamic_cast<GFluxI *>       (flux);
  GeomAnalyzerI * geomb = dynamic_cast<GeomAnalyzerI *>(geom);

  mcj.UseFluxDriver  (fluxb);
  mcj.UseGeomAnalyzer(geomb);

  //-- Configure the GENIE MC driver

  mcj.Configure();

  //-- Start generating events -here, just 1 for testing purposes-

  EventRecord * event = mcj.GenerateEvent();

  LOG("Main", pINFO) << *event;

  delete f1;
  delete flux;
  delete geom;

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  string base_dir = string( gSystem->Getenv("GENIE") );

  // default options
  gOptBuildSplines = false;  
  gOptRootGeom     = base_dir + string("/src/test/TestGeometry.root");

  char * argument = new char[128];

  while( argc>1 && (argv[1][0] == '-'))
  {
    if (argv[1][1] == 'f') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptRootGeom = string(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptRootGeom = string(argument);
      }
    }

    if (argv[1][1] == 's') gOptBuildSplines = true;

    argc--;
    argv++;
  }

  delete [] argument;
}
//___________________________________________________________________

