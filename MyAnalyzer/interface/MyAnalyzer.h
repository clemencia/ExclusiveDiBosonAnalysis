// -*- C++ -*-
//
// Package:    MyAnalyzer
// Class:      MyAnalyzer
//
/*
*\class MyAnalyzer MyAnalyzer.cc MyAnalyzer/src/MyAnalyzer.cc
    Description: [one line class summary]
 Implementation: [Notes on implementation] 
                                                                             */
//
// Original Author:  Clemencia Mora Herrera,32 2-A13,+41227676740,
//         Created:  Wed May 29 22:56:42 CEST 2013
// $Id$
//
//
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"


#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration
//

using namespace std;


class MyAnalyzer : public edm::EDAnalyzer {
public:
  explicit MyAnalyzer(const edm::ParameterSet&);
  ~MyAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------                                                                                                                       

  TTree * newtree;
  int iEvt;
  int RunNo;
  int LumiB;
  double PileUpWeight;


  bool isMC;
  std::vector<edm::InputTag> isoValInputTags;

  HLTConfigProvider hltConfig_;

  TH1F * h_pt_extra_track;
  edm::LumiReWeighting *LumiWeights;

  int n_muons;
  vector<bool>*   muon_id;
  vector<double>* muon_pt;
  vector<double>* muon_eta;
  vector<double>* muon_phi;
  vector<double>* muon_E;
  vector<int>*    muon_vtx_idx;
  vector<double>* muon_px;
  vector<double>* muon_py;
  vector<double>* muon_pz;
  vector<double>* muon_Q;
  vector<double>* muon_et;

  int n_electrons;
  vector<bool>* electron_id;
  vector<double>* electron_pt;
  vector<double>* electron_eta;
  vector<double>* electron_phi;
  vector<double>* electron_E;
  vector<int>*    electron_vtx_idx;
  vector<double>* electron_px;
  vector<double>* electron_py;
  vector<double>* electron_pz;
  vector<double>* electron_Q;
  vector<double>* electron_et;

  double MET;
  double MET_phi;
  double MET_px;
  double MET_py;
  double SumET;

  int n_vertices;
  vector<double>* vertex_x;
  vector<double>* vertex_y;
  vector<double>* vertex_z;
  vector<int>*    vertex_idx; // so far unused
  vector<int>*    vertex_extra_ntracks_ee;
  vector<int>*    vertex_extra_ntracks_mumu;
  vector<int>*    vertex_extra_ntracks_emu;
  vector<int>*    vertex_mumu_cand1_idx;
  vector<int>*    vertex_mumu_cand2_idx;
  vector<int>*    vertex_ee_cand1_idx;
  vector<int>*    vertex_ee_cand2_idx;
  vector<int>*    vertex_emu_candE_idx;
  vector<int>*    vertex_emu_candMu_idx;


  vector<int>* n_tracks_per_vtx;
  vector<vector<double> >* track_pt;
  vector<vector<double> >* track_eta;
  vector<vector<double> >* track_phi;
  vector<vector<double> >* track_E;
  vector<vector<int> >*    track_vtx_idx;
  vector<vector<double> >* track_px;
  vector<vector<double> >* track_py;
  vector<vector<double> >* track_pz;
  vector<vector<double> >* track_Q;

  int n_protons;
  vector<bool>*   proton_id;
  vector<double>* proton_pt;
  vector<double>* proton_eta;
  vector<double>* proton_phi;
  vector<double>* proton_E;
  vector<int>*    proton_vtx_idx;
  vector<double>* proton_px;
  vector<double>* proton_py;
  vector<double>* proton_pz;
  vector<double>* proton_Q;
  vector<double>* proton_et;

  double mllnunu;
  double sqrtsgamgam;

  int n_genparticles;
  vector<int>*   genpart_pdgID;
  vector<double>* genpart_pt;
  vector<double>* genpart_eta;
  vector<double>* genpart_phi;
  vector<double>* genpart_E;
  vector<double>* genpart_px;
  vector<double>* genpart_py;
  vector<double>* genpart_pz;
  vector<double>* genpart_Q;
  vector<double>* genpart_et;

  vector<int>*    trigger_ps;
  vector<string>* trigger_name;
  vector<int>*    trigger_decision;


 };
