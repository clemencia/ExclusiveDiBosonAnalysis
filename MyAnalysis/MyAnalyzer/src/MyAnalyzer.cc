#include"MyAnalysis/MyAnalyzer/interface/MyAnalyzer.h"

//
// constructors and destructor
//
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)

{

  isMC = iConfig.getParameter<bool>("ismc"); 

   //now do what ever initialization is needed
  muon_id        = new vector<bool>;
  muon_pt        = new vector<double>;
  muon_E         = new vector<double>;
  muon_vtx_idx   = new vector<int>;
  muon_px        = new vector<double>;
  muon_py        = new vector<double>;
  muon_pz        = new vector<double>;
  muon_Q         = new vector<double>;
  muon_et        = new vector<double>;
  muon_eta       = new vector<double>;
  muon_phi       = new vector<double>;

  electron_id      = new vector<bool>;
  electron_pt      = new vector<double>;
  electron_eta     = new vector<double>;
  electron_phi     = new vector<double>;
  electron_E       = new vector<double>;
  electron_vtx_idx = new vector<int>;
  electron_px      = new vector<double>;
  electron_py      = new vector<double>;
  electron_pz      = new vector<double>;
  electron_Q       = new vector<double>;
  electron_et      = new vector<double>;

  vertex_x        = new vector<double>;
  vertex_y        = new vector<double>;
  vertex_z        = new vector<double>;
  vertex_idx      = new vector<int>;
  vertex_ntracks  = new vector<int>;
  vertex_mumu_idx = new vector<int>;
  vertex_ee_idx   = new vector<int>;
  vertex_emu_idx  = new vector<int>;

  n_tracks_per_vtx= new vector<int>;
  track_pt        = new vector<vector<double>* >;
  track_eta       = new vector<vector<double>* >;
  track_phi       = new vector<vector<double>* >;
  track_E         = new vector<vector<double>* >;
  track_vtx_idx   = new vector<vector<int>* >;
  track_px        = new vector<vector<double>* >;
  track_py        = new vector<vector<double>* >;
  track_pz        = new vector<vector<double>* >;
  track_Q         = new vector<vector<double>* >;

  proton_id        = new vector<bool>;
  proton_pt        = new vector<double>;
  proton_eta       = new vector<double>;
  proton_phi       = new vector<double>;
  proton_E         = new vector<double>;
  proton_vtx_idx   = new vector<int>;
  proton_px        = new vector<double>;
  proton_py        = new vector<double>;
  proton_pz        = new vector<double>;
  proton_Q         = new vector<double>;
  proton_et        = new vector<double>;

  genpart_pdgID    = new vector<int>;
  genpart_pt       = new vector<double>;
  genpart_eta      = new vector<double>;
  genpart_phi      = new vector<double>;
  genpart_E        = new vector<double>;
  genpart_px       = new vector<double>;
  genpart_py       = new vector<double>;
  genpart_pz       = new vector<double>;
  genpart_Q        = new vector<double>;
  genpart_et       = new vector< double>;

  trigger_ps       = new vector<int>;
  trigger_name     = new vector<string>;
  trigger_decision = new vector<int>;

  newtree=fs->make<TTree>("NewNtuple","NewNtuple");
  
  h_pt_extra_track=fs->make<TH1F>("h_pt_extra_track","Extra Tracks; number of tracks;  events",5000, 0., 5000.);


  newtree->Branch("RunNo", &RunNo,"RunNo/i");
  newtree->Branch("iEvt",  &iEvt,"iEvt/i");
  newtree->Branch("LumiB", &LumiB,"LumiB/i");
  newtree->Branch("PileUpWeight",&PileUpWeight,"PileUpWeight/d");

  newtree->Branch("MET",    &MET,"MET/d");
  newtree->Branch("MET_phi",&MET_phi,"MET_phi/d");
  newtree->Branch("MET_px", &MET_px,"MET_px/d");
  newtree->Branch("MET_py", &MET_py,"MET_py/d");
  newtree->Branch("SumET",  &SumET,"SumET/d");

  newtree->Branch("n_muons",&n_muons,"n_muons/i");
  newtree->Branch("muon_id",&muon_id);
  newtree->Branch("muon_pt",&muon_pt);
  newtree->Branch("muon_px",&muon_px);
  newtree->Branch("muon_py",&muon_py);
  newtree->Branch("muon_pz",&muon_pz);
  newtree->Branch("muon_et",&muon_et);
  newtree->Branch("muon_E",&muon_E);
  newtree->Branch("muon_Q",&muon_Q);
  newtree->Branch("muon_eta",&muon_eta);
  newtree->Branch("muon_phi",&muon_phi);
  newtree->Branch("muon_vtx_idx",&muon_vtx_idx);

  newtree->Branch("n_electrons",&n_electrons,"n_electrons/i");
  newtree->Branch("electron_id",&electron_id);
  newtree->Branch("electron_pt",&electron_pt);
  newtree->Branch("electron_eta",&electron_eta);
  newtree->Branch("electron_phi",&electron_phi);
  newtree->Branch("electron_E",&electron_E);
  newtree->Branch("electron_vtx_idx",&electron_vtx_idx);
  newtree->Branch("electron_px",&electron_px);
  newtree->Branch("electron_py",&electron_py);
  newtree->Branch("electron_pz",&electron_pz);
  newtree->Branch("electron_Q",&electron_Q);
  newtree->Branch("electron_et",&electron_et);

  newtree->Branch("n_protons",&n_protons,"n_protons/i");
  newtree->Branch("proton_id",&proton_id);
  newtree->Branch("proton_pt",&proton_pt);
  newtree->Branch("proton_eta",&proton_eta);
  newtree->Branch("proton_phi",&proton_phi);
  newtree->Branch("proton_E",&proton_E);
  newtree->Branch("proton_vtx_idx",&proton_vtx_idx);
  newtree->Branch("proton_px",&proton_px);
  newtree->Branch("proton_py",&proton_py);
  newtree->Branch("proton_pz",&proton_pz);
  newtree->Branch("proton_Q",&proton_Q);
  newtree->Branch("proton_et",&proton_et);


  newtree->Branch("n_vertices",&n_vertices,"n_vertices/i");
  newtree->Branch("vertex_ntracks",&vertex_ntracks);
  newtree->Branch("vertex_mumu_idx",&vertex_mumu_idx);
  newtree->Branch("vertex_ee_idx",&vertex_ee_idx);
  newtree->Branch("vertex_emu_idx",&vertex_emu_idx);
  newtree->Branch("vertex_idx",&vertex_idx);
  newtree->Branch("vertex_x",&vertex_x);
  newtree->Branch("vertex_y",&vertex_y);
  newtree->Branch("vertex_z",&vertex_z);


  newtree->Branch("n_tracks_per_vtx",&n_tracks_per_vtx,"n_tracks_per_vtx/i");
  newtree->Branch("track_pt",&track_pt);
  newtree->Branch("track_eta",&track_eta);
  newtree->Branch("track_phi",&track_phi);
  newtree->Branch("track_E",&track_E);
  newtree->Branch("track_vtx_idx",&track_vtx_idx);
  newtree->Branch("track_px",&track_px);
  newtree->Branch("track_py",&track_py);
  newtree->Branch("track_pz",&track_pz);
  newtree->Branch("track_Q",&track_Q);

  newtree->Branch("n_genparticles",&n_genparticles,"n_genparticles/i");
  newtree->Branch("genpart_pdgID",&genpart_pdgID);
  newtree->Branch("genpart_pt",&genpart_pt);
  newtree->Branch("genpart_eta",&genpart_eta);
  newtree->Branch("genpart_phi",&genpart_phi);
  newtree->Branch("genpart_E",&genpart_E);
  newtree->Branch("genpart_px",&genpart_px);
  newtree->Branch("genpart_py",&genpart_py);
  newtree->Branch("genpart_pz",&genpart_pz);
  newtree->Branch("genpart_Q",&genpart_Q);
  newtree->Branch("genpart_et",&genpart_et);

  newtree->Branch("trigger_ps",      &trigger_ps);
  newtree->Branch("trigger_name",    &trigger_name);
  newtree->Branch("trigger_decision",&trigger_decision);


  newtree->Branch("mllnunu",&mllnunu,"mllnunu/d");
  newtree->Branch("sqrtsgamgam",&sqrtsgamgam,"sqrtsgamgam/d");



}

MyAnalyzer::~MyAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete muon_id        ;
  delete muon_pt        ;
  delete muon_E         ;
  delete muon_vtx_idx   ;
  delete muon_px        ;
  delete muon_py        ;
  delete muon_pz        ;
  delete muon_Q         ;
  delete muon_et        ;
  delete muon_eta       ;
  delete muon_phi       ;

  delete electron_id      ;
  delete electron_pt      ;
  delete electron_eta     ;
  delete electron_phi     ;
  delete electron_E       ;
  delete electron_vtx_idx ;
  delete electron_px      ;
  delete electron_py      ;
  delete electron_pz      ;
  delete electron_Q       ;
  delete electron_et      ;

  delete vertex_x        ;
  delete vertex_y        ;
  delete vertex_z        ;
  delete vertex_idx      ;
  delete vertex_ntracks  ;
  delete vertex_mumu_idx ;
  delete vertex_ee_idx   ;
  delete vertex_emu_idx  ;

  delete n_tracks_per_vtx;
  delete track_pt        ;
  delete track_eta       ;
  delete track_phi       ;
  delete track_E         ;
  delete track_vtx_idx   ;
  delete track_px        ;
  delete track_py        ;
  delete track_pz        ;
  delete track_Q         ;

  delete proton_id        ;
  delete proton_pt        ;
  delete proton_eta       ;
  delete proton_phi       ;
  delete proton_E         ;
  delete proton_vtx_idx   ;
  delete proton_px        ;
  delete proton_py        ;
  delete proton_pz        ;
  delete proton_Q         ;
  delete proton_et        ;

  delete genpart_pdgID    ;
  delete genpart_pt       ;
  delete genpart_eta      ;
  delete genpart_phi      ;
  delete genpart_E        ;
  delete genpart_px       ;
  delete genpart_py       ;
  delete genpart_pz       ;
  delete genpart_Q        ;
  delete genpart_et       ;

  delete trigger_ps       ;
  delete trigger_name     ;
  delete trigger_decision ;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ////////////////////////////////////////////////////////////////////////////////
   //Fill basic event variables
   ////////////////////////////////////////////////////////////////////////////////
   //Event information
   RunNo = iEvent.id().run();
   iEvt  = iEvent.id().event();
   LumiB = iEvent.luminosityBlock();
   PileUpWeight=1.;
   int nextra_tracks = -999;
   int nextra_tracks_mumu = -999;
   int nextra_tracks_ee = -999;


 
   ////////////////////////////////////////////////////////////////////////////////
   //        Fill MET info
   ////////////////////////////////////////////////////////////////////////////////

   Handle< edm::View<reco::PFMET> > metHandle;
   iEvent.getByLabel("pfType1CorrectedMet",metHandle);
   edm::View<reco::PFMET> mets = *metHandle;

   //edm::View<reco::PFMET>::const_iterator iMet;
   // no need to loop, there's only one met for event
   //for (iMet = mets.begin(); iMet != mets.end(); ++iMet) {

   reco::PFMET iMet = mets.front();
   MET    = iMet.et();
   MET_px = iMet.px();
   MET_py = iMet.py();
   MET_phi= iMet.phi();
   SumET  = iMet.sumEt();
   //}

   ////////////////////////////////////////////////////////////////////////////////
   //         Fill Trigger info
   ////////////////////////////////////////////////////////////////////////////////
   Handle<TriggerResults> hltResults;
   iEvent.getByLabel(InputTag("TriggerResults","","HLT"),hltResults);
   const TriggerNames & trigNames = iEvent.triggerNames(*hltResults);

   for(unsigned int i=0; i<trigNames.size();i++){
     (*trigger_ps).push_back(hltConfig_.prescaleValue(iEvent,iSetup,trigNames.triggerName(i)));
     (*trigger_name).push_back(trigNames.triggerName(i));
     (*trigger_decision).push_back(hltResults->accept(i));
   }

   if(isMC){

     ////////////////////////////////////////////////////////////////////////////////
     // for MC fill Gen particle variables
     ////////////////////////////////////////////////////////////////////////////////
     
     Handle<reco::GenParticleCollection> genP;
     iEvent.getByLabel("genParticles",genP);

     //Loop over gen particles with pt and eta cuts
     reco::GenParticleCollection::const_iterator mcIter = 0 ;
     int nphotons=0;
     double pho1px,pho1py,pho1pz,pho1E;
     double pho2px,pho2py,pho2pz,pho2E;
     math::XYZLorentzVector fourp_pho1,fourp_pho2;
     for (mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
       if (mcIter->pdgId()==22 && mcIter->status==-1){ // status of incoming photons? was -1 in LHE file
	 nphotons++;
	 if (nphotons<2){
	   pho1px = mcIter->px();
	   pho1py = mcIter->py();
	   pho1pz = mcIter->pz();
	   pho1E = mcIter->energy();
	   fourp_pho1.SetPxPyPxE(pho1px,pho1py,pho1pz,pho1E);
	 }
	 else if (nphotons ==2){
	   pho2px = mcIter->px();
	   pho2py = mcIter->py();
	   pho2pz = mcIter->pz();
	   pho2E = mcIter->energy();
	   fourp_pho2.SetPxPyPxE(pho2px,pho2py,pho2pz,pho2E);
	   sqrtqgamgam = (fourp_pho1+fourp_pho2).M();
	   ///make  Center of Mass energy between 2 photons ?e1+e2? or inv mass?
	 }
	 else{
	   cout<<"More than two photons!"<<endl;
	 }
       }// end incoming photons

       if (mcIter->pt()>20 && fabs(mcIter->eta())<2.5){ 
	 (*genpart_pt).push_back(mcIter->pt());
	 (*genpart_px).push_back(mcIter->px());
	 (*genpart_py).push_back(mcIter->py());
	 (*genpart_pz).push_back(mcIter->pz());
	 (*genpart_et).push_back(mcIter->et());
	 (*genpart_E).push_back(mcIter->energy());
	 (*genpart_Q).push_back(mcIter->charge());
	 (*genpart_eta).push_back(mcIter->eta());
	 (*genpart_phi).push_back(mcIter->phi());
	 (*genpart_pdgID).push_back(mcIter->pdgId());
       }	   
     }//end of looking at gen particles
   }//end of specifying that it needs to be MC

}//end analyze



// ------------ method called once each job just before starting event loop  ------------
void 
MyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MyAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyAnalyzer);
