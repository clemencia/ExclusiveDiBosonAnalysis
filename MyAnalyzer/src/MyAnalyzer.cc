#include"ExclusiveDiBosonAnalysis/MyAnalyzer/interface/MyAnalyzer.h"

//
// constructors and destructor
//
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)

{

  isMC = iConfig.getParameter<bool>("ismc"); 
  isoValInputTags = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

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
  vertex_idx      = new vector<int>; //so far unused
  vertex_extra_ntracks_ee  = new vector<int>;
  vertex_extra_ntracks_mumu  = new vector<int>;
  vertex_extra_ntracks_emu  = new vector<int>;
  vertex_mumu_cand1_idx = new vector<int>;
  vertex_mumu_cand2_idx = new vector<int>;
  vertex_ee_cand1_idx   = new vector<int>;
  vertex_ee_cand2_idx   = new vector<int>;
  vertex_emu_candE_idx   = new vector<int>;
  vertex_emu_candMu_idx  = new vector<int>;

  n_tracks_per_vtx= new vector<int>;
  track_pt        = new vector<vector<double> >;
  track_eta       = new vector<vector<double> >;
  track_phi       = new vector<vector<double> >;
  track_E         = new vector<vector<double> >;
  track_vtx_idx   = new vector<vector<int> >;
  track_px        = new vector<vector<double> >;
  track_py        = new vector<vector<double> >;
  track_pz        = new vector<vector<double> >;
  track_Q         = new vector<vector<double> >;

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

  edm::Service<TFileService> fs;

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

  //get generated ougoing protons?
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
  newtree->Branch("vertex_extra_ntracks_ee",&vertex_extra_ntracks_ee);
  newtree->Branch("vertex_extra_ntracks_mumu",&vertex_extra_ntracks_mumu);
  newtree->Branch("vertex_extra_ntracks_emu",&vertex_extra_ntracks_emu);
  newtree->Branch("vertex_mumu_cand1_idx",&vertex_mumu_cand1_idx);
  newtree->Branch("vertex_mumu_cand2_idx",&vertex_mumu_cand2_idx);
  newtree->Branch("vertex_ee_cand1_idx",&vertex_ee_cand1_idx);
  newtree->Branch("vertex_ee_cand2_idx",&vertex_ee_cand2_idx);
  newtree->Branch("vertex_emu_candE_idx",&vertex_emu_candE_idx);
  newtree->Branch("vertex_emu_candMu_idx",&vertex_emu_candMu_idx);
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
  //maybe add also ptll
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
  delete vertex_extra_ntracks_ee  ;
  delete vertex_extra_ntracks_mumu  ;
  delete vertex_extra_ntracks_emu  ;
  delete vertex_mumu_cand1_idx ;
  delete vertex_mumu_cand2_idx ;
  delete vertex_ee_cand1_idx   ;
  delete vertex_ee_cand2_idx   ;
  delete vertex_emu_candE_idx  ;
  delete vertex_emu_candMu_idx   ;

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
     reco::GenParticleCollection::const_iterator mcIter ;
     int nphotons=0;
     double pho1px,pho1py,pho1pz,pho1E;
     double pho2px,pho2py,pho2pz,pho2E;
     math::XYZTLorentzVector fourp_pho1,fourp_pho2;
     for (mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
       if (mcIter->pdgId()==22 && mcIter->status()==-1){ // status of incoming photons? was -1 in LHE file
	 nphotons++;
	 if (nphotons<2){
	   pho1px = mcIter->px();
	   pho1py = mcIter->py();
	   pho1pz = mcIter->pz();
	   pho1E = mcIter->energy();
	   fourp_pho1.SetPxPyPzE(pho1px,pho1py,pho1pz,pho1E);
	 }
	 else if (nphotons ==2){
	   pho2px = mcIter->px();
	   pho2py = mcIter->py();
	   pho2pz = mcIter->pz();
	   pho2E = mcIter->energy();
	   fourp_pho2.SetPxPyPzE(pho2px,pho2py,pho2pz,pho2E);
	   sqrtsgamgam = (fourp_pho1+fourp_pho2).M();
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


   // beam spot
   edm::Handle<reco::BeamSpot> beamspot_h;
   //   iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
   iEvent.getByLabel("offlineBeamSpot", beamspot_h);
   const reco::BeamSpot &beamSpot = *(beamspot_h.product());


   // vertices
   edm::Handle<reco::VertexCollection> vtxs;
   iEvent.getByLabel("offlinePrimaryVertices", vtxs);
   const reco::VertexCollection* vertexs = vtxs.product();
   reco::VertexCollection::const_iterator vertex_i;


   ////////////////////////////////////////////////////////////////////////////////
   //Muon Selection and filling of ntuple
   ////////////////////////////////////////////////////////////////////////////////

   edm::Handle< edm::View<reco::Muon> > muonHandle;
   iEvent.getByLabel("muons", muonHandle);
   edm::View<reco::Muon> muons = *muonHandle;
   edm::View<reco::Muon>::const_iterator iMuon;

   //loop over muons and store muons that pass id cuts, including pt>20 GeV and |eta|<2.4
   int n_MuonsPassingCuts = 0;
   vector<int> muonCharge;
   vector<double> muonTrackPt;
   vector<double> muonTrackEta;
   vector<double> muonTrackPhi;
   for (iMuon = muons.begin(); iMuon != muons.end(); ++iMuon) {
     bool pass_globalMuon = false;
     if(iMuon->isGlobalMuon())
       pass_globalMuon=true;

     bool pass_pfMuon = false;
     if(iMuon->isPFMuon())
       pass_pfMuon=true;

     bool pass_chi2 = false;
     bool pass_MuonChamberHits = false;
     if(iMuon->globalTrack().isNonnull()){
       if(iMuon->globalTrack()->normalizedChi2()< 10)
	 pass_chi2=true;
       if(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
	 pass_MuonChamberHits=true;
     }

     bool pass_MuonStations = false;
     if(iMuon->numberOfMatchedStations() > 1)
       pass_MuonStations=true;

     bool pass_NPxlHits = false;
     bool pass_NtrackerLayers = false;     
     if(iMuon->innerTrack().isNonnull()){
       if(iMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)
	 pass_NPxlHits = true;
       if(iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)
	 pass_NtrackerLayers = true;
       //       if( fabs(iMuon->innerTrack()->dxy(beamSpot.position())) < 0.2){d0cut = true;}
     }

     bool pass_d0cut = false;
     bool pass_dzcut = false;
     double d0vtx = 999.;
     double dzvtx = 999.;
     if(iMuon->muonBestTrack().isNonnull()){
       if (vtxs->size() > 0) {
	 reco::VertexRef vtx(vtxs, 0);    
	 d0vtx = iMuon->muonBestTrack()->dxy(vtx->position());
	 dzvtx = iMuon->muonBestTrack()->dz(vtx->position());
       } 
       else {
	 d0vtx = iMuon->muonBestTrack()->dxy();
	 dzvtx = iMuon->muonBestTrack()->dz();
       }
       if( fabs(d0vtx) < 0.2)
	 pass_d0cut = true;
       if( fabs(dzvtx) < 0.5)
	 pass_dzcut = true;
     }

     
     if(iMuon->pt()>20&&fabs(iMuon->eta())<2.4){
       bool pass_muonId =false;
       if(pass_d0cut &&	  pass_dzcut &&	  pass_NtrackerLayers &&
	  pass_NPxlHits &&  pass_MuonStations &&  pass_MuonChamberHits &&
	  pass_chi2 &&  pass_pfMuon &&   pass_globalMuon)
	 pass_muonId=true;
       
       
       (*muon_id).push_back(pass_muonId);
       (*muon_pt).push_back(iMuon->pt());
       (*muon_px).push_back(iMuon->px());
       (*muon_py).push_back(iMuon->py());
       (*muon_pz).push_back(iMuon->pz());
       (*muon_et).push_back(iMuon->et());
       (*muon_E).push_back(iMuon->energy());
       (*muon_Q).push_back(iMuon->charge());
       (*muon_eta).push_back(iMuon->eta());
       (*muon_phi).push_back(iMuon->phi());
       if(iMuon->innerTrack().isNonnull()){
	 muonTrackPt.push_back(iMuon->innerTrack()->pt());
	 muonTrackEta.push_back(iMuon->innerTrack()->eta());
	 muonTrackPhi.push_back(iMuon->innerTrack()->phi());
       }
       else{
	 muonTrackPt.push_back(-999.);
	 muonTrackEta.push_back(-999.);
	 muonTrackPhi.push_back(-999.);	 
	 cout<<"There is something wrong with the Muon track"<<endl;
       }
       muonCharge.push_back(iMuon->charge());
       n_MuonsPassingCuts++;
     }//end eta and pt selection

   }// end muons
   n_muons=n_MuonsPassingCuts;

   ////////////////////////////////////////////////////////////////////////////////
   //Electron Selection and filling of ntuple
   ////////////////////////////////////////////////////////////////////////////////

   // New 2012 electron ID variables
   // conversions 
   edm::Handle<reco::ConversionCollection> conversions_h; 
   iEvent.getByLabel("allConversions", conversions_h); 

   // iso deposits 
   std::vector< edm::Handle< edm::ValueMap<double> > >  isoVals(isoValInputTags.size()); 
   for (size_t j = 0; j < isoValInputTags.size(); ++j) { 
     iEvent.getByLabel(isoValInputTags[j], isoVals[j]); 
   } 

   // rho for isolation 
   edm::Handle<double> rhoIso_h; 
   iEvent.getByLabel("kt6PFJetsForIsolation","rho", rhoIso_h); 
   double rhoIso = *(rhoIso_h.product()); 

   // electrons
   edm::Handle<reco::GsfElectronCollection> els_h;
   iEvent.getByLabel("gsfElectrons", els_h);
   
   double checkduplicates[4] = {0,0,0,0};


   unsigned int n = els_h->size();
   int n_ElsPassingCuts = 0;
   vector<int> electronCharge;
   vector<double> electronTrackPt;
   vector<double> electronTrackEta;
   vector<double> electronTrackPhi;

   for(unsigned int i = 0; i < n; ++i) {

     // get reference to electron
     reco::GsfElectronRef ele(els_h, i);

     float pt            = ele->pt();
     float eta           = ele->superCluster()->eta();
     //cout << "ele= " << long(&ele) << endl;
     double iso_ch = (*(isoVals)[0])[ele];
     double iso_em = (*(isoVals)[1])[ele];
     double iso_nh = (*(isoVals)[2])[ele];
     bool pass_mediumEid=false;
     if(PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele ))       
       pass_mediumEid = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, 
						    ele, 
						    conversions_h, 
						    beamSpot, 
						    vtxs, 
						    iso_ch, 
						    iso_em, 
						    iso_nh, 
						    rhoIso);
 

     //     if(passmediumEid&&pt>20&&fabs(eta)<2.4){
     bool electronId = false;
     if(pass_mediumEid)
       electronId = true;

     if(pt>20&&fabs(eta)<2.4){  
       bool duplicate = false;
       if(pt==checkduplicates[0] && 
	  checkduplicates[1]==eta && 
	  checkduplicates[2]== ele->phi()&&
	  checkduplicates[3]==ele->charge()){
	 cout<<"Duplicate Electron"<<endl;
	 duplicate = true;
       }
       else{
	 checkduplicates[0] = pt;
	 checkduplicates[1] = eta;
	 checkduplicates[2] = ele->phi(); 
	 checkduplicates[3] = ele->charge();
       }

       
       if(!duplicate){
	 
	 (*electron_id).push_back(electronId);
	 (*electron_pt).push_back(ele->pt());
	 (*electron_px).push_back(ele->px());
	 (*electron_py).push_back(ele->py());
	 (*electron_pz).push_back(ele->pz());
	 (*electron_et).push_back(ele->et());
	 (*electron_E).push_back(ele->energy());
	 (*electron_Q).push_back(ele->charge());
	 (*electron_eta).push_back(ele->eta());
	 (*electron_phi).push_back(ele->phi());
	 n_ElsPassingCuts++;
	 electronCharge.push_back(ele->charge());
	 if(ele->closestCtfTrackRef().isNonnull()){
	   electronTrackPt.push_back(ele->closestCtfTrackRef()->pt());
	   electronTrackEta.push_back(ele->closestCtfTrackRef()->eta());
	   electronTrackPhi.push_back(ele->closestCtfTrackRef()->phi());
	 }
	 else{
	   electronTrackPt.push_back(-9999.);
	   electronTrackEta.push_back(-9999.);
	   electronTrackPhi.push_back(-9999.);
	 }

       }//end of requiring not a duplicate (this is for FAST SIM)
     }//end of pt and eta requirements
   }// end of loop in electrons
   n_electrons=n_ElsPassingCuts;

   /////////////////////////////////////////////////////////////////////////////////////////////////
   //Store vertex information for events that pass e e criteria (if electrons come from same vertex)
   /////////////////////////////////////////////////////////////////////////////////////////////////
   int count_Nvertices_match_ee = 0;
   int count_Nvertices_match_mumu = 0;
   int count_Nvertices_match_emu = 0;
   for(vertex_i = vertexs->begin(); vertex_i<vertexs->end(); vertex_i++){
     double vertexNtracks = vertex_i->tracksSize();
     double vertexPos_z = vertex_i->z();
     if(fabs(vertexPos_z)<24 ){
       if(vertexNtracks<=75){
	 for(int l = 0;l<n_ElsPassingCuts;l++){
	   for(int k = l+1;k<n_ElsPassingCuts;k++){
       
	     if((electronCharge[l]*electronCharge[k])<0){  
	       bool pass_electron_assoc = false;
	       bool pass_electron2_assoc = false;
	       for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		 if( fabs((*vertex_Tracks)->pt()-electronTrackPt[l])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-electronTrackEta[l])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-electronTrackPhi[l])<0.001){
		   if (pass_electron_assoc)
		     cout<<"this is weird, this electron already had a track associated to it"<<endl;
		   pass_electron_assoc = true;
		   continue; // go to next TRACK
		   //// I do this just in case, I don't want to associate the same track to diff leptons
		 }
		 
		 if( fabs((*vertex_Tracks)->pt()-electronTrackPt[k])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-electronTrackEta[k])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-electronTrackPhi[k])<0.001){
		   if (pass_electron2_assoc)
		     cout<<"this is weird, the electron already had an associated track"<<endl;
		   pass_electron2_assoc = true;
		   continue; // go to next TRACK
		   
		 }

	       }//Loop over tracks to see if the leptons are there
	       
	       if(pass_electron_assoc&&pass_electron2_assoc){
		 count_Nvertices_match_ee++;		 
		 nextra_tracks_ee = vertexNtracks-2;
		 (*vertex_extra_ntracks_ee).push_back(nextra_tracks_ee);
		 (*vertex_ee_cand1_idx).push_back(l);
		 (*vertex_ee_cand2_idx).push_back(k);
		 //		 cout<<"This event passes mu+ mu- criteria and should be vetoed"<<endl;
	       }
	       
	     }//end of if statement for charge
	   }//end loop electron 1


	 }//end loop electron 2

 	 for(int l = 0;l<n_MuonsPassingCuts;l++){
 	   for(int k = l+1;k<n_MuonsPassingCuts;k++){
	     
 	     if((muonCharge[l]*muonCharge[k])<0){
	       //	     if(vertex_ntracks<=17){	     if(vertex_ntracks<=75){
	       
 	       bool pass_muon_assoc = false;
	       bool pass_muon2_assoc = false;
 	       for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		 if( fabs((*vertex_Tracks)->pt()-muonTrackPt[l])<0.001 && 
 		     fabs((*vertex_Tracks)->eta()-muonTrackEta[l])<0.001 && 
 		     fabs((*vertex_Tracks)->phi()-muonTrackPhi[l])<0.001){
 		   pass_muon_assoc = true;
		   /// can I do this? // muon_vtx_idx->push_back(??);
		   continue; // go to next TRACK
		   //// I do this just in case, I don't want to associate the same track to diff leptons

 		 }
		 
 		 if( fabs((*vertex_Tracks)->pt()-muonTrackPt[k])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-muonTrackEta[k])<0.001 && 
 		     fabs((*vertex_Tracks)->phi()-muonTrackPhi[k])<0.001){
 		   pass_muon2_assoc = true;
		   continue; // go to next TRACK
 		 }
 	       }//Loop over tracks to see if the leptons are there
	       
 	       if(pass_muon_assoc && pass_muon2_assoc){
		 count_Nvertices_match_mumu++;
 		 //veto_event=true;		 
 		 nextra_tracks_mumu = vertexNtracks-2;
 		 (*vertex_extra_ntracks_mumu).push_back(nextra_tracks_mumu);
 		 (*vertex_mumu_cand1_idx).push_back(l);
 		 (*vertex_mumu_cand2_idx).push_back(k);
 		 //		 cout<<"This event passes mu+ mu- criteria and should be vetoed"<<endl;
 	       }
	       
 	     }//end of if statement for charge
 	   }//end loop muon1
 	 }//end loop muon2


 	 for(int l = 0;l<n_MuonsPassingCuts;l++){
 	   for(int k = 0;k<n_ElsPassingCuts;k++){
	     
 	     if((muonCharge[l]*electronCharge[k])<0){
	       //	     if(vertex_ntracks<=17){	     if(vertex_ntracks<=75){
	       
 	       bool pass_muon_assoc = false;
	       bool pass_electron_assoc = false;
 	       for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		 if( fabs((*vertex_Tracks)->pt()-muonTrackPt[l])<0.001 && 
 		     fabs((*vertex_Tracks)->eta()-muonTrackEta[l])<0.001 && 
 		     fabs((*vertex_Tracks)->phi()-muonTrackPhi[l])<0.001){
 		   pass_muon_assoc = true;
		   /// can I do this? // muon_vtx_idx->push_back(??);
		   continue; // go to next TRACK
		   //// I do this just in case, I don't want to associate the same track to diff leptons

 		 }
		 
 		 if( fabs((*vertex_Tracks)->pt()-electronTrackPt[k])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-electronTrackEta[k])<0.001 && 
 		     fabs((*vertex_Tracks)->phi()-electronTrackPhi[k])<0.001){
 		   pass_electron_assoc = true;
		   continue; // go to next TRACK
 		 }
 	       }//Loop over tracks to see if the leptons are there
	       
 	       if(pass_muon_assoc && pass_electron_assoc){
		 count_Nvertices_match_emu++;
 		 //veto_event=true;		 
 		 nextra_tracks_emu = vertexNtracks-2;
 		 (*vertex_extra_ntracks_emu).push_back(nextra_tracks_emu);
 		 (*vertex_emu_candMu_idx).push_back(l);
 		 (*vertex_emu_candE_idx).push_back(k);
 		 //		 cout<<"This event passes mu+ mu- criteria and should be vetoed"<<endl;
 	       }
	       
 	     }//end of if statement for charge
 	   }//end loop electron
	 }//end loop muon
       }//end of vertex Ntracks<75 requirementa
     }//end vertex z
   }// end of loop over vertices

   if(count_Nvertices_match_ee>1){
     cout<<"More than one matching ee pair vertex in the event"<<endl;
   }
   if(count_Nvertices_match_mumu>1){
     cout<<"More than one matching mumu pair vertex in the event"<<endl;
   }
   if(count_Nvertices_match_ee+count_Nvertices_match_mumu>1){
     cout<<"More than one matching same flavour ll pair vertex in the event"<<endl;
   }
   if(count_Nvertices_match_emu>1){
     cout<<"More than one matching emu pair vertex in the event"<<endl;
   }
   if(count_Nvertices_match_emu+count_Nvertices_match_ee+count_Nvertices_match_mumu>1){
     cout<<"More than one matching pair vertex of any kind in the event"<<endl;
   }

          
   newtree->Fill();


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
