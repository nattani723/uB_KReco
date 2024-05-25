#include "CCKaonAnalyzerRebuild_module.h"
#include "TrackRebuilder/ReconstructionOrchestrator.cc"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif

//const int kMaxTracks=20;

using namespace std;
namespace kaon_reconstruction{

  // ======================== Local Function Definition to get the reco origin ======================                                                                                                           
  art::Ptr<simb::MCTruth>TrackIDToMCTruth(art::Event const & evt, std::string _geant_producer, int geant_track_id)
  {
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(evt, _geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth) {
      if (iter.first->TrackId() == geant_track_id) {
	return iter.second;
      }
    }

    art::Ptr<simb::MCTruth> null_ptr;
    return null_ptr;
  }


//========================================================================
CCKaonAnalyzerRebuild::CCKaonAnalyzerRebuild(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
  fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora")),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","pandora")),
  //fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoraKalmanShowercali")),
  fPIDLabel                 (pset.get< std::string >("PIDLabel","pandorapid")),
  fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
  fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")), 
  fHitShowerAssns            (pset.get< std::string >("HitShowerAssn","pandora")), 
  m_pfp_producer            (pset.get< std::string >("pfp_producer","pandora")),
  fPFParticleLabel          (pset.get< std::string >("PFParticleLabel", "pandora")),
  fSpacePointproducer       (pset.get< std::string >("SpacePointproducer", "pandora")),
  fPandoraLabel             (pset.get< std::string >("PandoraLabel", "pandora")),
  //fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");
//  m_pandoraLabel            (pset.get< std::string >("PandoraLabel")),
//  m_is_verbose              (pset.get<bool>("Verbose",false)),
  isMC                      (pset.get< bool >("IsMC",true))

  //reco_track_dEdx(nullptr),
  //reco_track_ResRan(nullptr)
{



  // set dedx pdf parameters
  llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
  llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
  llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

  llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
  llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
  llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

  llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
  llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
  llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);


  llr_pid_calculator_k.set_dedx_binning(0, protonmuon_parameters_k.dedx_edges_pl_0);
  llr_pid_calculator_k.set_par_binning(0, protonmuon_parameters_k.parameters_edges_pl_0);
  llr_pid_calculator_k.set_lookup_tables(0, protonmuon_parameters_k.dedx_pdf_pl_0);

  llr_pid_calculator_k.set_dedx_binning(1, protonmuon_parameters_k.dedx_edges_pl_1);
  llr_pid_calculator_k.set_par_binning(1, protonmuon_parameters_k.parameters_edges_pl_1);
  llr_pid_calculator_k.set_lookup_tables(1, protonmuon_parameters_k.dedx_pdf_pl_1);

  llr_pid_calculator_k.set_dedx_binning(2, protonmuon_parameters_k.dedx_edges_pl_2);
  llr_pid_calculator_k.set_par_binning(2, protonmuon_parameters_k.parameters_edges_pl_2);
  llr_pid_calculator_k.set_lookup_tables(2, protonmuon_parameters_k.dedx_pdf_pl_2);
}
 
//========================================================================
CCKaonAnalyzerRebuild::~CCKaonAnalyzerRebuild()
{
  //destructor
}

//========================================================================
void CCKaonAnalyzerRebuild::endSubRun(const art::SubRun &subrun)
{
  //if (!m_isData)
  //{
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    m_pot = subrun.getByLabel("generator", potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    // -- std::cout << "[CCKaonAnalyzerRebuild::endSubRun] Storing POT info!" << std::endl;
  //}

  m_run = subrun.run();
  m_subrun = subrun.subRun();
  fSubrunTree->Fill();
}

//========================================================================
void CCKaonAnalyzerRebuild::beginJob()
{

  //initialiseCanvas();

  art::ServiceHandle<art::TFileService> tfs;

  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");

  //gInterpreter->GenerateDictionary("vector<vector<vector<Float_t>> >", "vector");
  fEventTree->Branch("event", &event, "event/I");
  fEventTree->Branch("run", &run, "run/I");
  fEventTree->Branch("subrun", &subrun, "surbrun/I");

  fEventTree->Branch("IsKaon", &IsKaon, "IsKaon/I");
  fEventTree->Branch("IsSingleKaon", &IsSingleKaon, "IsSingleKaon/I");
  fEventTree->Branch("IsAssociatedKaon", &IsAssociatedKaon, "IsAssociatedKaon/I");
  fEventTree->Branch("IsMuBR", &IsMuBR, "IsMuBR/I");
  fEventTree->Branch("IsPiBR", &IsPiBR, "IsPiBR/I");
  fEventTree->Branch("prip", &prip, "prip/I");
  fEventTree->Branch("prip_k_dau", &prip_k_dau, "prip_k_dau/I");


  fEventTree->Branch("process",&Process);
  fEventTree->Branch("endprocess",&EndProcess);

  fEventTree->Branch("event_weight" ,&event_weight);

  fEventTree->Branch("evtwgt_funcname" ,&evtwgt_funcname);
  fEventTree->Branch("evtwgt_weight" ,&evtwgt_weight);
  fEventTree->Branch("evtwgt_nweight" ,&evtwgt_nweight);

  fEventTree->Branch("true_nu_energy", &true_nu_energy, "true_nu_energy/F");
  fEventTree->Branch("true_nu_pdg", &true_nu_pdg, "true_nu_pdg/I");
  fEventTree->Branch("true_nu_mode", &true_nu_mode, "true_nu_mode/I");
  fEventTree->Branch("true_nu_ccnc", &true_nu_ccnc, "true_nu_ccnc/I");
  fEventTree->Branch("true_nu_vtx_x", &true_nu_vtx_x, "true_nu_vtx_x/F");
  fEventTree->Branch("true_nu_vtx_y", &true_nu_vtx_y, "true_nu_vtx_y/F");
  fEventTree->Branch("true_nu_vtx_z", &true_nu_vtx_z, "true_nu_vtx_z/F");
  fEventTree->Branch("true_nu_vtx_inTPC", &true_nu_vtx_inTPC, "true_nu_vtx_inTPC/O");
  fEventTree->Branch("true_nu_vtx_in5cmTPC", &true_nu_vtx_in5cmTPC, "true_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("true_nu_vtx_inCCInclusiveTPC", &true_nu_vtx_inCCInclusiveTPC, "true_nu_vtx_inCCInclusiveTPC/O");

  fEventTree->Branch("true_lepton_pdg", &true_lepton_pdg, "true_lepton_pdg/I");
  fEventTree->Branch("true_lepton_p", &true_lepton_p, "true_lepton_p/F");
  fEventTree->Branch("true_lepton_ke", &true_lepton_ke, "true_lepton_ke/F");
  fEventTree->Branch("true_lepton_theta", &true_lepton_theta, "true_lepton_theta/F");
  fEventTree->Branch("true_lepton_costheta", &true_lepton_costheta, "true_lepton_costheta/F");
  fEventTree->Branch("true_lepton_phi", &true_lepton_phi, "true_lepton_phi/F");

  fEventTree->Branch("true_nkaons", &true_nkaons, "true_nkaons/I");

  fEventTree->Branch("true_kaon_length", &true_kaon_length, "true_kaon_length/F");
  fEventTree->Branch("true_kaon_p", &true_kaon_p, "true_kaon_p/F");
  fEventTree->Branch("true_kaon_ke", &true_kaon_ke, "true_kaon_ke/F");
  fEventTree->Branch("true_kaon_theta", &true_kaon_theta, "true_kaon_theta/F");
  fEventTree->Branch("true_kaon_costheta", &true_kaon_costheta, "true_kaon_costheta/F");
  fEventTree->Branch("true_kaon_phi", &true_kaon_phi, "true_kaon_phi/F");
  fEventTree->Branch("true_kaon_ccmuon_angle", &true_kaon_ccmuon_angle, "true_kaon_ccmuon_angle/F");
  fEventTree->Branch("true_kaon_ccmuon_cosangle", &true_kaon_ccmuon_cosangle, "true_kaon_ccmuon_cosangle/F");
  fEventTree->Branch("true_kaon_end_process", &true_kaon_end_process, "true_kaon_end_process/I");
  fEventTree->Branch("true_kaon_end_ke", &true_kaon_end_ke, "true_kaon_end_ke/F");
  fEventTree->Branch("true_kaon_start_x", &true_kaon_start_x, "true_kaon_start_x/F");
  fEventTree->Branch("true_kaon_start_y", &true_kaon_start_y, "true_kaon_start_y/F");
  fEventTree->Branch("true_kaon_start_z", &true_kaon_start_z, "true_kaon_start_z/F");
  fEventTree->Branch("true_kaon_end_x", &true_kaon_end_x, "true_kaon_end_x/F");
  fEventTree->Branch("true_kaon_end_y", &true_kaon_end_y, "true_kaon_end_y/F");
  fEventTree->Branch("true_kaon_end_z", &true_kaon_end_z, "true_kaon_end_z/F");
  fEventTree->Branch("true_kaon_end_inTPC", &true_kaon_end_inTPC, "true_kaon_end_inTPC/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC", &true_kaon_end_in5cmTPC, "true_kaon_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC", &true_kaon_end_inCCInclusiveTPC, "true_kaon_end_inCCInclusiveTPC/O");


  fEventTree->Branch("true_dau_muon_length", &true_dau_muon_length, "true_dau_muon_length/F");
  fEventTree->Branch("true_dau_muon_p", &true_dau_muon_p, "true_dau_muon_p/F");
  fEventTree->Branch("true_dau_muon_ke", &true_dau_muon_ke, "true_dau_muon_ke/F");
  fEventTree->Branch("true_dau_muon_theta", &true_dau_muon_theta, "true_dau_muon_theta/F");
  fEventTree->Branch("true_dau_muon_costheta", &true_dau_muon_costheta, "true_dau_muon_costheta/F");
  fEventTree->Branch("true_dau_muon_phi", &true_dau_muon_phi, "true_dau_muon_phi/F");
  fEventTree->Branch("true_dau_muon_ccmuon_angle", &true_dau_muon_ccmuon_angle, "true_dau_muon_ccmuon_angle/F");
  fEventTree->Branch("true_dau_muon_ccmuon_cosangle", &true_dau_muon_ccmuon_cosangle, "true_dau_muon_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_muon_end_process", &true_dau_muon_end_process, "true_dau_muon_end_process/I");
  fEventTree->Branch("true_dau_muon_end_ke", &true_dau_muon_end_ke, "true_dau_muon_end_ke/F");
  fEventTree->Branch("true_dau_muon_start_x", &true_dau_muon_start_x, "true_dau_muon_start_x/F");
  fEventTree->Branch("true_dau_muon_start_y", &true_dau_muon_start_y, "true_dau_muon_start_y/F");
  fEventTree->Branch("true_dau_muon_start_z", &true_dau_muon_start_z, "true_dau_muon_start_z/F");
  fEventTree->Branch("true_dau_muon_end_x", &true_dau_muon_end_x, "true_dau_muon_end_x/F");
  fEventTree->Branch("true_dau_muon_end_y", &true_dau_muon_end_y, "true_dau_muon_end_y/F");
  fEventTree->Branch("true_dau_muon_end_z", &true_dau_muon_end_z, "true_dau_muon_end_z/F");

  fEventTree->Branch("true_dau_pip_length", &true_dau_pip_length, "true_dau_pip_length/F");
  fEventTree->Branch("true_dau_pip_p", &true_dau_pip_p, "true_dau_pip_p/F");
  fEventTree->Branch("true_dau_pip_ke", &true_dau_pip_ke, "true_dau_pip_ke/F");
  fEventTree->Branch("true_dau_pip_theta", &true_dau_pip_theta, "true_dau_pip_theta/F");
  fEventTree->Branch("true_dau_pip_costheta", &true_dau_pip_costheta, "true_dau_pip_costheta/F");
  fEventTree->Branch("true_dau_pip_phi", &true_dau_pip_phi, "true_dau_pip_phi/F");
  fEventTree->Branch("true_dau_pip_ccmuon_angle", &true_dau_pip_ccmuon_angle, "true_dau_pip_ccmuon_angle/F");
  fEventTree->Branch("true_dau_pip_ccmuon_cosangle", &true_dau_pip_ccmuon_cosangle, "true_dau_pip_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_pip_end_process", &true_dau_pip_end_process, "true_dau_pip_end_process/I");
  fEventTree->Branch("true_dau_pip_end_ke", &true_dau_pip_end_ke, "true_dau_pip_end_ke/F");
  fEventTree->Branch("true_dau_pip_start_x", &true_dau_pip_start_x, "true_dau_pip_start_x/F");
  fEventTree->Branch("true_dau_pip_start_y", &true_dau_pip_start_y, "true_dau_pip_start_y/F");
  fEventTree->Branch("true_dau_pip_start_z", &true_dau_pip_start_z, "true_dau_pip_start_z/F");
  fEventTree->Branch("true_dau_pip_end_x", &true_dau_pip_end_x, "true_dau_pip_end_x/F");
  fEventTree->Branch("true_dau_pip_end_y", &true_dau_pip_end_y, "true_dau_pip_end_y/F");
  fEventTree->Branch("true_dau_pip_end_z", &true_dau_pip_end_z, "true_dau_pip_end_z/F");


  fEventTree->Branch("true_dau_pin_length", &true_dau_pin_length, "true_dau_pin_length/F");
  fEventTree->Branch("true_dau_pin_p", &true_dau_pin_p, "true_dau_pin_p/F");
  fEventTree->Branch("true_dau_pin_ke", &true_dau_pin_ke, "true_dau_pin_ke/F");
  fEventTree->Branch("true_dau_pin_theta", &true_dau_pin_theta, "true_dau_pin_theta/F");
  fEventTree->Branch("true_dau_pin_costheta", &true_dau_pin_costheta, "true_dau_pin_costheta/F");
  fEventTree->Branch("true_dau_pin_phi", &true_dau_pin_phi, "true_dau_pin_phi/F");
  fEventTree->Branch("true_dau_pin_ccmuon_angle", &true_dau_pin_ccmuon_angle, "true_dau_pin_ccmuon_angle/F");
  fEventTree->Branch("true_dau_pin_ccmuon_cosangle", &true_dau_pin_ccmuon_cosangle, "true_dau_pin_ccmuon_cosangle/F");
  fEventTree->Branch("true_dau_pin_end_process", &true_dau_pin_end_process, "true_dau_pin_end_process/I");
  fEventTree->Branch("true_dau_pin_end_ke", &true_dau_pin_end_ke, "true_dau_pin_end_ke/F");
  fEventTree->Branch("true_dau_pin_start_x", &true_dau_pin_start_x, "true_dau_pin_start_x/F");
  fEventTree->Branch("true_dau_pin_start_y", &true_dau_pin_start_y, "true_dau_pin_start_y/F");
  fEventTree->Branch("true_dau_pin_start_z", &true_dau_pin_start_z, "true_dau_pin_start_z/F");
  fEventTree->Branch("true_dau_pin_end_x", &true_dau_pin_end_x, "true_dau_pin_end_x/F");
  fEventTree->Branch("true_dau_pin_end_y", &true_dau_pin_end_y, "true_dau_pin_end_y/F");
  fEventTree->Branch("true_dau_pin_end_z", &true_dau_pin_end_z, "true_dau_pin_end_z/F");

  fEventTree->Branch("cheat_num_hits", &cheat_num_hits, "cheat_num_hits[20][10]/I");
  fEventTree->Branch("cheat_ini_hits", &cheat_ini_hits, "cheat_ini_hits[20][10]/I");
  fEventTree->Branch("cheat_closest_distance", &cheat_closest_distance, "cheat_closest_distance[20][10]/F");

  fEventTree->Branch("cheat_peak_pdg", &cheat_peak_pdg, "cheat_peak_pdg[20][10]/F");
  fEventTree->Branch("cheat_peak_theta", &cheat_peak_theta, "cheat_peak_theta[20][10]/F");
  fEventTree->Branch("cheat_pip_trkln", &cheat_pip_trkln, "cheat_pip_trkln/F");
  fEventTree->Branch("cheat_mup_trkln", &cheat_mup_trkln, "cheat_mup_trkln/F");
  fEventTree->Branch("best_peak_trkln", &best_peak_trkln, "best_peak_trkln[20][10]/F");
  fEventTree->Branch("cheat_peak_phi", &cheat_peak_phi, "cheat_peak_phi[20][10]/F");

  fEventTree->Branch("n_recoRebDauTracks", &n_recoRebDauTracks, "n_recoRebDauTracks[20]/I"); 
  fEventTree->Branch("rebdautrack_length", &rebdautrack_length, "rebdautrack_length[20][10]/F"); 
  fEventTree->Branch("rebdautracktrue_length", &rebdautracktrue_length, "rebdautracktrue_length[20][10]/F"); 
  fEventTree->Branch("rebdautracktruedir_length", &rebdautracktruedir_length, "rebdautracktruedir_length[20][10]/F"); 
  fEventTree->Branch("rebdautrack_pdg", &rebdautrack_pdg, "rebdautrack_pdg[20][10]/F"); 
  fEventTree->Branch("best_peak_x", &best_peak_x, "best_peak_x[20][10]/F");
  fEventTree->Branch("best_peak_y", &best_peak_y, "best_peak_y[20][10]/F");
  fEventTree->Branch("best_peak_z", &best_peak_z, "best_peak_z[20][10]/F");
  fEventTree->Branch("best_peak_x_true", &best_peak_x_true, "best_peak_x_true[20][10]/F");
  fEventTree->Branch("best_peak_y_true", &best_peak_y_true, "best_peak_y_true[20][10]/F");
  fEventTree->Branch("best_peak_z_true", &best_peak_z_true, "best_peak_z_true[20][10]/F");


  fEventTree->Branch("true_length", &true_length, "true_length[10]/F");
  fEventTree->Branch("true_p", &true_p, "true_p[10]/F");
  fEventTree->Branch("true_ke", &true_ke, "true_ke[10]/F");
  fEventTree->Branch("true_theta", &true_theta, "true_theta[10]/F");
  fEventTree->Branch("true_costheta", &true_costheta, "true_costheta[10]/F");
  fEventTree->Branch("true_phi", &true_phi, "true_phi[10]/F");
  fEventTree->Branch("true_ccmuon_angle", &true_ccmuon_angle, "true_ccmuon_angle[10]/F");
  fEventTree->Branch("true_ccmuon_cosangle", &true_ccmuon_cosangle, "true_ccmuon_cosangle[10]/F");
  fEventTree->Branch("true_end_process", &true_end_process, "true_end_process[10]/I");
  fEventTree->Branch("true_end_ke", &true_end_ke, "true_end_ke[10]/F");
  fEventTree->Branch("true_end_x", &true_end_x, "true_end_x[10]/F");
  fEventTree->Branch("true_end_y", &true_end_y, "true_end_y[10]/F");
  fEventTree->Branch("true_end_z", &true_end_z, "true_end_z[10]/F");
  fEventTree->Branch("true_end_inTPC", &true_end_inTPC, "true_end_inTPC[10]/O");
  fEventTree->Branch("true_end_in5cmTPC", &true_end_in5cmTPC, "true_end_in5cmTPC[10]/O");
  fEventTree->Branch("true_end_inCCInclusiveTPC", &true_end_inCCInclusiveTPC, "true_end_inCCInclusiveTPC[10]/O");

  fEventTree->Branch("true_kaon_ndaughters", &true_kaon_ndaughters, "true_kaon_ndaughters/I");
  fEventTree->Branch("true_kaon_ndaughters_decay", &true_kaon_ndaughters_decay, "true_kaon_ndaughters_decay/I");
  fEventTree->Branch("true_kaon_ndaughters_inelastic", &true_kaon_ndaughters_inelastic, "true_kaon_ndaughters_inelastic/I");
  fEventTree->Branch("true_kaon_ndecmup", &true_kaon_ndecmup, "true_kaon_ndecmup/I");
  fEventTree->Branch("true_kaon_ndecpip", &true_kaon_ndecpip, "true_kaon_ndecpip/I");
  fEventTree->Branch("true_kaon_ninekap", &true_kaon_ninekap, "true_kaon_ninekap/I");
  fEventTree->Branch("true_kaon_ninepip", &true_kaon_ninepip, "true_kaon_ninepip/I");
  fEventTree->Branch("true_kaon_ninepro", &true_kaon_ninepro, "true_kaon_ninepro/I");
  fEventTree->Branch("true_kaon_daughter_length", &true_kaon_daughter_length, "true_kaon_daughter_length/F");
  fEventTree->Branch("true_kaon_daughter_p", &true_kaon_daughter_p, "true_kaon_daughter_p/F");
  fEventTree->Branch("true_kaon_daughter_ke", &true_kaon_daughter_ke, "true_kaon_daughter_ke/F");
  fEventTree->Branch("true_kaon_daughter_theta", &true_kaon_daughter_theta, "true_kaon_daughter_theta/F");
  fEventTree->Branch("true_kaon_daughter_costheta", &true_kaon_daughter_costheta, "true_kaon_daughter_costheta/F");
  fEventTree->Branch("true_kaon_daughter_angle", &true_kaon_daughter_angle, "true_kaon_daughter_angle/F");
  fEventTree->Branch("true_kaon_daughter_cosangle", &true_kaon_daughter_cosangle, "true_kaon_daughter_cosangle/F");
  fEventTree->Branch("true_kaon_daughter_pdg", &true_kaon_daughter_pdg, "true_kaon_daughter_pdg/I");
  fEventTree->Branch("true_kaon_daughter_start_x", &true_kaon_daughter_start_x, "true_kaon_daughter_start_x/F");
  fEventTree->Branch("true_kaon_daughter_start_y", &true_kaon_daughter_start_y, "true_kaon_daughter_start_y/F");
  fEventTree->Branch("true_kaon_daughter_start_z", &true_kaon_daughter_start_z, "true_kaon_daughter_start_z/F");
  fEventTree->Branch("true_kaon_daughter_end_x", &true_kaon_daughter_end_x, "true_kaon_daughter_end_x/F");
  fEventTree->Branch("true_kaon_daughter_end_y", &true_kaon_daughter_end_y, "true_kaon_daughter_end_y/F");
  fEventTree->Branch("true_kaon_daughter_end_z", &true_kaon_daughter_end_z, "true_kaon_daughter_end_z/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC", &true_kaon_daughter_end_inTPC, "true_kaon_daughter_endInTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC", &true_kaon_daughter_end_in5cmTPC, "true_kaon_daughter_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC", &true_kaon_daughter_end_inCCInclusiveTPC, "true_kaon_daughter_end_inCCInclusiveTPC/O");

  fEventTree->Branch("true_nhyperons", &true_nhyperons, "true_nhyperons/I");



  fEventTree->Branch("reco_nu_cc_filter", &reco_nu_cc_filter, "reco_nu_cc_filter/O");

  fEventTree->Branch("reco_nu_vtx_x", &reco_nu_vtx_x, "reco_nu_vtx_x/F");
  fEventTree->Branch("reco_nu_vtx_y", &reco_nu_vtx_y, "reco_nu_vtx_y/F");
  fEventTree->Branch("reco_nu_vtx_z", &reco_nu_vtx_z, "reco_nu_vtx_z/F");
  fEventTree->Branch("reco_nu_vtx_inTPC", &reco_nu_vtx_inTPC, "reco_nu_vtx_inTPC/O");
  fEventTree->Branch("reco_nu_vtx_in5cmTPC", &reco_nu_vtx_in5cmTPC, "reco_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_nu_vtx_inCCInclusiveTPC", &reco_nu_vtx_inCCInclusiveTPC, "reco_nu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_nu_ndaughters", &reco_nu_ndaughters, "reco_nu_ndaughters/I"          );
  fEventTree->Branch("reco_nu_cc_nmue", &reco_nu_cc_nmue, "reco_nu_cc_nmue/I"          );

  fEventTree->Branch("reco_ccmu_vtx_x", &reco_ccmu_vtx_x, "reco_ccmu_vtx_x/F");
  fEventTree->Branch("reco_ccmu_vtx_y", &reco_ccmu_vtx_y, "reco_ccmu_vtx_y/F");
  fEventTree->Branch("reco_ccmu_vtx_z", &reco_ccmu_vtx_z, "reco_ccmu_vtx_z/F");
  fEventTree->Branch("reco_ccmu_vtx_inTPC",&reco_ccmu_vtx_inTPC, "reco_ccmu_vtx_inTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_in5cmTPC", &reco_ccmu_vtx_in5cmTPC, "reco_ccmu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_inCCInclusiveTPC", &reco_ccmu_vtx_inCCInclusiveTPC, "reco_ccmu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_pdg", &reco_ccmu_true_pdg, "reco_ccmu_true_pdg/I");
  fEventTree->Branch("reco_ccmu_true_origin", &reco_ccmu_true_origin, "reco_ccmu_true_origin/I");
  fEventTree->Branch("reco_ccmu_true_primary", &reco_ccmu_true_primary, "reco_ccmu_true_primary/O");
  fEventTree->Branch("reco_ccmu_true_end_inTPC", &reco_ccmu_true_end_inTPC, "reco_ccmu_true_end_inTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_in5cmTPC", &reco_ccmu_true_end_in5cmTPC, "reco_ccmu_true_end_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_inCCInclusiveTPC", &reco_ccmu_true_end_inCCInclusiveTPC, "reco_ccmu_true_end_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_length", &reco_ccmu_true_length, "reco_ccmu_true_length/F");

  fEventTree->Branch("reco_ntracks", &reco_ntracks, "reco_ntracks/I");
  fEventTree->Branch("reco_nshowers", &reco_ntracks, "reco_nshowers/I");
  fEventTree->Branch("nTracks", &nTracks, "nTracks/I");
  fEventTree->Branch("nShowers", &nShowers, "nShowers/I");
  //  fEventTree->Branch("m_reco_track_sliceId", &m_reco_track_sliceId, "m_reco_track_sliceId[20]/F");
  //fEventTree->Branch("m_reco_track_is_nuslice", &m_reco_track_is_nuslice, "m_reco_track_is_nuslice[20]/F");

  fEventTree->Branch("reco_track_start_x", &reco_track_start_x, "reco_track_start_x[20]/F");
  fEventTree->Branch("reco_track_start_y", &reco_track_start_y, "reco_track_start_y[20]/F");
  fEventTree->Branch("reco_track_start_z", &reco_track_start_z, "reco_track_start_z[20]/F");
  fEventTree->Branch("reco_track_end_x", &reco_track_end_x, "reco_track_end_x[20]/F");
  fEventTree->Branch("reco_track_end_y", &reco_track_end_y, "reco_track_end_y[20]/F");
  fEventTree->Branch("reco_track_end_z", &reco_track_end_z, "reco_track_end_z[20]/F");


  fEventTree->Branch("reco_track_distance", &reco_track_distance, "reco_track_distance[20]/F");
  fEventTree->Branch("reco_track_nhits0", &reco_track_nhits0, "reco_track_nhits0[20]/I");
  fEventTree->Branch("reco_track_nhits1", &reco_track_nhits1, "reco_track_nhits1[20]/I");
  fEventTree->Branch("reco_track_nhits2", &reco_track_nhits2, "reco_track_nhits2[20]/I");

  fEventTree->Branch("reco_track_kin0", &reco_track_kin0, "reco_track_kin0[20]/F");
  fEventTree->Branch("reco_track_kin1", &reco_track_kin1, "reco_track_kin1[20]/F");
  fEventTree->Branch("reco_track_kin2", &reco_track_kin2, "reco_track_kin2[20]/F");
  //fEventTree->Branch("reco_track_dEdx" ,&reco_track_dEdx);
  //fEventTree->Branch("reco_track_ResRan" ,&reco_track_ResRan);
  //

  fEventTree->Branch("reco_track_dEdx_pl0", &reco_track_dEdx_pl0,"reco_track_dEdx_pl0[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl0", &reco_track_ResRan_pl0,"reco_track_ResRan_pl0[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl1", &reco_track_dEdx_pl1,"reco_track_dEdx_pl1[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl1", &reco_track_ResRan_pl1,"reco_track_ResRan_pl1[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl2", &reco_track_dEdx_pl2,"reco_track_dEdx_pl2[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl2", &reco_track_ResRan_pl2,"reco_track_ResRan_pl2[20][2000]");

  fEventTree->Branch("reco_track_length", &reco_track_length, "reco_track_length[20]/F");
  fEventTree->Branch("reco_track_theta", &reco_track_theta, "reco_track_theta[20]/F");
  fEventTree->Branch("reco_track_phi", &reco_track_phi, "reco_track_phi[20]/F");
  fEventTree->Branch("reco_track_dir", &reco_track_dir, "reco_track_dir[20]/O");

  fEventTree->Branch("reco_track_P_vtx", &reco_track_P_vtx, "reco_track_P_vtx[20]/F");
  fEventTree->Branch("reco_track_P_str", &reco_track_P_str, "reco_track_P_str[20]/F");
  fEventTree->Branch("reco_track_P_end", &reco_track_P_end, "reco_track_P_end[20]/F");

  fEventTree->Branch("reco_track_chi2ka_pl0", &reco_track_chi2ka_pl0, "reco_track_chi2ka_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl0", &reco_track_chi2pr_pl0, "reco_track_chi2pr_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl0", &reco_track_chi2pi_pl0, "reco_track_chi2pi_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl0", &reco_track_chi2mu_pl0, "reco_track_chi2mu_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl1", &reco_track_chi2ka_pl1, "reco_track_chi2ka_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl1", &reco_track_chi2pr_pl1, "reco_track_chi2pr_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl1", &reco_track_chi2pi_pl1, "reco_track_chi2pi_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl1", &reco_track_chi2mu_pl1, "reco_track_chi2mu_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl2", &reco_track_chi2ka_pl2, "reco_track_chi2ka_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl2", &reco_track_chi2pr_pl2, "reco_track_chi2pr_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl2", &reco_track_chi2pi_pl2, "reco_track_chi2pi_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl2", &reco_track_chi2mu_pl2, "reco_track_chi2mu_pl2[20]/F");

  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl0", &reco_track_Bragg_fwd_ka_pl0, "reco_track_Bragg_fwd_ka_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl0", &reco_track_Bragg_fwd_pr_pl0, "reco_track_Bragg_fwd_pr_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl0", &reco_track_Bragg_fwd_pi_pl0, "reco_track_Bragg_fwd_pi_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl0", &reco_track_Bragg_fwd_mu_pl0, "reco_track_Bragg_fwd_mu_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl1", &reco_track_Bragg_fwd_ka_pl1, "reco_track_Bragg_fwd_ka_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl1", &reco_track_Bragg_fwd_pr_pl1, "reco_track_Bragg_fwd_pr_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl1", &reco_track_Bragg_fwd_pi_pl1, "reco_track_Bragg_fwd_pi_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl1", &reco_track_Bragg_fwd_mu_pl1, "reco_track_Bragg_fwd_mu_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl2", &reco_track_Bragg_fwd_ka_pl2, "reco_track_Bragg_fwd_ka_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl2", &reco_track_Bragg_fwd_pr_pl2, "reco_track_Bragg_fwd_pr_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl2", &reco_track_Bragg_fwd_pi_pl2, "reco_track_Bragg_fwd_pi_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl2", &reco_track_Bragg_fwd_mu_pl2, "reco_track_Bragg_fwd_mu_pl2[20]/F");


  fEventTree->Branch("reco_track_MIP_pl0", &reco_track_MIP_pl0, "reco_track_MIP_pl0[20]/F");
  fEventTree->Branch("reco_track_MIP_pl1", &reco_track_MIP_pl1, "reco_track_MIP_pl1[20]/F");
  fEventTree->Branch("reco_track_MIP_pl2", &reco_track_MIP_pl2, "reco_track_MIP_pl2[20]/F");


  fEventTree->Branch("reco_track_chi2ka_3pl", &reco_track_chi2ka_3pl, "reco_track_chi2ka_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pr_3pl", &reco_track_chi2pr_3pl, "reco_track_chi2pr_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pi_3pl", &reco_track_chi2pi_3pl, "reco_track_chi2pi_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2mu_3pl", &reco_track_chi2mu_3pl, "reco_track_chi2mu_3pl[20]/F");
  fEventTree->Branch("reco_track_likepr_3pl", &reco_track_likepr_3pl, "reco_track_likepr_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_3pl", &reco_track_llrpid_3pl, "reco_track_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_total_llrpid_3pl", &reco_track_total_llrpid_3pl, "reco_track_total_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_k_3pl", &reco_track_llrpid_k_3pl, "reco_track_llrpid_k_3pl[20]/F");
  fEventTree->Branch("reco_track_vtx_inTPC", &reco_track_vtx_inTPC, "reco_track_vtx_inTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_in5cmTPC", &reco_track_vtx_in5cmTPC, "reco_track_vtx_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_inCCInclusiveTPC", &reco_track_vtx_inCCInclusiveTPC, "reco_track_vtx_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_end_inTPC", &reco_track_end_inTPC, "reco_track_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_end_in5cmTPC", &reco_track_end_in5cmTPC, "reco_track_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_end_inCCInclusiveTPC", &reco_track_end_inCCInclusiveTPC, "reco_track_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_pdg", &reco_track_true_pdg, "reco_track_true_pdg[20]/I");
  fEventTree->Branch("reco_track_true_origin", &reco_track_true_origin, "reco_track_true_origin[20]/I");
  fEventTree->Branch("reco_track_true_primary", &reco_track_true_primary, "reco_track_true_primary[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC", &reco_track_true_end_inTPC, "reco_track_true_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC", &reco_track_true_end_in5cmTPC, "reco_track_true_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC", &reco_track_true_end_inCCInclusiveTPC, "reco_track_true_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_length", &reco_track_true_length, "reco_track_true_length[20]/F");


  fEventTree->Branch("reco_track_match_e", &reco_track_match_e, "reco_track_match_e[20][10]/F");
  fEventTree->Branch("reco_track_match_hit", &reco_track_match_hit, "reco_track_match_hit[20][10]/F");
  fEventTree->Branch("reco_track_match_epdg", &reco_track_match_epdg, "reco_track_match_epdg[20][10]/I");
  fEventTree->Branch("reco_track_match_hitpdg", &reco_track_match_hitpdg, "reco_track_match_hitpdg[20][10]/I");

  fEventTree->Branch("reco_track_daughter_match_e", &reco_track_daughter_match_e, "reco_track_daughter_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_hit", &reco_track_daughter_match_hit, "reco_track_daughter_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_epdg", &reco_track_daughter_match_epdg, "reco_track_daughter_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_match_hitpdg", &reco_track_daughter_match_hitpdg, "reco_track_daughter_match_hitpdg[20][20][10]/I");

  fEventTree->Branch("reco_track_daughter_shower_match_e", &reco_track_daughter_shower_match_e, "reco_track_daughter_shower_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_hit", &reco_track_daughter_shower_match_hit, "reco_track_daughter_shower_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_epdg", &reco_track_daughter_shower_match_epdg, "reco_track_daughter_shower_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_shower_match_hitpdg", &reco_track_daughter_shower_match_hitpdg, "reco_track_daughter_shower_match_hitpdg[20][20][10]/I");


  fEventTree->Branch("reco_track_true_pdg_sh", &reco_track_true_pdg_sh, "reco_track_true_pdg_sh[20]/I");
  fEventTree->Branch("reco_track_true_origin_sh", &reco_track_true_origin_sh, "reco_track_true_origin_sh[20]/I");
  fEventTree->Branch("reco_track_true_primary_sh", &reco_track_true_primary_sh, "reco_track_true_primary_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC_sh", &reco_track_true_end_inTPC_sh, "reco_track_true_end_inTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC_sh", &reco_track_true_end_in5cmTPC_sh, "reco_track_true_end_in5cmTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC_sh", &reco_track_true_end_inCCInclusiveTPC_sh, "reco_track_true_end_inCCInclusiveTPC[20]_sh/O");
  fEventTree->Branch("reco_track_true_length_sh", &reco_track_true_length_sh, "reco_track_true_length_sh[20]/F");

  fEventTree->Branch("reco_track_daughter_true_pdg_sh", &reco_track_daughter_true_pdg_sh, "reco_track_daughter_true_pdg_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin_sh", &reco_track_daughter_true_origin_sh, "reco_track_daughter_true_origin_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary_sh", &reco_track_daughter_true_primary_sh, "reco_track_daughter_true_primary[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC_sh", &reco_track_daughter_true_end_inTPC_sh, "reco_track_daughter_true_end_inTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC_sh", &reco_track_daughter_true_end_in5cmTPC_sh, "reco_track_daughter_true_end_in5cmTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC_sh", &reco_track_daughter_true_end_inCCInclusiveTPC_sh, "reco_track_daughter_true_end_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length_sh", &reco_track_daughter_true_length_sh, "reco_track_daughter_true_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother_sh", &reco_track_daughter_true_mother_sh, "reco_track_daughter_true_mother_sh[20][20]/I");

  fEventTree->Branch("reco_track_daughter_start_x", &reco_track_daughter_start_x, "reco_track_daughter_start_x[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_y", &reco_track_daughter_start_y, "reco_track_daughter_start_y[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_z", &reco_track_daughter_start_z, "reco_track_daughter_start_z[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_end_x", &reco_track_daughter_end_x, "reco_track_daughter_end_x[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_y", &reco_track_daughter_end_y, "reco_track_daughter_end_y[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_z", &reco_track_daughter_end_z, "reco_track_daughter_end_z[20][20]/F");

  fEventTree->Branch("reco_shower_ndaughters", &reco_shower_ndaughters, "reco_shower_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance_sh", &reco_track_daughter_distance_sh, "reco_track_daughter_distance_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance_sh", &reco_track_daughter_vtx_distance_sh, "reco_track_daughter_vtx_distance_sh[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter_sh", &reco_angle_track_daughter_sh, "reco_angle_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_angle_daughter_track_daughter_sh", &reco_angle_daughter_track_daughter_sh, "reco_angle_daughter_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_length_sh", &reco_track_daughter_length_sh, "reco_track_daughter_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta_sh", &reco_track_daughter_theta_sh, "reco_track_daughter_theta_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi_sh", &reco_track_daughter_phi_sh, "reco_track_daughter_phi_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_open_angle_sh", &reco_track_daughter_open_angle_sh, "reco_track_daughter_open_angle_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl0_sh", &reco_track_daughter_dedx_pl0_sh, "reco_track_daughter_dedx_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl1_sh", &reco_track_daughter_dedx_pl1_sh, "reco_track_daughter_dedx_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl2_sh", &reco_track_daughter_dedx_pl2_sh, "reco_track_daughter_dedx_pl2_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl0_sh", &reco_track_daughter_energy_pl0_sh, "reco_track_daughter_energy_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl1_sh", &reco_track_daughter_energy_pl1_sh, "reco_track_daughter_energy_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl2_sh", &reco_track_daughter_energy_pl2_sh, "reco_track_daughter_energy_pl2_sh[20][20]/F");

  fEventTree->Branch("reco_track_daughter_vtx_inTPC_sh", &reco_track_daughter_vtx_inTPC_sh, "reco_track_daughter_vtx_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC_sh", &reco_track_daughter_vtx_in5cmTPC_sh, "reco_track_daughter_vtx_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC_sh", &reco_track_daughter_vtx_inCCInclusiveTPC_sh, "reco_track_daughter_vtx_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC_sh", &reco_track_daughter_end_inTPC_sh, "reco_track_daughter_end_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC_sh", &reco_track_daughter_end_in5cmTPC_sh, "reco_track_daughter_end_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC_sh", &reco_track_daughter_end_inCCInclusiveTPC_sh, "reco_track_daughter_end_inCCInclusiveTPC_sh[20][20]/O");

  fEventTree->Branch("reco_track_ndaughters", &reco_track_ndaughters, "reco_track_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance", &reco_track_daughter_distance, "reco_track_daughter_distance[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance", &reco_track_daughter_vtx_distance, "reco_track_daughter_vtx_distance[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter", &reco_angle_track_daughter, "reco_angle_track_daughter[20][20]/F");
  fEventTree->Branch("reco_track_daughter_nhits0", &reco_track_daughter_nhits0, "reco_track_daughter_nhits0[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits1", &reco_track_daughter_nhits1, "reco_track_daughter_nhits1[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits2", &reco_track_daughter_nhits2, "reco_track_daughter_nhits2[20][20]/I");

  //fEventTree->Branch("reco_track_daughter_dEdx", &reco_track_daughter_dEdx);
  //fEventTree->Branch("reco_track_daughter_ResRan", &reco_track_daughter_ResRan);
  //
  fEventTree->Branch("reco_track_daughter_dEdx_pl0", &reco_track_daughter_dEdx_pl0, "reco_track_daughter_dEdx_pl0[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl0", &reco_track_daughter_ResRan_pl0, "reco_track_daughter_ResRan_pl0[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl1", &reco_track_daughter_dEdx_pl1, "reco_track_daughter_dEdx_pl1[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl1", &reco_track_daughter_ResRan_pl1, "reco_track_daughter_ResRan_pl1[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl2", &reco_track_daughter_dEdx_pl2, "reco_track_daughter_dEdx_pl2[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl2", &reco_track_daughter_ResRan_pl2, "reco_track_daughter_ResRan_pl2[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_length", &reco_track_daughter_length, "reco_track_daughter_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta", &reco_track_daughter_theta, "reco_track_daughter_theta[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi", &reco_track_daughter_phi, "reco_track_daughter_phi[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl0", &reco_track_daughter_chi2ka_pl0, "reco_track_daughter_chi2ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl0", &reco_track_daughter_chi2pr_pl0, "reco_track_daughter_chi2pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl0", &reco_track_daughter_chi2pi_pl0, "reco_track_daughter_chi2pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl0", &reco_track_daughter_chi2mu_pl0, "reco_track_daughter_chi2mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl1", &reco_track_daughter_chi2ka_pl1, "reco_track_daughter_chi2ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl1", &reco_track_daughter_chi2pr_pl1, "reco_track_daughter_chi2pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl1", &reco_track_daughter_chi2pi_pl1, "reco_track_daughter_chi2pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl1", &reco_track_daughter_chi2mu_pl1, "reco_track_daughter_chi2mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl2", &reco_track_daughter_chi2ka_pl2, "reco_track_daughter_chi2ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl2", &reco_track_daughter_chi2pr_pl2, "reco_track_daughter_chi2pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl2", &reco_track_daughter_chi2pi_pl2, "reco_track_daughter_chi2pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl2", &reco_track_daughter_chi2mu_pl2, "reco_track_daughter_chi2mu_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_3pl", &reco_track_daughter_chi2ka_3pl, "reco_track_daughter_chi2ka_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_3pl", &reco_track_daughter_chi2pr_3pl, "reco_track_daughter_chi2pr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_3pl", &reco_track_daughter_chi2pi_3pl, "reco_track_daughter_chi2pi_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_3pl", &reco_track_daughter_chi2mu_3pl, "reco_track_daughter_chi2mu_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_likepr_3pl", &reco_track_daughter_likepr_3pl, "reco_track_daughter_likepr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_3pl", &reco_track_daughter_llrpid_3pl, "reco_track_daughter_llrpid_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_k_3pl", &reco_track_daughter_llrpid_k_3pl, "reco_track_daughter_llrpid_k_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC", &reco_track_daughter_vtx_inTPC, "reco_track_daughter_vtx_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC", &reco_track_daughter_vtx_in5cmTPC, "reco_track_daughter_vtx_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC", &reco_track_daughter_vtx_inCCInclusiveTPC, "reco_track_daughter_vtx_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC", &reco_track_daughter_end_inTPC, "reco_track_daughter_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC", &reco_track_daughter_end_in5cmTPC, "reco_track_daughter_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC", &reco_track_daughter_end_inCCInclusiveTPC, "reco_track_daughter_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_pdg", &reco_track_daughter_true_pdg, "reco_track_daughter_true_pdg[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin", &reco_track_daughter_true_origin, "reco_track_daughter_true_origin[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary", &reco_track_daughter_true_primary, "reco_track_daughter_true_primary[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC", &reco_track_daughter_true_end_inTPC, "reco_track_daughter_true_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC", &reco_track_daughter_true_end_in5cmTPC, "reco_track_daughter_true_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC", &reco_track_daughter_true_end_inCCInclusiveTPC, "reco_track_daughter_true_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length", &reco_track_daughter_true_length, "reco_track_daughter_true_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother", &reco_track_daughter_true_mother, "reco_track_daughter_true_mother[20][20]/I");


  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl0", &reco_track_daughter_Bragg_fwd_ka_pl0, "reco_track_daughter_Bragg_fwd_ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl0", &reco_track_daughter_Bragg_fwd_pr_pl0, "reco_track_daughter_Bragg_fwd_pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl0", &reco_track_daughter_Bragg_fwd_pi_pl0, "reco_track_daughter_Bragg_fwd_pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl0", &reco_track_daughter_Bragg_fwd_mu_pl0, "reco_track_daughter_Bragg_fwd_mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl1", &reco_track_daughter_Bragg_fwd_ka_pl1, "reco_track_daughter_Bragg_fwd_ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl1", &reco_track_daughter_Bragg_fwd_pr_pl1, "reco_track_daughter_Bragg_fwd_pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl1", &reco_track_daughter_Bragg_fwd_pi_pl1, "reco_track_daughter_Bragg_fwd_pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl1", &reco_track_daughter_Bragg_fwd_mu_pl1, "reco_track_daughter_Bragg_fwd_mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl2", &reco_track_daughter_Bragg_fwd_ka_pl2, "reco_track_daughter_Bragg_fwd_ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl2", &reco_track_daughter_Bragg_fwd_pr_pl2, "reco_track_daughter_Bragg_fwd_pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl2", &reco_track_daughter_Bragg_fwd_pi_pl2, "reco_track_daughter_Bragg_fwd_pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl2", &reco_track_daughter_Bragg_fwd_mu_pl2, "reco_track_daughter_Bragg_fwd_mu_pl2[20][20]/F");


  fEventTree->Branch("reco_track_daughter_MIP_pl0", &reco_track_daughter_MIP_pl0, "reco_track_daughter_MIP_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl1", &reco_track_daughter_MIP_pl1, "reco_track_daughter_MIP_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl2", &reco_track_daughter_MIP_pl2, "reco_track_daughter_MIP_pl2[20][20]/F");



  fEventTree->Branch("k_can_trkid", &k_can_trkid,"k_can_trkid/I");
  fEventTree->Branch("mu_can_trkid", &mu_can_trkid,"mu_can_trkid/I");
  fEventTree->Branch("k_mu_can_dis", &k_mu_can_dis,"k_mu_can_dis/F");
  fEventTree->Branch("k_mu_open_angle", &k_mu_open_angle,"k_mu_open_angle/F");
  fEventTree->Branch("k_vtx_dis", &k_vtx_dis,"k_vtx_dis/F");
  //fEventTree->Branch("k_geant_ID", &k_geant_ID,"k_geant_ID/I");
  //fEventTree->Branch("k_origin", &k_origin,"k_origin/I");
  //fEventTree->Branch("k_pdg", &k_pdg,"k_pdg/I");
  //fEventTree->Branch("k_isPri", &k_isPri,"k_isPri/I");
  //fEventTree->Branch("k_endE", &k_endE,"k_endE/F");
  //fEventTree->Branch("k_ness", &k_ness,"k_ness/I");
  //fEventTree->Branch("kaon_vtx_dis", &kaon_vtx_dis,"kaon_vtx_dis/F");
  //fEventTree->Branch("k_plen", &k_plen,"k_plen/F");
  //fEventTree->Branch("k_phi", &k_phi,"k_phi/F");
  //fEventTree->Branch("k_theta", &k_theta,"k_theta/F");
  //fEventTree->Branch("k_in_5_TPC", &k_in_5_TPC,"k_in_5_TPC/I");
  //fEventTree->Branch("k_in_CC_TPC", &k_in_CC_TPC,"k_in_CC_TPC/I");
  //fEventTree->Branch("k_hit", &k_hit,"k_hit/I");
  //fEventTree->Branch("k_range", &k_range,"k_range/F");
  //fEventTree->Branch("k_KE", &k_KE,"k_KE/F");
  //fEventTree->Branch("k_large_dedx", &k_large_dedx,"k_large_dedx/F");
  //fEventTree->Branch("k_small_dedx", &k_small_dedx,"k_small_dedx/F");
  //fEventTree->Branch("k_chi_p", &k_chi_p,"k_chi_p/F");
  //fEventTree->Branch("k_chi_k", &k_chi_k,"k_chi_k/F");
  //fEventTree->Branch("k_chi_pi", &k_chi_pi,"k_chi_pi/F");
  //fEventTree->Branch("k_chi_mu", &k_chi_mu,"k_chi_mu/F");
  //fEventTree->Branch("k_p_max", &k_p_max,"k_p_max/F");
  //fEventTree->Branch("k_k_max", &k_k_max,"k_k_max/F");
  //fEventTree->Branch("k_pi_max", &k_pi_max,"k_pi_max/F");
  //fEventTree->Branch("k_mu_max", &k_mu_max,"k_mu_max/F");
  //fEventTree->Branch("k_mip_max", &k_mip_max,"k_mip_max/F");
  //fEventTree->Branch("k_L1_ratio", &k_L1_ratio,"k_L1_ratio/F");
  //fEventTree->Branch("k_LL1", &k_LL1,"k_LL1/F");
  //fEventTree->Branch("k_L2_ratio", &k_L2_ratio,"k_L2_ratio/F");
  //fEventTree->Branch("k_LL2", &k_LL2,"k_LL2/F");
  //fEventTree->Branch("k_Lp_ratio", &k_Lp_ratio,"k_Lp_ratio/F");
  //fEventTree->Branch("k_LLp", &k_LLp,"k_LLp/F");
  //fEventTree->Branch("k_Lk_ratio", &k_Lk_ratio,"k_Lk_ratio/F");
  //fEventTree->Branch("k_LLk", &k_LLk,"k_LLk/F");
  //fEventTree->Branch("k_pida_mean", &k_pida_mean,"k_pida_mean/F");
  //fEventTree->Branch("k_pida_med", &k_pida_med,"k_pida_med/F");
  //fEventTree->Branch("k_kde", &k_kde,"k_kde/F");
  //fEventTree->Branch("k_trm_dqdx", &k_trm_dqdx,"k_trm_dqdx/F");
  //fEventTree->Branch("k_trm_dedx", &k_trm_dedx,"k_trm_dedx/F");
  //fEventTree->Branch("k_dedx",k_dedx,"k_dedx[3][3000]/F");
  //fEventTree->Branch("k_rr",k_rr,"k_rr[3][3000]/F");
  //fEventTree->Branch("mu_pdg", &mu_pdg,"mu_pdg/I");
  //fEventTree->Branch("mu_isDec", &mu_isDec,"mu_isDec/I");
  //fEventTree->Branch("mu_origin", &mu_origin,"mu_origin/I");
  //fEventTree->Branch("mu_k_is_Mother", &mu_k_is_Mother,"mu_k_is_Mother/I");
  //fEventTree->Branch("mu_mom_k_inelas", &mu_mom_k_inelas,"mu_mom_k_inelas/I");
  //fEventTree->Branch("mu_ness", &mu_ness,"mu_ness/I");
  //fEventTree->Branch("mu_plen", &mu_plen,"mu_plen/F");
  //fEventTree->Branch("mu_phi", &mu_phi,"mu_phi/F");
  //fEventTree->Branch("mu_theta", &mu_theta,"mu_theta/F");
  //fEventTree->Branch("mu_in_5_TPC", &mu_in_5_TPC,"mu_in_5_TPC/I");
  //fEventTree->Branch("mu_in_CC_TPC", &mu_in_CC_TPC,"mu_in_CC_TPC/I");
  //fEventTree->Branch("mu_KE", &mu_KE,"mu_KE/F");
  //fEventTree->Branch("mu_hit", &mu_hit,"mu_hit/I");
  //fEventTree->Branch("mu_range", &mu_range,"mu_range/F");
  //fEventTree->Branch("mu_large_dedx", &mu_large_dedx,"mu_large_dedx/F");
  //fEventTree->Branch("mu_small_dedx", &mu_small_dedx,"mu_small_dedx/F");
  //fEventTree->Branch("mu_chi_p", &mu_chi_p,"mu_chi_p/F");
  //fEventTree->Branch("mu_chi_k", &mu_chi_k,"mu_chi_k/F");
  //fEventTree->Branch("mu_chi_pi", &mu_chi_pi,"mu_chi_pi/F");
  //fEventTree->Branch("mu_chi_mu", &mu_chi_mu,"mu_chi_mu/F");
  //fEventTree->Branch("mu_p_max", &mu_p_max,"mu_p_max/F");
  //fEventTree->Branch("mu_k_max", &mu_k_max,"mu_k_max/F");
  //fEventTree->Branch("mu_pi_max", &mu_pi_max,"mu_pi_max/F");
  //fEventTree->Branch("mu_mu_max", &mu_mu_max,"mu_mu_max/F");
  //fEventTree->Branch("mu_mip_max", &mu_mip_max,"mu_mip_max/F");
  //fEventTree->Branch("mu_L1_ratio", &mu_L1_ratio,"mu_L1_ratio/F");
  //fEventTree->Branch("mu_LL1", &mu_LL1,"mu_LL1/F");
  //fEventTree->Branch("mu_L2_ratio", &mu_L2_ratio,"mu_L2_ratio/F");
  //fEventTree->Branch("mu_LL2", &mu_LL2,"mu_LL2/F");
  //fEventTree->Branch("mu_Lp_ratio", &mu_Lp_ratio,"mu_Lp_ratio/F");
  //fEventTree->Branch("mu_LLp", &mu_LLp,"mu_LLp/F");
  //fEventTree->Branch("mu_Lk_ratio", &mu_Lk_ratio,"mu_Lk_ratio/F");
  //fEventTree->Branch("mu_LLk", &mu_LLk,"mu_LLk/F");
  //fEventTree->Branch("mu_pida_mean", &mu_pida_mean,"mu_pida_mean/F");
  //fEventTree->Branch("mu_pida_med", &mu_pida_med,"mu_pida_med/F");
  //fEventTree->Branch("mu_kde", &mu_kde,"mu_kde/F");
  //fEventTree->Branch("mu_trm_dqdx", &mu_trm_dqdx,"mu_trm_dqdx/F");
  //fEventTree->Branch("mu_trm_dedx", &mu_trm_dedx,"mu_trm_dedx/F");
  //fEventTree->Branch("mu_dedx",mu_dedx,"mu_dedx[3][3000]/F");
  //fEventTree->Branch("mu_rr",mu_rr,"mu_rr[3][3000]/F");
  //fEventTree->Branch("mu_mom_process",&mu_mom_process);
  //fEventTree->Branch("cc_mu_trkid", &cc_mu_trkid,"cc_mu_trkid/I");
  //fEventTree->Branch("cc_mu_tlen", &cc_mu_tlen,"cc_mu_tlen/F");
  //fEventTree->Branch("cc_mu_phi", &cc_mu_phi,"cc_mu_phi/F");
  //fEventTree->Branch("cc_mu_theta", &cc_mu_theta,"cc_mu_theta/F");
  //fEventTree->Branch("cc_mu_range", &cc_mu_range,"cc_mu_range/F");
  //fEventTree->Branch("cc_mu_KE", &cc_mu_KE,"cc_mu_KE/F");
  //fEventTree->Branch("cc_mu_hit", &cc_mu_hit,"cc_mu_hit/I");
  //fEventTree->Branch("cc_mu_large_dedx", &cc_mu_large_dedx,"cc_mu_large_dedx/F");
  //fEventTree->Branch("cc_mu_small_dedx", &cc_mu_small_dedx,"cc_mu_small_dedx/F");
  //fEventTree->Branch("cc_dis_vtx", &cc_dis_vtx,"cc_dis_vtx/F");
  //fEventTree->Branch("cc_mu_pdg", &cc_mu_pdg,"cc_mu_pdg/I");
  //fEventTree->Branch("cc_mu_chi_p", &cc_mu_chi_p,"cc_mu_chi_p/F");
  //fEventTree->Branch("cc_mu_chi_k", &cc_mu_chi_k,"cc_mu_chi_k/F");
  //fEventTree->Branch("cc_mu_chi_pi", &cc_mu_chi_pi,"cc_mu_chi_pi/F");
  //fEventTree->Branch("cc_mu_chi_mu", &cc_mu_chi_mu,"cc_mu_chi_mu/F");
  //fEventTree->Branch("cc_mu_p_max", &cc_mu_p_max,"cc_mu_p_max/F");
  //fEventTree->Branch("cc_mu_k_max", &cc_mu_k_max,"cc_mu_k_max/F");
  //fEventTree->Branch("cc_mu_pi_max", &cc_mu_pi_max,"cc_mu_pi_max/F");
  //fEventTree->Branch("cc_mu_mu_max", &cc_mu_mu_max,"cc_mu_mu_max/F");
  //fEventTree->Branch("cc_mu_mip_max", &cc_mu_mip_max,"cc_mu_mip_max/F");
  //fEventTree->Branch("cc_mu_L1_ratio", &cc_mu_L1_ratio,"cc_mu_L1_ratio/F");
  //fEventTree->Branch("cc_mu_LL1", &cc_mu_LL1,"cc_mu_LL1/F");
  //fEventTree->Branch("cc_mu_L2_ratio", &cc_mu_L2_ratio,"cc_mu_L2_ratio/F");
  //fEventTree->Branch("cc_mu_LL2", &cc_mu_LL2,"cc_mu_LL2/F");
  //fEventTree->Branch("cc_mu_Lp_ratio", &cc_mu_Lp_ratio,"cc_mu_Lp_ratio/F");
  //fEventTree->Branch("cc_mu_LLp", &cc_mu_LLp,"cc_mu_LLp/F");
  //fEventTree->Branch("cc_mu_Lk_ratio", &cc_mu_Lk_ratio,"cc_mu_Lk_ratio/F");
  //fEventTree->Branch("cc_mu_LLk", &cc_mu_LLk,"cc_mu_LLk/F");
  //fEventTree->Branch("cc_mu_pida_mean", &cc_mu_pida_mean,"cc_mu_pida_mean/F");
  //fEventTree->Branch("cc_mu_pida_med", &cc_mu_pida_med,"cc_mu_pida_med/F");
  //fEventTree->Branch("cc_mu_kde", &cc_mu_kde,"cc_mu_kde/F");
  //fEventTree->Branch("cc_mu_trm_dqdx", &cc_mu_trm_dqdx,"cc_mu_trm_dqdx/F");
  //fEventTree->Branch("cc_mu_trm_dedx", &cc_mu_trm_dedx,"cc_mu_trm_dedx/F");
  //fEventTree->Branch("longest_trkid", &longest_trkid,"longest_trkid/I");
  //fEventTree->Branch("longest_trklen", &longest_trklen,"longest_trklen/F");
  //fEventTree->Branch("Tr_pri_mu_pdg", &Tr_pri_mu_pdg,"Tr_pri_mu_pdg/I");
  //fEventTree->Branch("pri_Mu_is", &pri_Mu_is,"pri_Mu_is/I");
  //fEventTree->Branch("Tr_pri_st_k_is", &Tr_pri_st_k_is,"Tr_pri_st_k_is/I");
  //fEventTree->Branch("Tr_K_Inelas", &Tr_K_Inelas,"Tr_K_Inelas/I");
  //fEventTree->Branch("Tr_k_plen", &Tr_k_plen,"Tr_k_plen/F");
  //fEventTree->Branch("Tr_k_theta", &Tr_k_theta,"Tr_k_theta/F");
  //fEventTree->Branch("Tr_k_phi", &Tr_k_phi,"Tr_k_phi/F");
  //fEventTree->Branch("Tr_dec_mu_is", &Tr_dec_mu_is,"Tr_dec_mu_is/I");
  //fEventTree->Branch("Tr_dec_mu_pi_pdg", &Tr_dec_mu_pi_pdg,"Tr_dec_mu_pi_pdg/I");
  //fEventTree->Branch("Tr_mu_plen", &Tr_mu_plen,"Tr_mu_plen/F");
  //fEventTree->Branch("Tr_k_endE", &Tr_k_endE,"Tr_k_endE/F");
  //fEventTree->Branch("Tr_mu_theta", &Tr_mu_theta,"Tr_mu_theta/F");
  //fEventTree->Branch("Tr_mu_phi", &Tr_mu_phi,"Tr_mu_phi/F");
  //fEventTree->Branch("Tr_k_inTPC", &Tr_k_inTPC,"Tr_k_inTPC/I");
  //fEventTree->Branch("Tr_mu_inTPC", &Tr_mu_inTPC,"Tr_mu_inTPC/I");
  //fEventTree->Branch("Tr_k_in_5_TPC", &Tr_k_in_5_TPC,"Tr_k_in_5_TPC/I");
  //fEventTree->Branch("Tr_k_in_CC_TPC", &Tr_k_in_CC_TPC,"Tr_k_in_CC_TPC/I");
  //fEventTree->Branch("Tr_mu_in_5_TPC", &Tr_mu_in_5_TPC,"Tr_mu_in_5_TPC/I");
  //fEventTree->Branch("Tr_mu_in_CC_TPC", &Tr_mu_in_CC_TPC,"Tr_mu_in_CC_TPC/I");
  //fEventTree->Branch("Tr_kmu_open_ang", &Tr_kmu_open_ang,"Tr_kmu_open_ang/F");
  //fEventTree->Branch("vtx_5cm_mult", &vtx_5cm_mult,"vtx_5cm_mult/I");
  //fEventTree->Branch("k_start_dedx", &k_start_dedx,"k_start_dedx/F");
  //fEventTree->Branch("k_end_dedx", &k_end_dedx,"k_end_dedx/F");
  //fEventTree->Branch("mu_start_dedx", &mu_start_dedx,"mu_start_dedx/F");
  //fEventTree->Branch("mu_end_dedx", &mu_end_dedx,"mu_end_dedx/F");
  fEventTree->Branch("cut_1", &cut_1,"cut_1/I");
  fEventTree->Branch("cut_2", &cut_2,"cut_2/I");
  fEventTree->Branch("cut_3", &cut_3,"cut_3/I");
  fEventTree->Branch("cut_4", &cut_4,"cut_4/I");
  fEventTree->Branch("cut_5", &cut_5,"cut_5/I");
  fEventTree->Branch("cut_6", &cut_6,"cut_6/I");
  fEventTree->Branch("cut_7", &cut_7,"cut_7/I");
  fEventTree->Branch("cut_8", &cut_8,"cut_8/I");
  fEventTree->Branch("cut_9", &cut_9,"cut_9/I");
  fEventTree->Branch("cut_10", &cut_10,"cut_10/I");
  fEventTree->Branch("cut_11", &cut_11,"cut_11/I");
  fEventTree->Branch("cut_12", &cut_12,"cut_12/I");
  fEventTree->Branch("cut_13", &cut_13,"cut_13/I");

  fSubrunTree = tfs->make<TTree>("subruns", "SubRun Tree");
  fSubrunTree->Branch("run", &m_run, "run/i");
  fSubrunTree->Branch("subRun", &m_subrun, "subRun/i");
  //if (!m_isData)
  fSubrunTree->Branch("pot", &m_pot, "pot/F");
    
}

void CCKaonAnalyzerRebuild::analyze( const art::Event& evt){
  reset();  
  run=evt.run();
  subrun=evt.subRun();
  event=evt.id().event(); 
  cout << "Run " << run;
  cout << " Subrun " << subrun;
  cout << " Event " << event << endl;

  art::Handle< std::vector<simb::MCTruth> > mctruths;
  std::vector< art::Ptr<simb::MCParticle> > ptList;

  
  // get MCTruth     
  evt.getByLabel("generator", mctruths);
  if (mctruths->size()!=1) {
    //std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
    return;
  }
  simb::MCTruth mctruth = mctruths->at(0);
  
  
  // get MCParticles
  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
    art::fill_ptr_vector(ptList, mcParticleHandle); 
  }

  
  //-----------------------------


  if (isMC) {

    // get MCTruth
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
      return;
    }
    simb::MCTruth mctruth = mctruths->at(0);

    // true neutrino information
    true_nu_energy = mctruth.GetNeutrino().Nu().E();
    true_nu_pdg = mctruth.GetNeutrino().Nu().PdgCode();
    true_nu_ccnc = mctruth.GetNeutrino().CCNC();//0=CC,1=NC
    true_nu_mode = mctruth.GetNeutrino().Mode();//0=QE,1=RES,2=DIS
    true_nu_vtx_x = mctruth.GetNeutrino().Nu().Vx();
    true_nu_vtx_y = mctruth.GetNeutrino().Nu().Vy();
    true_nu_vtx_z = mctruth.GetNeutrino().Nu().Vz();
    true_nu_vtx_inTPC = isInsideVolume("TPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",mctruth.GetNeutrino().Nu().Position().Vect());
    true_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",mctruth.GetNeutrino().Nu().Position().Vect());

    true_lepton_pdg = mctruth.GetNeutrino().Lepton().PdgCode();
    true_lepton_p = mctruth.GetNeutrino().Lepton().P();
    true_lepton_ke = mctruth.GetNeutrino().Lepton().E()-mctruth.GetNeutrino().Lepton().Mass();
    true_lepton_theta = mctruth.GetNeutrino().Lepton().Momentum().Theta();
    true_lepton_costheta = mctruth.GetNeutrino().Lepton().Momentum().CosTheta();
    true_lepton_phi = mctruth.GetNeutrino().Lepton().Momentum().Phi();
    TVector3 true_lepton_pvector = mctruth.GetNeutrino().Lepton().Momentum().Vect();

    // find the the highest momentum k+
    double true_kaon_maxp = -1;
    int true_kaon_maxp_id = -1;
    true_nkaons = 0;
    true_nprimary = 0;
    true_nhyperons = 0;

    int nparts = 0; 
    for (auto const& pPart : ptList) {
      if(pPart->Process()=="primary" && int(mctruths->at(0).Origin())==1){
	TLorentzVector mcstart, mcend;
	unsigned int pstarti, pendi;
	true_length[nparts] = length(*pPart, mcstart, mcend, pstarti, pendi);
	true_p[nparts]  = pPart->P();
	true_ke[nparts] = pPart->E()-pPart->Mass();
	true_theta[nparts] = pPart->Momentum().Theta();
	true_costheta[nparts] = pPart->Momentum().CosTheta();
	true_phi[nparts] = pPart->Momentum().Phi();
	true_ccmuon_angle[nparts] = pPart->Momentum().Angle(true_lepton_pvector);
	true_ccmuon_cosangle[nparts] = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
	//true_pvector = pPart->Momentum().Vect();
	true_end_ke[nparts] = pPart->EndE()-pPart->Mass();
	true_end_x[nparts] = pPart->EndX();
	true_end_y[nparts] = pPart->EndY();
	true_end_z[nparts] = pPart->EndZ();
	true_end_inTPC[nparts] = isInsideVolume("TPC",pPart->EndPosition().Vect());
	true_end_in5cmTPC[nparts] = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
	true_end_inCCInclusiveTPC[nparts] = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
      }

      if (pPart->Process()=="primary" && pPart->PdgCode()==321) {
	Process = pPart->Process();
	EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          if (true_kaon_maxp < pPart->P()) {
            true_kaon_maxp = pPart->P();
            true_kaon_maxp_id = pPart->TrackId();
            true_nkaons++;
          }
        }  
      }

      if (pPart->Process()=="primary") {
	//Process = pPart->Process();
	//EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
            true_nprimary++;
        }  
      }

      if (pPart->Process()=="primary" && (pPart->PdgCode()==3212 || pPart->PdgCode()==3122 || pPart->PdgCode()==3222 || pPart->PdgCode()==3112)) {
	//Process = pPart->Process();
	//EndProcess = pPart->EndProcess();
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          //if (true_hyperon_maxp < pPart->P()) {
            //true_hyperon_maxp = pPart->P();
            //true_hyperon_maxp_id = pPart->TrackId();
            true_nhyperons++;
	    //}
        }  
      }
      nparts++;
    }

    // true signal interaction: nu CC with at least one k+ // selection with true info
    if (true_nu_pdg==14 && true_nu_ccnc==0 && true_nkaons>0) {

      vector<int> decay_daughters;
      vector<int> decay_muplus_id;
      vector<int> decay_piplus_id;
      vector<int> inelastic_daughters;
      vector<int> inelastic_kaplus_id;
      TVector3 true_kaon_pvector;
      TVector3 true_dau_muon_pvector;
      TVector3 true_dau_pip_pvector;
      TVector3 true_dau_pin_pvector;

      //cout << "true_kaon_maxp_id: " << true_kaon_maxp_id << endl;

      for (auto const& pPart : ptList) {

        // primary kaon+ (w/ highest momentum)
        if (pPart->TrackId()==true_kaon_maxp_id) {
          TLorentzVector mcstart, mcend;
          unsigned int pstarti, pendi;
          true_kaon_length = length(*pPart, mcstart, mcend, pstarti, pendi);
          true_kaon_p  = pPart->P();
          true_kaon_ke = pPart->E()-pPart->Mass();
          true_kaon_theta = pPart->Momentum().Theta();
          true_kaon_costheta = pPart->Momentum().CosTheta();
          true_kaon_phi = pPart->Momentum().Phi();
          true_kaon_ccmuon_angle = pPart->Momentum().Angle(true_lepton_pvector);
          true_kaon_ccmuon_cosangle = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
          true_kaon_pvector = pPart->Momentum().Vect();
          true_kaon_end_ke = pPart->EndE()-pPart->Mass();

	  true_kaon_start_x = pPart->Vx();
	  true_kaon_start_y = pPart->Vy();
	  true_kaon_start_z = pPart->Vz();
          true_kaon_end_x = pPart->EndX();
          true_kaon_end_y = pPart->EndY();
          true_kaon_end_z = pPart->EndZ();
          true_kaon_end_inTPC = isInsideVolume("TPC",pPart->EndPosition().Vect());
          true_kaon_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
          true_kaon_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
        }

        // find kaon+ daughters
        else if (pPart->Mother()==true_kaon_maxp_id) {

	  //cout << "pPart->TrackId(): " << pPart->TrackId() << ", it's PDG: " << pPart->PdgCode() << endl;

          // decay daughters
          if (pPart->Process()=="Decay") {
            decay_daughters.push_back(pPart->PdgCode());
	    //cout << "decay_daughters: " << pPart->PdgCode() << endl;
            if (pPart->PdgCode()==-13) {
              decay_muplus_id.push_back(pPart->TrackId());

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_muon_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_muon_p  = pPart->P();
	      true_dau_muon_ke = pPart->E()-pPart->Mass();
	      true_dau_muon_theta = pPart->Momentum().Theta();
	      true_dau_muon_costheta = pPart->Momentum().CosTheta();
	      true_dau_muon_phi = pPart->Momentum().Phi();
	      true_dau_muon_pvector = pPart->Momentum().Vect();
	      true_dau_muon_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_muon_start_x = pPart->Vx();
	      true_dau_muon_start_y = pPart->Vy();
	      true_dau_muon_start_z = pPart->Vz();
	      true_dau_muon_end_x = pPart->EndX();
	      true_dau_muon_end_y = pPart->EndY();
	      true_dau_muon_end_z = pPart->EndZ();

            }
            else if (pPart->PdgCode()==211) {
              decay_piplus_id.push_back(pPart->TrackId());

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_pip_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_pip_p  = pPart->P();
	      true_dau_pip_ke = pPart->E()-pPart->Mass();
	      true_dau_pip_theta = pPart->Momentum().Theta();
	      true_dau_pip_costheta = pPart->Momentum().CosTheta();
	      true_dau_pip_phi = pPart->Momentum().Phi();
	      true_dau_pip_pvector = pPart->Momentum().Vect();
	      true_dau_pip_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_pip_start_x = pPart->Vx();
	      true_dau_pip_start_y = pPart->Vy();
	      true_dau_pip_start_z = pPart->Vz();
	      true_dau_pip_end_x = pPart->EndX();
	      true_dau_pip_end_y = pPart->EndY();
	      true_dau_pip_end_z = pPart->EndZ();

            }
	    else if (pPart->PdgCode()==111) {

	      TLorentzVector mcstart, mcend;
	      unsigned int pstarti, pendi;
	      true_dau_pin_length = length(*pPart, mcstart, mcend, pstarti, pendi);
	      true_dau_pin_p  = pPart->P();
	      true_dau_pin_ke = pPart->E()-pPart->Mass();
	      true_dau_pin_theta = pPart->Momentum().Theta();
	      true_dau_pin_costheta = pPart->Momentum().CosTheta();
	      true_dau_pin_phi = pPart->Momentum().Phi();
	      true_dau_pin_pvector = pPart->Momentum().Vect();
	      true_dau_pin_end_ke = pPart->EndE()-pPart->Mass();
	      
	      true_dau_pin_start_x = pPart->Vx();
	      true_dau_pin_start_y = pPart->Vy();
	      true_dau_pin_start_z = pPart->Vz();
	      true_dau_pin_end_x = pPart->EndX();
	      true_dau_pin_end_y = pPart->EndY();
	      true_dau_pin_end_z = pPart->EndZ();

    cout << "true_track_length: " << true_dau_pin_length << endl;

	    }
            //cout << "Decay " << pPart->PdgCode() << endl;
          }

          // inelastic interaction daughters
          if (pPart->Process()=="kaon+Inelastic") {
            inelastic_daughters.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==321) {
              inelastic_kaplus_id.push_back(pPart->TrackId());
            }
            //cout << "Inelastic " << pPart->PdgCode() << endl;
          }

        }

      }//MC particles loop

      true_kaon_ndaughters = decay_daughters.size() + inelastic_daughters.size();
      true_kaon_ndaughters_decay = decay_daughters.size();
      true_kaon_ndaughters_inelastic = inelastic_daughters.size();

      int n_decay_muplus     = decay_muplus_id.size();
      int n_decay_piplus     = decay_piplus_id.size();
      int n_decay_numu       = count(decay_daughters.begin(),decay_daughters.end(),14);
      int n_decay_pi0        = count(decay_daughters.begin(),decay_daughters.end(),111);
      int n_inelastic_kaplus = inelastic_kaplus_id.size();
      int n_inelastic_piplus = count(inelastic_daughters.begin(),inelastic_daughters.end(),211);
      int n_inelastic_proton = count(inelastic_daughters.begin(),inelastic_daughters.end(),2212);

      true_kaon_ndecmup = n_decay_muplus;
      true_kaon_ndecpip = n_decay_piplus;
      true_kaon_ninekap = n_inelastic_kaplus;
      true_kaon_ninepip = n_inelastic_piplus;
      true_kaon_ninepro = n_inelastic_proton;
    
      // define kaon end process
      int daughter_id = -1;
      if (decay_daughters.size()) {
        //kaon+ -> muon+ numu
        if (decay_daughters.size()==2 && n_decay_muplus==1 && n_decay_numu==1) {
          true_kaon_end_process = 0;
          daughter_id = decay_muplus_id.front();
        }
        //kaon+ -> pion+ + pi0
        else if (decay_daughters.size()==2 && n_decay_piplus==1 && n_decay_pi0==1) {
          true_kaon_end_process = 1;
          daughter_id = decay_piplus_id.front();
        }
        //kaon+ -> muon+ nopion+
        else if (n_decay_muplus==1 && n_decay_piplus==0) {
          true_kaon_end_process = 2;
          daughter_id = decay_muplus_id.front();
        }
        //kaon+ -> pion+ and nomuon+
        else if (n_decay_muplus==0 && n_decay_piplus==1) {
          true_kaon_end_process = 3;
          daughter_id = decay_piplus_id.front();
        }
        //other decays
        else {
          true_kaon_end_process = 4;
        }
      }

      if (inelastic_daughters.size()) {
        //inelasitic with kaon+
        if (n_inelastic_kaplus==1) {
          true_kaon_end_process = 5;
          daughter_id = inelastic_kaplus_id.front();
        }
        //other inelastic
        else {
          true_kaon_end_process = 6;
        }
      }

      //this shouldn't happen but it does! rarely
      if (decay_daughters.size() && inelastic_daughters.size()) {
        // -- cout << "Primary kaon+ decay and inelastic at the same time!!!" << endl;
        true_kaon_end_process = 7;
      }

      cout << "True kaon end process " << true_kaon_end_process << endl;

      //check kaon daughter if it exists
      //only if there is one pion+, one muon+ or one kaon+
      if (daughter_id!=-1) {
        for (auto const& pPart : ptList) {
      
          // kaon+ daughter
          if (pPart->TrackId()==daughter_id) {
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            true_kaon_daughter_length = length(*pPart, mcstart, mcend, pstarti, pendi);
            true_kaon_daughter_p = pPart->P();
            true_kaon_daughter_ke = pPart->E()-pPart->Mass();
            true_kaon_daughter_theta = pPart->Momentum().Theta();
            true_kaon_daughter_costheta = pPart->Momentum().CosTheta();
            true_kaon_daughter_angle = pPart->Momentum().Angle(true_kaon_pvector);
            true_kaon_daughter_cosangle = TMath::Cos(pPart->Momentum().Angle(true_kaon_pvector));
            true_kaon_daughter_pdg = pPart->PdgCode();

	    true_kaon_daughter_start_x = pPart->Vx();
	    true_kaon_daughter_start_y = pPart->Vy();
	    true_kaon_daughter_start_z = pPart->Vz();
            true_kaon_daughter_end_x = pPart->EndX();
            true_kaon_daughter_end_y = pPart->EndY();
            true_kaon_daughter_end_z = pPart->EndZ();
            true_kaon_daughter_end_inTPC = isInsideVolume("TPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_in5cmTPC = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
            break;
          }
      
        }//MC particles loop
      }//kaon daughter exists

    }//is true signal
    
  }//isMC

  // Check if event passed the NuCC inclusive filter
  reco_nu_cc_filter = false;
  string process(isMC ? "OverlayFiltersPostStage2" : "DataFiltersPostStage2");
  art::InputTag trigResInputTag("TriggerResults","",process.data()); // the last is the name of process where the filters were run
  art::ValidHandle<art::TriggerResults> trigRes = evt.getValidHandle<art::TriggerResults>(trigResInputTag);
  fhicl::ParameterSet pset;
  if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
  std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
  if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
  for (size_t itp=0;itp<trigRes->size();itp++) {
    //cout << "Filter name " << trigger_path_names.at(itp) << " decision=" << trigRes->at(itp).accept() << endl;
    if (trigger_path_names.at(itp)=="NuCC") {
      reco_nu_cc_filter = trigRes->at(itp).accept();
    }
  }


  art::Handle< vector< recob::PFParticle > > pfphandle;
  vector< art::Ptr< recob::PFParticle > > pfpvector;
  vector< art::Ptr< recob::PFParticle > > pfpvectorPandora;

  if(evt.getByLabel(fPFParticleLabel,pfphandle)) art::fill_ptr_vector(pfpvector,pfphandle);
  if(!pfpvector.size()) return;

  art::FindManyP< larpandoraobj::PFParticleMetadata > pfpmeta(pfphandle,evt,"pandora");
  art::FindManyP<recob::Slice> nuslice(pfpvector, evt, "pandora");

  size_t pfpkey = 999999;

  for(const art::Ptr< recob::PFParticle > &pfp : pfpvector){

    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfpmeta.at(pfp.key()));

    if (pfParticleMetadataList.empty()) continue;


    for (const art::Ptr<larpandoraobj::PFParticleMetadata> &pfpmd : pfParticleMetadataList ) {
      auto pfParticlePropertiesMap = pfpmd->GetPropertiesMap();
      //for (const larpandoraobj::PFParticleMetadata pfpmd : pfParticleMetadataList ) {
      //auto pfParticlePropertiesMap = pfpmd.GetPropertiesMap();
      if (pfParticlePropertiesMap.empty()) continue; //cout << pfParticlePropertiesMap.size() << endl;
      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	if(!(it->first == "TrackScore")) continue;
	//cout << "PFParticle properties: " << it->first << " Key: " << pfpmd.key() << " Value: " << it->second << endl;
	pfpkey = pfp.key();
	pfpvectorPandora.push_back(pfp);
      }
    }// for all metadata items in the particle metadata
    if(pfp.key() == pfpkey){
      /*
      cout << "PFP ID: " <<  pfp->Self() <<endl;
      cout << "PFP IsPrimary: " << pfp->IsPrimary() <<endl;
      cout << "PFP parent: " << pfp->Parent() <<endl;
      cout << "PFP PDG: " << pfp->PdgCode() <<endl;
      cout << "PFP key: " << pfp.key() << "=====================\n" <<endl;
      */
    }

    auto slices = nuslice.at(pfp.key());

  }// for entries in list


  if (!reco_nu_cc_filter) {
    // -- std::cout << "Event didn't pass NuCC filter" << std::endl;
    fEventTree->Fill();
    return;
  }

  // -- std::cout << "Event passed NuCC filter" << std::endl;
  cut_1=1;

  // Collect all recontructed particles
  art::Handle<std::vector<recob::PFParticle>> pfparticles;
  evt.getByLabel(m_pfp_producer, pfparticles);

  if (pfparticles->size()==0) {
    std::cout << "No PFParticles found" << std::endl;
    fEventTree->Fill();
    return;
  }
  // -- std::cout << "Number of PFParticles " << pfparticles->size() << std::endl;
  cut_2=1;	  

  // Get PFParticle associations
  art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
  if(!pfp_muon_assn.isValid()){
    // -- cout << "PFParticle-T0 associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_3=1;	  

  art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
  if(!pfparticleTrackAssn.isValid()){
    // -- cout << "PFParticle-Track associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_4=1;	  

  art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
  if(!pfparticleVertexAssn.isValid()){
    // -- cout << "PFParticle-Vertex associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  cut_5=1;	  

  // Find recontructed neutrino (there should be one)
  lar_pandora::PFParticleVector pfneutrinos(0);
  for (unsigned int i=0; i<pfparticles->size(); ++i) {

    art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

    if (pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
      pfneutrinos.push_back(pfparticle);
    }

  }

  if (pfneutrinos.size() != 1) {
    // --cout << "Number of neutrinos is not one" << endl;
    cout << "Number of neutrinos " << pfneutrinos.size() << endl;
    fEventTree->Fill();
    return;
  }
  cut_6=1;	  

  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

  reco_nu_vtx_x = pfparticleVertexAssn.at(pfnu.key()).front()->position().X();
  reco_nu_vtx_y = pfparticleVertexAssn.at(pfnu.key()).front()->position().Y();
  reco_nu_vtx_z = pfparticleVertexAssn.at(pfnu.key()).front()->position().Z();
  //--cout << "Neutrino vertex (x,y,z) = " << reco_nu_vtx_x << ", " << reco_nu_vtx_y << ", " << reco_nu_vtx_z << endl;
  reco_nu_vtx_inTPC = isInsideVolume("TPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
  reco_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
  reco_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);

  // Find CC muon and daughters
  lar_pandora::PFParticleVector pfmuons(0);
  vector<int> reco_nu_daughters_id(0);

  for (unsigned int i=0; i<pfparticles->size(); ++i) {

    art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

    // look at particles with neutrino parent and one associated track
    if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1) {

      art::Ptr<recob::Track> track = pfparticleTrackAssn.at(i).front();

      reco_nu_daughters_id.push_back(track.key());

      // CC muon has a T0 associated
      if (pfp_muon_assn.at(i).size()==1) {
        pfmuons.push_back(pfparticle);
      }

    }

  }

  reco_nu_ndaughters = reco_nu_daughters_id.size();
  //--std::cout << "Number of neutrino daughters with one associated track " << reco_nu_ndaughters << std::endl;
  reco_nu_cc_nmue = pfmuons.size();

  if (pfmuons.size()!=1) {
    //--std::cout << "Number of CC inclusive muons is not 1: " << pfmuons.size() << std::endl;
    fEventTree->Fill();
    return;
  }

  art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
  art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();

  reco_ccmu_vtx_x = pfparticleVertexAssn.at(pfmuon.key()).front()->position().X();
  reco_ccmu_vtx_y = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Y();
  reco_ccmu_vtx_z = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Z();
  reco_ccmu_vtx_inTPC = isInsideVolume("TPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);

  // get hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitlist, hitListHandle);
  }

  // get tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
    art::fill_ptr_vector(trackList, trackListHandle);
  }

  // get showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerList;
  if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
    art::fill_ptr_vector(showerList, showerListHandle);
  }

  // get spacepoints
  art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
  std::vector< art::Ptr<recob::SpacePoint> > spacepointVector;
  if(evt.getByLabel(fPandoraLabel,spacepointHandle)){
    art::fill_ptr_vector(spacepointVector,spacepointHandle);
  }
  
  art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
  art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
  if(!hits_from_tracks.isValid()){
    fEventTree->Fill();
    return;
  }

  art::FindOneP<recob::Hit> findSPToHit(spacepointVector, evt, fSpacePointproducer); 
  art::FindManyP<recob::Hit> findTrackToHit(trackList, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit> findShowerToHit(showerList, evt, fShowerModuleLabel);
  
  std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> hitToSpacePointMap;
  std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> spacepointToHitMap;

  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
  std::vector< art::Ptr<recob::SpacePoint> > spacepointFromRecoObject;
  std::vector< art::Ptr<recob::SpacePoint> > spacepointFromRecoObjectOld;
  std::vector< art::Ptr<recob::SpacePoint> > spacepointFromMu;
  std::vector< art::Ptr<recob::SpacePoint> > spacepointFromPi;
  std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
  
  //fill the SP-Hit map
  for (unsigned int itrk=0; itrk < trackList.size(); ++itrk) {
    
    const art::Ptr<recob::Track> track = trackList.at(itrk);

    if (track.key()==trkmuon.key()) continue;

    //skip primary track
    bool skip = false;
    for (int i_nutrk=0; i_nutrk<reco_nu_ndaughters; i_nutrk++) {
      if (int(track.key())==reco_nu_daughters_id[i_nutrk]) {
	skip=true;
	break;
      }
    }
    if (skip) continue;    

    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track.key());

    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
    
    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
      
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());
      
      if(spacepoint_vec.size()!=1) continue;
      art::Ptr<recob::SpacePoint> spacepoint = spacepoint_vec.at(0);
      art::Ptr<recob::Hit> hit = hits_from_track.at(i_h);
      
      spacepointFromRecoObject.push_back(spacepoint);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
      
      simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
      //std::cout << "hit.key(): " << hit.key() << endl;
      if(!mcparticle) continue;
      //if(hit.key()==1509) cout << "hit 1509 has pdg of " << mcparticle->PdgCode() << endl;
      //std::cout << "(mcparticle->PdgCode(): " << mcparticle->PdgCode() << endl;
      if(mcparticle->PdgCode()==-13) spacepointFromMu.push_back(spacepoint);
      if(mcparticle->PdgCode()==211) spacepointFromPi.push_back(spacepoint);

    }
    
  }
  
  
  for (unsigned int ishw=0; ishw < showerList.size(); ++ishw) {
    
    const art::Ptr<recob::Shower> shower = showerList.at(ishw);
    
    std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_showers.at(shower.key());
    
    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
    
    for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++){
      
      spacepoint_vec.clear();
      spacepoint_vec = spacepoint_per_hit.at(hits_from_shower[i_h].key());
      
      if(spacepoint_vec.size()!=1) continue;
      art::Ptr<recob::SpacePoint> spacepoint = spacepoint_vec.at(0);
      art::Ptr<recob::Hit> hit = hits_from_shower.at(i_h);
      
      spacepointFromRecoObject.push_back(spacepoint);
      spacepointToHitMap[spacepoint] = hit;
      hitToSpacePointMap[hit] = spacepoint;
      
      simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
      if(!mcparticle) continue;
      if(mcparticle->PdgCode()==-13) spacepointFromMu.push_back(spacepoint);
      if(mcparticle->PdgCode()==211) spacepointFromPi.push_back(spacepoint);
      
    }
    
  }
  
  
  // get track associations
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  if(!fmcal.isValid()){
    fEventTree->Fill();
    return;
  }

  art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fPIDLabel);
  if(!trackPIDAssn.isValid()){
    fEventTree->Fill();
    return;
  }



  //selection cuts
  std::vector<int> kaon_can_trkID;
  std::vector<int> muon_can_trkID;

  std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>> assocMCPartt;
  if (isMC) {
    assocMCPartt = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> (new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, evt, fHitTruthAssns));
  }

  // loop over nu daughters
  trkf::TrackMomentumCalculator trkmom;

  int ntracks = 0;

  //int NTracks=trackList.size();
  //int NShowers=showerList.size();

  for (int i=0; i<reco_nu_ndaughters; i++) {

    n_recoRebDauTracks[i] = 0;

    art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);

    // skip cc muon track
    if (ptrack.key()==trkmuon.key()) continue;

    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());

    simb::MCParticle const* mcparticle = truthMatchTrack(hits_from_track, particles_per_hit);
    //if(mcparticle && mcparticle->PdgCode()==321) cout << "This is Primary Kaon" << endl;
    //if(mcparticle->PdgCode()!=321) continue;
    //if(true_kaon_end_process!=0) continue;
    if(mcparticle) std::cout << mcparticle->PdgCode() << ": primary mcparticle->PdgCode()" << endl;
    //if(mcparticle) recoprimarttrack_pdg[itrk]mcparticle->PdgCode();    

    /*
    */

    /*
    // check track start and end
    TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
    TVector3 end(track.End().X(),track.End().Y(),track.End().Z());
    */
    
    ReconstructionOrchestrator orchestrator;
    ReconstructionOrchestrator orchestrator_cheatpi;
    ReconstructionOrchestrator orchestrator_cheatmu;

    orchestrator.runReconstruction(spacepointFromRecoObject, spacepointToHitMap, hitToSpacePointMap,ptrack, hits_from_track);
    //orchestrator.runReconstruction(spacepointFromMu, spacepointToHitMap, hitToSpacePointMap,ptrack, hits_from_track);

    std::vector<recob::Track> rebuildTrackList = orchestrator.getRebuildTrackList();
    std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists = orchestrator.getHitLists();

    for(unsigned int j=0; j < rebuildTrackList.size(); j++) {

      const std::vector<TVector3> peakDirectionVector = orchestrator.getPeakDirectionList();
      //if (!(j < peakDirectionVector.size())) break;
      best_peak_x[i][j] = peakDirectionVector.at(j).X();
      best_peak_y[i][j] = peakDirectionVector.at(j).Y();
      best_peak_z[i][j] = peakDirectionVector.at(j).Z();

      if(trackHitLists[j].empty()) continue;

      cout << "j: " << j << ", n_recoRebDauTracks[i]: " << n_recoRebDauTracks[i] << ", rebuildTrackList[j].Length(): " << rebuildTrackList[j].Length() << endl;
      rebdautrack_length[i][n_recoRebDauTracks[i]] = rebuildTrackList[j].Length();
      std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[j];

      //simb::MCParticle const* mcparticle;
      std::map<int,int> nhits_pdg_map;
      std::map<recob::Hit,int> hit_pdg_map;

      //simb::MCParticle const* mcparticle = truthMatchTrack(hits_from_track_rebuild,  particles_per_hit);
      //if(mcparticle) std::cout << mcparticle->PdgCode() << ": mcparticle->PdgCode()" << endl;

      for(unsigned int i_h=0; i_h<hits_from_track_rebuild.size(); i_h++){

	art::Ptr<recob::Hit> hit = hits_from_track_rebuild.at(i_h);

	simb::MCParticle const* mcparticle = truthMatchHit(hit, particles_per_hit);
	//std::cout << "hit.key(): " << hit.key() << endl;
	if(!mcparticle){
	  //cout << "THIS HAS NO MCPARTICLE!!!" << endl;
	  continue;
	}//else cout << "THIS HAS MCPARTICLE WITH PDG " << mcparticle->PdgCode()<< endl;
	//if(hit.key()==1509) cout << "hit 1509 has pdg of " << mcparticle->PdgCode() << endl;

 
	hit_pdg_map[(*hit)] = mcparticle->PdgCode();

	if(nhits_pdg_map.find(mcparticle->PdgCode()) == nhits_pdg_map.end()) nhits_pdg_map[mcparticle->PdgCode()] = 1;
	else nhits_pdg_map[mcparticle->PdgCode()] += 1;

      }

      vector<pair<int, int>> v;
      if(nhits_pdg_map.size()){
	for (map<int, int>::iterator it = nhits_pdg_map.begin(); it != nhits_pdg_map.end(); it++) {
	  v.push_back({ it->second, it->first });
	  //cout << it->second << " " << it->first << endl;
	}
	sort(v.rbegin(), v.rend());
	rebdautrack_pdg[i][n_recoRebDauTracks[i]] = v[0].second;
	std::cout<< "rebdau pdg: " <<  v[0].second << endl;
      }
      n_recoRebDauTracks[i]++;

    }
    
    if(true_kaon_end_process==0){
      cout << "this is cheat mu start" << endl;
      orchestrator_cheatmu.runReconstruction(spacepointFromMu, spacepointToHitMap, hitToSpacePointMap,ptrack, hits_from_track);
      cout << "this is cheat mu end" << endl;
      std::vector<recob::Track> rebuildTrackList_cheatmu = orchestrator_cheatmu.getRebuildTrackList(); 
      
      if(!rebuildTrackList_cheatmu.empty()){
        
	rebdautracktrue_length[i] = rebuildTrackList_cheatmu[0].Length();
        
	std::vector<TVector3> peakDirectionVector =  orchestrator_cheatmu.getPeakDirectionList();
        
	best_peak_x_true[i] = peakDirectionVector[0].X(); 
	best_peak_y_true[i] = peakDirectionVector[0].Y(); 
	best_peak_z_true[i] = peakDirectionVector[0].Z(); 
	
	orchestrator_cheatmu.runReconstruction(spacepointFromRecoObject, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track, peakDirectionVector);
	rebdautracktruedir_length[i] = orchestrator_cheatmu.getRebuildTrackList().at(0).Length();
      }
      
    }
    else if(true_kaon_end_process==1){
      cout << "this is cheat pi start" << endl;
      orchestrator_cheatmu.runReconstruction(spacepointFromPi, spacepointToHitMap, hitToSpacePointMap,ptrack, hits_from_track);
      cout << "this is cheat pi start" << endl;
      std::vector<recob::Track> rebuildTrackList_cheatpi = orchestrator_cheatpi.getRebuildTrackList(); 
      
      if(!rebuildTrackList_cheatpi.empty()){
        
	rebdautracktrue_length[i] = rebuildTrackList_cheatpi[0].Length();
        
	std::vector<TVector3> peakDirectionVector =  orchestrator_cheatpi.getPeakDirectionList();
        
	best_peak_x_true[i] = peakDirectionVector[0].X(); 
	best_peak_y_true[i] = peakDirectionVector[0].Y(); 
	best_peak_z_true[i] = peakDirectionVector[0].Z(); 
	
	orchestrator_cheatpi.runReconstruction(spacepointFromRecoObject, spacepointToHitMap, hitToSpacePointMap, ptrack, hits_from_track, peakDirectionVector);
	rebdautracktruedir_length[i] = orchestrator_cheatpi.getRebuildTrackList().at(0).Length();
      }
      
    }


    if (isMC) {
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      fillTrueMatching(hits_from_track, particles_per_hit, ntracks);
    }

    ntracks++;
    	
  }//NTrack loop
  
  reco_ntracks = ntracks;
    
  // find true matched particle for ccmuon
  if (isMC) {

    simb::MCParticle const* matched_mcparticle = NULL;
    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(trkmuon.key());
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle,evt,fHitTruthAssns);


    //spacepoint
    art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer); 
    //std::vector<recob::SpacePoint const*> spacepoint_vec;
    std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
      particle_vec.clear(); match_vec.clear();
      spacepoint_vec.clear();

      particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
      //spacepoint_per_hit.get(hits_from_track[i_h].key(), spacepoint_vec);

      //const std::vector< art::Ptr<recob::SpacePoint> > spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());
      spacepoint_vec = spacepoint_per_hit.at(hits_from_track[i_h].key());

      //for(size_t i_s=0; i_s<spacepoint_vec.size(); ++i_s){
	//cout << "3D XYZ " << spacepoint_vec[i_s]->XYZ()[0] << " " << spacepoint_vec[i_s]->XYZ()[1] << " " << spacepoint_vec[i_s]->XYZ()[2] << endl;
      //}  

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
        trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
        //trkide[ particle_vec[i_p]->TrackId() ] ++;
        tote += match_vec[i_p]->energy;
        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
          maxe = trkide[ particle_vec[i_p]->TrackId() ];
          matched_mcparticle = particle_vec[i_p];
        }
      }
    }

    if(matched_mcparticle){
      reco_ccmu_true_pdg = matched_mcparticle->PdgCode();
      const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
      reco_ccmu_true_origin=int(mc_truth->Origin());
      reco_ccmu_true_primary=matched_mcparticle->Process()=="primary";
      double x = matched_mcparticle->EndX();
      double y = matched_mcparticle->EndY();
      double z = matched_mcparticle->EndZ();
      reco_ccmu_true_end_inTPC = isInsideVolume("TPC",x,y,z);
      reco_ccmu_true_end_in5cmTPC = isInsideVolume("5cmTPC",x,y,z);
      reco_ccmu_true_end_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",x,y,z);
      TLorentzVector mcstart, mcend;
      unsigned int pstarti, pendi;
      reco_ccmu_true_length = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
    }
  }//isMC

  cout << "Filling tree ----------------" << endl;

  fEventTree->Fill();

} // end of analyze function

 /////////////////////////////////////////// Reset Function ///////////////////////////////

 void CCKaonAnalyzerRebuild::reset(){

   run=-9999;
   subrun=-9999;
   event=-9999;
   
   event_weight.clear();
   evtwgt_funcname.clear();
   evtwgt_weight.clear();
   evtwgt_nweight.clear();
   
   
   true_nu_energy=-9999;
   true_nu_pdg=-9999;
   true_nu_mode=-9999;
   true_nu_ccnc=-9999;
   true_nu_vtx_x=-9999;
   true_nu_vtx_y=-9999;
   true_nu_vtx_z=-9999;
   true_nu_vtx_inTPC=false;
   true_nu_vtx_in5cmTPC=false;
   true_nu_vtx_inCCInclusiveTPC=false;
   
   IsKaon=false;
   IsSingleKaon=false;
   IsAssociatedKaon=false;
   IsMuBR=false;
   IsPiBR=false;
   flag_prim_k=false;
   
   true_lepton_pdg=-9999;
   true_lepton_p=-9999;
   true_lepton_ke=-9999;
   true_lepton_theta=-9999;
   true_lepton_costheta=-9999;
   true_lepton_phi=-9999;
   
   true_nkaons=-9999;
   true_kaon_length=-9999;
   true_kaon_p=-9999;
   true_kaon_ke=-9999;
   true_kaon_theta=-9999;
   true_kaon_costheta=-9999;
   true_kaon_phi=-9999;
   true_kaon_ccmuon_angle=-9999;
   true_kaon_ccmuon_cosangle=-9999;
   true_kaon_end_process=-9999;
   true_kaon_end_ke=-9999;
   true_kaon_start_x=-9999;
   true_kaon_start_y=-9999;
   true_kaon_start_z=-9999;
   true_kaon_end_x=-9999;
   true_kaon_end_y=-9999;
   true_kaon_end_z=-9999;
   true_kaon_end_inTPC=false;
   true_kaon_end_in5cmTPC=false;
   true_kaon_end_inCCInclusiveTPC=false;
   
   true_dau_muon_length=-9999;
   true_dau_muon_p=-9999;
   true_dau_muon_ke=-9999;
   true_dau_muon_theta=-9999;
   true_dau_muon_costheta=-9999;
   true_dau_muon_phi=-9999;
   true_dau_muon_ccmuon_angle=-9999;
   true_dau_muon_ccmuon_cosangle=-9999;
   true_dau_muon_end_process=-9999;
   true_dau_muon_end_ke=-9999;
   true_dau_muon_start_x=-9999;
   true_dau_muon_start_y=-9999;
   true_dau_muon_start_z=-9999;
   true_dau_muon_end_x=-9999;
   true_dau_muon_end_y=-9999;
   true_dau_muon_end_z=-9999;
    
   true_dau_pip_length=-9999;
   true_dau_pip_p=-9999;
   true_dau_pip_ke=-9999;
   true_dau_pip_theta=-9999;
   true_dau_pip_costheta=-9999;
   true_dau_pip_phi=-9999;
   true_dau_pip_ccmuon_angle=-9999;
   true_dau_pip_ccmuon_cosangle=-9999;
   true_dau_pip_end_process=-9999;
   true_dau_pip_end_ke=-9999;
   true_dau_pip_start_x=-9999;
   true_dau_pip_start_y=-9999;
   true_dau_pip_start_z=-9999;
   true_dau_pip_end_x=-9999;
   true_dau_pip_end_y=-9999;
   true_dau_pip_end_z=-9999;
   
   
   true_dau_pin_length=-9999;
   true_dau_pin_p=-9999;
   true_dau_pin_ke=-9999;
   true_dau_pin_theta=-9999;
   true_dau_pin_costheta=-9999;
   true_dau_pin_phi=-9999;
   true_dau_pin_ccmuon_angle=-9999;
   true_dau_pin_ccmuon_cosangle=-9999;
   true_dau_pin_end_process=-9999;
   true_dau_pin_end_ke=-9999;
   true_dau_pin_start_x=-9999;
   true_dau_pin_start_y=-9999;
   true_dau_pin_start_z=-9999;
   true_dau_pin_end_x=-9999;
   true_dau_pin_end_y=-9999;
   true_dau_pin_end_z=-9999;


   for(int k=0; k<kMaxParticles; k++){ 
     true_length[k]=-9999;
     true_p[k]=-9999;
     true_ke[k]=-9999;
     true_theta[k]=-9999;
     true_costheta[k]=-9999;
     true_phi[k]=-9999;
     true_ccmuon_angle[k]=-9999;
     true_ccmuon_cosangle[k]=-9999;
     true_end_process[k]=-9999;
     true_end_ke[k]=-9999;
     true_end_x[k]=-9999;
     true_end_y[k]=-9999;
     true_end_z[k]=-9999;
     true_end_inTPC[k]=false;
     true_end_in5cmTPC[k]=false;
     true_end_inCCInclusiveTPC[k]=false;
   }
 
   true_kaon_ndaughters=-9999;
   true_kaon_ndaughters_decay=-9999;
   true_kaon_ndaughters_inelastic=-9999;
   true_kaon_ndecmup=-9999;
   true_kaon_ndecpip=-9999;
   true_kaon_ninekap=-9999;
   true_kaon_ninepip=-9999;
   true_kaon_ninepro=-9999;
   true_kaon_daughter_length=-9999;
   true_kaon_daughter_p=-9999;
   true_kaon_daughter_ke=-9999;
   true_kaon_daughter_theta=-9999;
   true_kaon_daughter_costheta=-9999;
   true_kaon_daughter_angle=-9999;
   true_kaon_daughter_cosangle=-9999;
   true_kaon_daughter_pdg=-9999;
   true_kaon_daughter_start_x=-9999;
   true_kaon_daughter_start_y=-9999;
   true_kaon_daughter_start_z=-9999;
   true_kaon_daughter_end_x=-9999;
   true_kaon_daughter_end_y=-9999;
   true_kaon_daughter_end_z=-9999;
   true_kaon_daughter_end_inTPC=false;
   true_kaon_daughter_end_in5cmTPC=false;
   true_kaon_daughter_end_inCCInclusiveTPC=false;
   
   true_nhyperons=-9999;
   
   
   cheat_pip_trkln = -9999;
   cheat_mup_trkln = -9999;

   reco_nu_cc_filter=false;
   
   reco_nu_vtx_x=-9999;
   reco_nu_vtx_y=-9999;
   reco_nu_vtx_z=-9999;
   reco_nu_vtx_inTPC=false;
   reco_nu_vtx_in5cmTPC=false;
   reco_nu_vtx_inCCInclusiveTPC=false;
   reco_nu_ndaughters=-9999;
   reco_nu_cc_nmue=-9999;
   
   reco_ccmu_vtx_x=-9999;
   reco_ccmu_vtx_y=-9999;
   reco_ccmu_vtx_z=-9999;
   reco_ccmu_vtx_inTPC=false;
   reco_ccmu_vtx_in5cmTPC=false;
   reco_ccmu_vtx_inCCInclusiveTPC=false;
   reco_ccmu_true_pdg=-9999;
   reco_ccmu_true_origin=-9999;
   reco_ccmu_true_primary=false;
   reco_ccmu_true_end_inTPC=false;
   reco_ccmu_true_end_in5cmTPC=false;
   reco_ccmu_true_end_inCCInclusiveTPC=false;
   reco_ccmu_true_length=-9999;
   
   //reco_track_dEdx.clear();
   //reco_track_ResRan.clear();
   
   //reco_track_daughter_dEdx.clear();
   //reco_track_daughter_ResRan.clear();
   //vector<Float_t> v;
   
   reco_ntracks=0;
   reco_nshowers=0;
   
   for(int k=0; k<kMaxTracks; k++){


     for(int l=0; l<kMaxParticles; l++){ 
       cheat_num_hits[k][l] = -9999;
       cheat_ini_hits[k][l] = -9999;
       cheat_closest_distance[k][l] = -9999;       

       cheat_peak_pdg[k][l] = -9999;
       cheat_peak_theta[k][l] = -9999;
       cheat_peak_phi[k][l] = -9999;
       best_peak_trkln[k][l] = -9999;

       rebdautrack_length[k][l] = -9999;
       rebdautrack_pdg[k][l] = -9999;

       best_peak_x[k][l] = -9999;
       best_peak_y[k][l] = -9999;
       best_peak_z[k][l] = -9999;
     }  
     n_recoRebDauTracks[k]=-9;
     rebdautracktrue_length[k] = -9999;
     rebdautracktruedir_length[k] = -9999;
     best_peak_x_true[k] = -9999;
     best_peak_y_true[k] = -9999;
     best_peak_z_true[k] = -9999;
         

     //m_reco_track_sliceId[k]=-9;
     //m_reco_track_is_nuslice[k]=-9;
     reco_track_start_x[k]=-9;
     reco_track_start_y[k]=-9;
     reco_track_start_z[k]=-9;
     reco_track_end_x[k]=-9;
     reco_track_end_y[k]=-9;
     reco_track_end_z[k]=-9;
     
     reco_track_distance[k]=-9;
     reco_track_nhits0[k]=-9;
     reco_track_nhits1[k]=-9;
     reco_track_nhits2[k]=-9;
     reco_track_kin0[k]=-9;
     reco_track_kin1[k]=-9;
     reco_track_kin2[k]=-9;
     
     //reco_track_dEdx_arr.at(k).clear();
     //reco_track_ResRan_arr.at(k).clear();
     
     reco_track_length[k]=-9;
     reco_track_theta[k]=-9;
     reco_track_phi[k]=-9;
     
     reco_track_P_vtx[k]=-9;
     reco_track_P_str[k]=-9;
     reco_track_P_end[k]=-9;
     
     reco_track_dir[k]=false;
     reco_track_chi2ka_pl0[k]=-9;
     reco_track_chi2pr_pl0[k]=-9;
     reco_track_chi2pi_pl0[k]=-9;
     reco_track_chi2mu_pl0[k]=-9;
     reco_track_chi2ka_pl1[k]=-9;
     reco_track_chi2pr_pl1[k]=-9;
     reco_track_chi2pi_pl1[k]=-9;
     reco_track_chi2mu_pl1[k]=-9;
     reco_track_chi2ka_pl2[k]=-9;
     reco_track_chi2pr_pl2[k]=-9;
     reco_track_chi2pi_pl2[k]=-9;
     reco_track_chi2mu_pl2[k]=-9;
     
     
     reco_track_Bragg_fwd_ka_pl0[k] = -9; 
     reco_track_Bragg_fwd_pr_pl0[k] = -9; 
     reco_track_Bragg_fwd_pi_pl0[k] = -9; 
     reco_track_Bragg_fwd_mu_pl0[k] = -9; 
     reco_track_Bragg_fwd_ka_pl1[k] = -9;  
     reco_track_Bragg_fwd_pr_pl1[k] = -9;  
     reco_track_Bragg_fwd_pi_pl1[k] = -9; 
     reco_track_Bragg_fwd_mu_pl1[k] = -9; 
     reco_track_Bragg_fwd_ka_pl2[k] = -9; 
     reco_track_Bragg_fwd_pr_pl2[k] = -9; 
     reco_track_Bragg_fwd_pi_pl2[k] = -9; 
     reco_track_Bragg_fwd_mu_pl2[k] = -9; 
     
     reco_track_MIP_pl0[k] = -9; 
     reco_track_MIP_pl1[k] = -9; 
     reco_track_MIP_pl2[k] = -9; 
     
     reco_track_chi2ka_3pl[k]=-99;
     reco_track_chi2pr_3pl[k]=-99;
     reco_track_chi2pi_3pl[k]=-99;
     reco_track_chi2mu_3pl[k]=-99;
     reco_track_likepr_3pl[k]=-99;
     reco_track_llrpid_3pl[k] = -9;
     reco_track_total_llrpid_3pl[k] = -9;
     reco_track_llrpid_k_3pl[k] = -9;
     reco_track_vtx_inTPC[k]=false;
     reco_track_vtx_in5cmTPC[k]=false;
     reco_track_vtx_inCCInclusiveTPC[k]=false;
     reco_track_end_inTPC[k]=false;
     reco_track_end_in5cmTPC[k]=false;
     reco_track_end_inCCInclusiveTPC[k]=false;
     reco_track_true_pdg[k]=-999;
     reco_track_true_origin[k]=-99;
     reco_track_true_primary[k]=false;
     reco_track_true_end_inTPC[k]=false;
     reco_track_true_end_in5cmTPC[k]=false;
     reco_track_true_end_inCCInclusiveTPC[k]=false;
     reco_track_true_length[k]=-9;
     reco_track_ndaughters[k]=0;
     reco_shower_ndaughters[k]=0;
  
     for(int l=0; l<kMaxMerge; l++){ 
       reco_track_match_e[k][l] = -999;
       reco_track_match_hit[k][l] = -999;
       reco_track_match_epdg[k][l] = -999;
       reco_track_match_hitpdg[k][l] = -999;
     }
     
     
       for(int d=0; d<2000; d++){
       reco_track_dEdx_pl0[k][d]=-9;
       reco_track_ResRan_pl0[k][d]=-9;
       reco_track_dEdx_pl1[k][d]=-9;
       reco_track_ResRan_pl1[k][d]=-9;
       reco_track_dEdx_pl2[k][d]=-9;
       reco_track_ResRan_pl2[k][d]=-9;
       }
      
      
     for(int m=0; m<kMaxShowers; m++){ 
       
       reco_track_daughter_distance_sh[k][m]=-9;
       reco_track_daughter_vtx_distance_sh[k][m]=-9;
       reco_angle_track_daughter_sh[k][m]=-9;
       reco_angle_daughter_track_daughter_sh[k][m]=-9;
       
       reco_track_daughter_length_sh[k][m]=-99;
       reco_track_daughter_theta_sh[k][m]=-99;
       reco_track_daughter_phi_sh[k][m]=-99;
       reco_track_daughter_open_angle_sh[k][m]=-99;
       reco_track_daughter_dedx_pl0_sh[k][m]=-99;
       reco_track_daughter_dedx_pl1_sh[k][m]=-99;
       reco_track_daughter_dedx_pl2_sh[k][m]=-99;
       reco_track_daughter_energy_pl0_sh[k][m]=-99;
       reco_track_daughter_energy_pl1_sh[k][m]=-99;
       reco_track_daughter_energy_pl2_sh[k][m]=-99;
       
       reco_track_true_pdg_sh[k]=-999;
       reco_track_true_origin_sh[k]=-99;
       reco_track_true_primary_sh[k]=false;
       reco_track_true_end_inTPC_sh[k]=false;
       reco_track_true_end_in5cmTPC_sh[k]=false;
       reco_track_true_end_inCCInclusiveTPC_sh[k]=false;
       reco_track_true_length_sh[k]=-9;
       
       reco_track_daughter_true_pdg_sh[k][m]=-999;
       reco_track_daughter_true_origin_sh[k][m]=-99;
       reco_track_daughter_true_primary_sh[k][m]=false;
       reco_track_daughter_true_end_inTPC_sh[k][m]=false;
       reco_track_daughter_true_end_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_true_end_inCCInclusiveTPC_sh[k][m]=false;
       reco_track_daughter_true_length_sh[k][m]=-9999;
       reco_track_daughter_true_mother_sh[k][m]=-9999;
       
       reco_track_daughter_vtx_inTPC_sh[k][m]=false;
       reco_track_daughter_vtx_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_vtx_inCCInclusiveTPC_sh[k][m]=false;
       reco_track_daughter_end_inTPC_sh[k][m]=false;
       reco_track_daughter_end_in5cmTPC_sh[k][m]=false;
       reco_track_daughter_end_inCCInclusiveTPC_sh[k][m]=false;
       
       for(int l=0; l<kMaxMerge; l++){ 
	 reco_track_daughter_shower_match_e[k][m][l] = -999;
	 reco_track_daughter_shower_match_hit[k][m][l] = -999;
	 reco_track_daughter_shower_match_epdg[k][m][l] = -999;
	 reco_track_daughter_shower_match_hitpdg[k][m][l] = -999;
       }
     }
   
     for(int m=0; m<kMaxTracks; m++){
       reco_track_daughter_start_x[k][m]=-9;
       reco_track_daughter_start_y[k][m]=-9;
       reco_track_daughter_start_z[k][m]=-9;
       reco_track_daughter_end_x[k][m]=-9;
       reco_track_daughter_end_y[k][m]=-9;
       reco_track_daughter_end_z[k][m]=-9;
       reco_track_daughter_distance[k][m]=-9;
       reco_track_daughter_vtx_distance[k][m]=-9;
       reco_angle_track_daughter[k][m]=-9;
       reco_track_daughter_nhits0[k][m]=-9;
       reco_track_daughter_nhits1[k][m]=-9;
       reco_track_daughter_nhits2[k][m]=-9;
     
       for(int d=0; d<2000; d++){
	 reco_track_daughter_dEdx_pl0[k][m][d]=-9;
	 reco_track_daughter_ResRan_pl0[k][m][d]=-9;
	 reco_track_daughter_dEdx_pl1[k][m][d]=-9;
	 reco_track_daughter_ResRan_pl1[k][m][d]=-9;
	 reco_track_daughter_dEdx_pl2[k][m][d]=-9;
	 reco_track_daughter_ResRan_pl2[k][m][d]=-9;
       }
         
       
       for(int l=0; l<kMaxMerge; l++){ 
	 reco_track_daughter_match_e[k][m][l] = -999;
	 reco_track_daughter_match_hit[k][m][l] = -999;
	 reco_track_daughter_match_epdg[k][m][l] = -999;
	 reco_track_daughter_match_hitpdg[k][m][l] = -999;
       }
       reco_track_daughter_length[k][m]=-99;
       reco_track_daughter_theta[k][m]=-99;
       reco_track_daughter_phi[k][m]=-99;
       reco_track_daughter_chi2ka_pl0[k][m]=-9;
       reco_track_daughter_chi2pr_pl0[k][m]=-9;
       reco_track_daughter_chi2pi_pl0[k][m]=-9;
       reco_track_daughter_chi2mu_pl0[k][m]=-9;
       reco_track_daughter_chi2ka_pl1[k][m]=-9;
       reco_track_daughter_chi2pr_pl1[k][m]=-9;
       reco_track_daughter_chi2pi_pl1[k][m]=-9;
       reco_track_daughter_chi2mu_pl1[k][m]=-9;
       reco_track_daughter_chi2ka_pl2[k][m]=-9;
       reco_track_daughter_chi2pr_pl2[k][m]=-9;
       reco_track_daughter_chi2pi_pl2[k][m]=-9;
       reco_track_daughter_chi2mu_pl2[k][m]=-9;
       reco_track_daughter_chi2ka_3pl[k][m]=-99;
       reco_track_daughter_chi2pr_3pl[k][m]=-99;
       reco_track_daughter_chi2pi_3pl[k][m]=-99;
       reco_track_daughter_chi2mu_3pl[k][m]=-99;
       reco_track_daughter_likepr_3pl[k][m]=-99;
       reco_track_daughter_llrpid_3pl[k][m]=-9;
       reco_track_daughter_llrpid_k_3pl[k][m]=-9;
       reco_track_daughter_vtx_inTPC[k][m]=false;
       reco_track_daughter_vtx_in5cmTPC[k][m]=false;
       reco_track_daughter_vtx_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_end_inTPC[k][m]=false;
       reco_track_daughter_end_in5cmTPC[k][m]=false;
       reco_track_daughter_end_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_true_pdg[k][m]=-999;
       reco_track_daughter_true_origin[k][m]=-99;
       reco_track_daughter_true_primary[k][m]=false;
       reco_track_daughter_true_end_inTPC[k][m]=false;
       reco_track_daughter_true_end_in5cmTPC[k][m]=false;
       reco_track_daughter_true_end_inCCInclusiveTPC[k][m]=false;
       reco_track_daughter_true_length[k][m]=-9999;
       reco_track_daughter_true_mother[k][m]=-9999;
       
       reco_track_daughter_Bragg_fwd_ka_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl0[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl1[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pr_pl1[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pi_pl1[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl1[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl2[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl2[k][m] = -999; 
       
       reco_track_daughter_MIP_pl0[k][m] = -999; 
       reco_track_daughter_MIP_pl1[k][m] = -999; 
       reco_track_daughter_MIP_pl2[k][m] = -999; 
        
        
     }
   }
 
   k_can_trkid=-999;
   mu_can_trkid=-999;
   k_mu_can_dis=-999;
   k_mu_open_angle=-999;
   k_vtx_dis=-999;
   //k_geant_ID=-9999;
   //k_origin=-9999;
   //k_pdg=-9999;
   //k_isPri=-9999;
   //k_endE=-9999;
   //k_ness=-9999;
   //kaon_vtx_dis=-9999;
   //k_plen=-9999;
   //k_phi=-9999;
   //k_theta=-9999;
   //k_in_5_TPC=-9999;
   //k_in_CC_TPC=-9999;
   //k_hit=-9999;
   //k_range=-9999;
   //k_KE=-9999;
   //k_large_dedx=-9999;
   //k_small_dedx=-9999;
   //k_chi_p=-9999;
   //k_chi_k=-9999;
   //k_chi_pi=-9999;
   //k_chi_mu=-9999;
   //k_p_max=-9999;
   //k_k_max=-9999;
   //k_pi_max=-9999;
   //k_mu_max=-9999;
   //k_mip_max=-9999;
   //k_L1_ratio=-9999;
   //k_LL1=-9999;
   //k_L2_ratio=-9999;
   //k_LL2=-9999;
   //k_Lp_ratio=-9999;
   //k_LLp=-9999;
   //k_Lk_ratio=-9999;
   //k_LLk=-9999;
   //k_pida_mean=-9999;
   //k_pida_med=-9999;
   //k_kde=-9999;
   //k_trm_dqdx=-9999;
   //k_trm_dedx=-9999;
   //mu_pdg=-9999;
   //mu_isDec=-9999;
   //mu_origin=-9999;
   //mu_k_is_Mother=-9999;
   //mu_mom_k_inelas=-9999;
   //mu_ness=-9999;
   //mu_plen=-9999;
   //mu_phi=-9999;
   //mu_theta=-9999;
   //mu_in_5_TPC=-9999;
   //mu_in_CC_TPC=-9999;
   //mu_KE=-9999;
   //mu_hit=-9999;
   //mu_range=-9999;
   //mu_large_dedx=-9999;
   //mu_small_dedx=-9999;
   //mu_chi_p=-9999;
   //mu_chi_k=-9999;
   //mu_chi_pi=-9999;
   //mu_chi_mu=-9999;
   //mu_p_max=-9999;
   //mu_k_max=-9999;
   //mu_pi_max=-9999;
   //mu_mu_max=-9999;
   //mu_mip_max=-9999;
   //mu_L1_ratio=-9999;
   //mu_LL1=-9999;
   //mu_L2_ratio=-9999;
   //mu_LL2=-9999;
   //mu_Lp_ratio=-9999;
   //mu_LLp=-9999;
   //mu_Lk_ratio=-9999;
   //mu_LLk=-9999;
   //mu_pida_mean=-9999;
   //mu_pida_med=-9999;
   //mu_kde=-9999;
   //mu_trm_dqdx=-9999;
   //mu_trm_dedx=-9999;
   //mu_mom_process="NA";
   //cc_mu_trkid=-9999;
   //cc_mu_tlen=-9999;
   //cc_mu_phi=-9999;
   //cc_mu_pdg=-9999;
   //cc_mu_theta=-9999;
   //cc_mu_range=-9999;
   //cc_mu_KE=-9999;
   //cc_mu_hit=-9999;
   //cc_mu_large_dedx=-9999;
   //cc_mu_small_dedx=-9999;
   //cc_dis_vtx=-9999;
   //cc_mu_chi_p=-9999;
   //cc_mu_chi_k=-9999;
   //cc_mu_chi_pi=-9999;
   //cc_mu_chi_mu=-9999;
   //cc_mu_p_max=-9999;
   //cc_mu_k_max=-9999;
   //cc_mu_pi_max=-9999;
   //cc_mu_mu_max=-9999;
   //cc_mu_mip_max=-9999;
   //cc_mu_L1_ratio=-9999;
   //cc_mu_LL1=-9999;
   //cc_mu_L2_ratio=-9999;
   //cc_mu_LL2=-9999;
   //cc_mu_Lp_ratio=-9999;
   //cc_mu_LLp=-9999;
   //cc_mu_Lk_ratio=-9999;
   //cc_mu_LLk=-9999;
   //cc_mu_pida_mean=-9999;
   //cc_mu_pida_med=-9999;
   //cc_mu_kde=-9999;
   //cc_mu_trm_dqdx=-9999;
      //cc_mu_trm_dedx=-9999;
      //longest_trkid=-9999;
      //longest_trklen=-9999;
      //pri_Mu_is=-9999;
      //Tr_pri_mu_pdg=-9999;
      //Tr_pri_mu_is=-9999;
      //Tr_pri_st_k_is=-9999;
      //Tr_K_Inelas=-9999;
      //Tr_k_plen=-9999;
      //Tr_k_endE=-9999;
      //Tr_k_theta=-9999;
      //Tr_k_phi=-9999;
      //Tr_dec_mu_is=-9999;
      //Tr_dec_mu_pi_pdg=-9999;
      //Tr_mu_plen=-9999;
      //Tr_mu_theta=-9999;
      //Tr_mu_phi=-9999;
      //Tr_k_inTPC=-9999;
      //Tr_mu_inTPC=-9999;
      //Tr_k_in_5_TPC=-9999;
      //Tr_k_in_CC_TPC=-9999;
      //Tr_mu_in_5_TPC=-9999;
      //Tr_mu_in_CC_TPC=-9999;
      //Tr_kmu_open_ang=-9999;
      //vtx_5cm_mult=-9999;
      //k_start_dedx=-9999;
      //k_end_dedx=-9999;
      //mu_start_dedx=-9999;
      //mu_end_dedx=-9999;
      cut_1=-9;
      cut_2=-9;
      cut_3=-9;
      cut_4=-9;
      cut_5=-9;
      cut_6=-9;
      cut_7=-9;
      cut_8=-9;
      cut_9=-9;
      cut_10=-9;
      cut_11=-9;
      cut_12=-9;
      cut_13=-9;
      //kinelas_has_traks=-9999;
      //kinelas_reco_trkID=-9999;
      //kinelas_tlen=-9999;
      //True_kinelas_KE=-9999;
      //True_kinelas_tlen=-9999;
      
      //for(int j=0; j<3; j++){
      //  for(int k=0; k<3000; k++){
      //    k_dedx[j][k]=-9999;
      //    k_rr[j][k]=-9999;
      //    mu_dedx[j][k]=-9999;
      //    mu_rr[j][k]=-9999;
      //  }
      //}
      //
      //delete [] reco_track_daughter_dEdx_a;
      
      rv0.clear();
      dv0.clear();

      rv1.clear();
      dv1.clear();

      rv2.clear();
      dv2.clear();

 }

 /////////////////////////////////////////////////////////////////////////////////////////
 
bool CCKaonAnalyzerRebuild::isInsideVolume(string volume, double x, double y, double z)
{
  if (volume=="TPC") {
    if (x>0 && x<256.35 && y>-116.5 && y<116.5 && z>0 && z<1036.8) {
      return true;
    }
  }
  else if (volume=="5cmTPC") {
    if (x>5 && x<251.35 && y>-111.5 && y<111.5 && z>5 && z<1031.8) {
      return true;
    }
  }
  else if (volume=="CCInclusiveTPC") {
    if (x>10 && x<246.35 && y>-106.5 && y<106.5 && z>10 && z<986.8) {
      return true;
    }
  }
  return false;
}

double CCKaonAnalyzerRebuild::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  art::ServiceHandle<geo::Geometry> geom;
  double bnd[6] = {0.,2.*geom->DetHalfWidth(),-geom->DetHalfHeight(),geom->DetHalfHeight(),0.,geom->DetLength()};
  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      if(first){
	start = p.Position(i);
	first = false;
	starti = i;
      }else{
	disp -= p.Position(i).Vect();
	result += disp.Mag();
      }
      disp = p.Position(i).Vect();
      end = p.Position(i);
      endi = i;
    }
  }
  return result;
}

//void CCKaonAnalyzerRebuild::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, const recob::Track trk, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart, int track_i, int daughter_i)
//{
void CCKaonAnalyzerRebuild::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i, int daughter_i)
{

  int hits_p0=0;
  int hits_p1=0;
  int hits_p2=0;
  float kin_p0=0;
  float kin_p1=0;
  float kin_p2=0;
    //vector<vector<Float_t>> dvv;
  //vector<vector<Float_t>> rvv;
  //double fADCtoE[]={232,249,243.7};



  double llr_pid_total = 0;
  double llr_pid_total_k = 0;

  //for(unsigned int ical=0; ical<calos.size(); ++ical){
  for (auto const &calo : calos) {

    //if(!calo) continue;
    //if(!calo->PlaneID().isValid) continue;
    int planenum = calo->PlaneID().Plane;
    //if(planenum<0||planenum>2) continue; 
    if(planenum==0) {
            hits_p0=calo->dEdx().size();
            kin_p0=calo->KineticEnergy();
    }

    if(planenum==1) {
            hits_p1=calo->dEdx().size();
            kin_p1=calo->KineticEnergy();
    }

    if(planenum==2) {
            hits_p2=calo->dEdx().size();
            kin_p2=calo->KineticEnergy();
    }

    //auto calo = calos[ical];
    auto const &plane = calo->PlaneID().Plane;
    auto const &dedx_values = calo->dEdx();
    //auto const &dqdx_values = calo->dQdx();

    if(planenum==0){
      for(auto& dedx : calo->dEdx()){
	cout << "calo->dEdx(): " << dedx << endl;
      }
    }


    if(planenum==0) dv0 = calo->dEdx();
    if(planenum==0) rv0 = calo->ResidualRange();

    if(planenum==1) dv1 = calo->dEdx();
    if(planenum==1) rv1 = calo->ResidualRange();

    if(planenum==2) dv2 = calo->dEdx();
    if(planenum==2) rv2 = calo->ResidualRange();

    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();
    //auto const &xyz_v = calo->XYZ();

    std::vector<std::vector<float>> par_values;
    par_values.push_back(rr);
    par_values.push_back(pitch);
/*
    std::vector<float> dqdx_values_corrected;//, dedx_values_corrected;

    if (!isMC){
        dqdx_values_corrected = dqdx_values;
    }
    else{
        dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk, assocMCPart, true, 0.1, false);
    }
    //dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, value, assocMCPart, true, 0.1, false);

    for (size_t i = 0; i < dqdx_values_corrected.size(); i++){
    float aux_dedx;
    aux_dedx = ModBoxCorrection(dqdx_values_corrected[i]*fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
    dedx_values_corrected.push_back(aux_dedx);
    calo_energy += aux_dedx * pitch[i];
    }

    if(planenum==2){ 
            for(auto& n:dqdx_values){cout << "dQdx before: " << n <<endl;}
            for(auto& n:dqdx_values_corrected){cout << "dQdx after: " << n <<endl;}
    }
*/
    if (calo->ResidualRange().size() == 0) continue;

    //llr_pid_total += llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
    //llr_pid_total_k += llr_pid_calculator_k.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
    llr_pid_total += llr_pid_calculator.LLR_many_hits_one_plane(dedx_values, par_values, plane);
    llr_pid_total_k += llr_pid_calculator_k.LLR_many_hits_one_plane(dedx_values, par_values, plane);

  }

  double llr_pid_score = atan(llr_pid_total / 100.) * 2 / 3.14159266;
  double llr_pid_score_k = atan(llr_pid_total_k / 100.) * 2 / 3.14159266;



  if (daughter_i<0) {
    reco_track_nhits0[track_i] = hits_p0;
    reco_track_nhits1[track_i] = hits_p1;
    reco_track_nhits2[track_i] = hits_p2;
    reco_track_kin0[track_i] = kin_p0;
    reco_track_kin1[track_i] = kin_p1;
    reco_track_kin2[track_i] = kin_p2;
    //cout << "Kin p0: " << kin_p0 << " Kin p1: " << kin_p1 << " Kin p2: " << kin_p2 <<endl;
    reco_track_llrpid_3pl[track_i] = llr_pid_score;
    reco_track_total_llrpid_3pl[track_i] = llr_pid_total;
    reco_track_llrpid_k_3pl[track_i] = llr_pid_score_k;

    //reco_track_dEdx.push_back(dv);
    //unsigned int ti = track_i;
    ////reco_track_ResRan.resize(kMaxTracks, std::vector<Float_t>(hits_p2));
    ////reco_track_ResRan[track_i] = rv;
    //reco_track_ResRan.push_back(rv);
    
    //if(dv.size()!=0){
    ////reco_track_dEdx.resize(kMaxTracks, std::vector<Float_t>(hits_p2));
    ////reco_track_dEdx[track_i] = dv;
    //reco_track_ResRan[track_i] = rv;
    //} else{    
    //reco_track_dEdx_arr = dvv;
    //cout << "dEdx has nothing: " << reco_track_dEdx_arr.at(track_i).at(0);
    //}
    
    for(unsigned int d=0; d<dv0.size(); d++){       
       reco_track_dEdx_pl0[track_i][d] = dv0[d];
       reco_track_ResRan_pl0[track_i][d] = rv0[d];
    }
    for(unsigned int d=0; d<dv1.size(); d++){       
       reco_track_dEdx_pl1[track_i][d] = dv1[d];
       reco_track_ResRan_pl1[track_i][d] = rv1[d];
    }
    for(unsigned int d=0; d<dv2.size(); d++){       
       reco_track_dEdx_pl2[track_i][d] = dv2[d];
       reco_track_ResRan_pl2[track_i][d] = rv2[d];
    }


  }
  else {
    reco_track_daughter_nhits0[track_i][daughter_i] = hits_p0;
    reco_track_daughter_nhits1[track_i][daughter_i] = hits_p1;
    reco_track_daughter_nhits2[track_i][daughter_i] = hits_p2;
    reco_track_daughter_llrpid_3pl[track_i][daughter_i] = llr_pid_score;
    reco_track_daughter_llrpid_k_3pl[track_i][daughter_i] = llr_pid_score_k;
  
    cout << "llr_pid_score is " << reco_track_daughter_llrpid_3pl[track_i][daughter_i] << endl;

    //reco_track_daughter_dEdx = vector<vector<vector<Float_t>>>(kMaxTracks,vector<vector<Float_t>>(kMaxTracks,vector<Float_t>(hits_p2,-9999)));

    //dvv.resize(kMaxTracks,vector<Float_t>(hits_p2));
    //reco_track_daughter_dEdx.resize(kMaxTracks, dvv);
    //reco_track_daughter_dEdx[track_i][daughter_i] = dv;

    //rvv.resize(kMaxTracks,vector<Float_t>(hits_p2));
    //reco_track_daughter_ResRan.resize(kMaxTracks, rvv);
    //reco_track_daughter_ResRan[track_i][daughter_i] = rv;

    for(unsigned int d=0; d<dv0.size(); d++){       
       reco_track_daughter_dEdx_pl0[track_i][daughter_i][d] = dv0[d];
       reco_track_daughter_ResRan_pl0[track_i][daughter_i][d] = rv0[d];
    }
    for(unsigned int d=0; d<dv1.size(); d++){       
       reco_track_daughter_dEdx_pl1[track_i][daughter_i][d] = dv1[d];
       reco_track_daughter_ResRan_pl1[track_i][daughter_i][d] = rv1[d];
    }
    for(unsigned int d=0; d<dv2.size(); d++){       
       reco_track_daughter_dEdx_pl2[track_i][daughter_i][d] = dv2[d];
       reco_track_daughter_ResRan_pl2[track_i][daughter_i][d] = rv2[d];
    }


  }

}

/*
float GetLocalEFieldMag(const float x, const float y, const float z){

   const detinfo::DetectorProperties* detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
   auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

   double E_field_nominal = detprop->Efield();        // Electric Field in the drift region in KV/cm

   //correct Efield for SCE
   geo::Vector_t E_field_offsets = {0.,0.,0.};
   E_field_offsets = sce->GetCalEfieldOffsets(geo::Point_t{x,y, z});
   TVector3 E_field_vector = {E_field_nominal*(1 + E_field_offsets.X()), E_field_nominal*E_field_offsets.Y(), E_field_nominal*E_field_offsets.Z()};
   float E_field = E_field_vector.Mag();
   return E_field;
}


double CCKaonAnalyzerRebuild::ModBoxCorrection(const double dQdx, const float x, const float y, const float z) {
              
double rho = 1.383;//detprop->Density();            // LAr density in g/cm^3
double Wion = 23.6/1e6;//util::kGeVToElectrons;  // 23.6 eV = 1e, Wion in MeV/e
                          
auto E_field = GetLocalEFieldMag(x,y,z); // kV / cm
                                
double fModBoxA = 0.930;
double fModBoxB = 0.212;
                                            
double Beta = fModBoxB / (rho * E_field);
double Alpha = fModBoxA;
double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;
                                                              
return dEdx;
                                                                    
}

*/
void CCKaonAnalyzerRebuild::fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i, int daughter_i)
{

  // check PID
  if (trackPID.size()==0) return;

  double chi2ka[3] = {-9.999,-9.999,-9.999};
  double chi2pr[3] = {-9.999,-9.999,-9.999};
  double chi2pi[3] = {-9.999,-9.999,-9.999};
  double chi2mu[3] = {-9.999,-9.999,-9.999};
  Float_t Bragg_fwd_ka[3] = {-1,-1,-1};
  Float_t Bragg_fwd_pr[3] = {-1,-1,-1};
  Float_t Bragg_fwd_pi[3] = {-1,-1,-1};
  Float_t Bragg_fwd_mu[3] = {-1,-1,-1};
  Float_t No_Bragg[3] = {-1,-1,-1};
  double likepr_3pl = -1;

  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
  for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){


    anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
    if (AlgScore.fAlgName == "Chi2") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF) {
        for (int pl=0; pl<3; pl++) {
          if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
            if (AlgScore.fAssumedPdg==321)  chi2ka[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==2212) chi2pr[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==211)  chi2pi[pl]=AlgScore.fValue;
            if (AlgScore.fAssumedPdg==13)   chi2mu[pl]=AlgScore.fValue;
          }
        }
      }
    }
    if (AlgScore.fAlgName == "ThreePlaneProtonPID") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood) {
        if (AlgScore.fPlaneMask==UBPID::uB_SinglePlaneGetBitset(2)) {
          if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward) {
            if (AlgScore.fAssumedPdg==2212) {
              likepr_3pl = AlgScore.fValue;
            }
          }
        }
      }
    }


    if (AlgScore.fAlgName == "BraggPeakLLH") {
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood) {
          //if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward) {
           
              
            for(int pl=0; pl<3; pl++){

            if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
            if (TMath::Abs(AlgScore.fAssumedPdg)==2212) Bragg_fwd_pr[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==321) Bragg_fwd_ka[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==13) Bragg_fwd_mu[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==211) Bragg_fwd_pi[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==0) No_Bragg[pl] = AlgScore.fValue;
            }

            }
          //}
      }
    }

  }

  //double delta_z = end.Z()-pos.Z();
  //double delta_y = end.Y()-pos.Y();
  //double theta0 = TMath::ATan2(delta_z,delta_y) - TMath::Pi()/3;
  //double theta1 = TMath::ATan2(delta_z,delta_y) + TMath::Pi()/3.;
  //double theta2 = TMath::ATan2(delta_z,delta_y);
  double theta0 = angle_y - TMath::Pi()/3;
  double theta1 = angle_y + TMath::Pi()/3.;
  double theta2 = angle_y;
  double wpl0 = TMath::Power(TMath::Sin(theta0),2) >= 0.05 ? 1 : 0;
  double wpl1 = TMath::Power(TMath::Sin(theta1),2) >= 0.05 ? 1 : 0;
  double wpl2 = TMath::Power(TMath::Sin(theta2),2) >= 0.05 ? 1 : 0;
  double chi2ka_3pl = (wpl0*chi2ka[0] + wpl1*chi2ka[1] + wpl2*chi2ka[2])/(wpl0 + wpl1 + wpl2);
  double chi2pr_3pl = (wpl0*chi2pr[0] + wpl1*chi2pr[1] + wpl2*chi2pr[2])/(wpl0 + wpl1 + wpl2);
  double chi2pi_3pl = (wpl0*chi2pi[0] + wpl1*chi2pi[1] + wpl2*chi2pi[2])/(wpl0 + wpl1 + wpl2);
  double chi2mu_3pl = (wpl0*chi2mu[0] + wpl1*chi2mu[1] + wpl2*chi2mu[2])/(wpl0 + wpl1 + wpl2);

  if (daughter_i<0) {
    reco_track_chi2ka_pl0[track_i] = chi2ka[0];
    reco_track_chi2pr_pl0[track_i] = chi2pr[0];
    reco_track_chi2pi_pl0[track_i] = chi2pi[0];
    reco_track_chi2mu_pl0[track_i] = chi2mu[0];
    reco_track_chi2ka_pl1[track_i] = chi2ka[1];
    reco_track_chi2pr_pl1[track_i] = chi2pr[1];
    reco_track_chi2pi_pl1[track_i] = chi2pi[1];
    reco_track_chi2mu_pl1[track_i] = chi2mu[1];
    reco_track_chi2ka_pl2[track_i] = chi2ka[2];
    reco_track_chi2pr_pl2[track_i] = chi2pr[2];
    reco_track_chi2pi_pl2[track_i] = chi2pi[2];
    reco_track_chi2mu_pl2[track_i] = chi2mu[2];


    reco_track_Bragg_fwd_ka_pl0[track_i] = Bragg_fwd_ka[0];
    reco_track_Bragg_fwd_pr_pl0[track_i] = Bragg_fwd_pr[0];
    reco_track_Bragg_fwd_pi_pl0[track_i] = Bragg_fwd_pi[0];
    reco_track_Bragg_fwd_mu_pl0[track_i] = Bragg_fwd_mu[0];
    reco_track_Bragg_fwd_ka_pl1[track_i] = Bragg_fwd_ka[1];
    reco_track_Bragg_fwd_pr_pl1[track_i] = Bragg_fwd_pr[1];
    reco_track_Bragg_fwd_pi_pl1[track_i] = Bragg_fwd_pi[1];
    reco_track_Bragg_fwd_mu_pl1[track_i] = Bragg_fwd_mu[1];
    reco_track_Bragg_fwd_ka_pl2[track_i] = Bragg_fwd_ka[2];
    reco_track_Bragg_fwd_pr_pl2[track_i] = Bragg_fwd_pr[2];
    reco_track_Bragg_fwd_pi_pl2[track_i] = Bragg_fwd_pi[2];
    reco_track_Bragg_fwd_mu_pl2[track_i] = Bragg_fwd_mu[2];

    std::cout << "Track: " << track_i << "\t\tBraggPeakLLH: " << reco_track_Bragg_fwd_pr_pl2[track_i] <<std::endl;

    reco_track_MIP_pl0[track_i] = No_Bragg[0];
    reco_track_MIP_pl1[track_i] = No_Bragg[1];
    reco_track_MIP_pl2[track_i] = No_Bragg[2];

    reco_track_chi2ka_3pl[track_i] = chi2ka_3pl;
    reco_track_chi2pr_3pl[track_i] = chi2pr_3pl;
    reco_track_chi2pi_3pl[track_i] = chi2pi_3pl;
    reco_track_chi2mu_3pl[track_i] = chi2mu_3pl;
    reco_track_likepr_3pl[track_i] = likepr_3pl;
  }
  else {
    reco_track_daughter_chi2ka_pl0[track_i][daughter_i] = chi2ka[0];
    reco_track_daughter_chi2pr_pl0[track_i][daughter_i] = chi2pr[0];
    reco_track_daughter_chi2pi_pl0[track_i][daughter_i] = chi2pi[0];
    reco_track_daughter_chi2mu_pl0[track_i][daughter_i] = chi2mu[0];
    reco_track_daughter_chi2ka_pl1[track_i][daughter_i] = chi2ka[1];
    reco_track_daughter_chi2pr_pl1[track_i][daughter_i] = chi2pr[1];
    reco_track_daughter_chi2pi_pl1[track_i][daughter_i] = chi2pi[1];
    reco_track_daughter_chi2mu_pl1[track_i][daughter_i] = chi2mu[1];
    reco_track_daughter_chi2ka_pl2[track_i][daughter_i] = chi2ka[2];
    reco_track_daughter_chi2pr_pl2[track_i][daughter_i] = chi2pr[2];
    reco_track_daughter_chi2pi_pl2[track_i][daughter_i] = chi2pi[2];
    reco_track_daughter_chi2mu_pl2[track_i][daughter_i] = chi2mu[2];
    reco_track_daughter_chi2ka_3pl[track_i][daughter_i] = chi2ka_3pl;
    reco_track_daughter_chi2pr_3pl[track_i][daughter_i] = chi2pr_3pl;
    reco_track_daughter_chi2pi_3pl[track_i][daughter_i] = chi2pi_3pl;
    reco_track_daughter_chi2mu_3pl[track_i][daughter_i] = chi2mu_3pl;
    reco_track_daughter_likepr_3pl[track_i][daughter_i] = likepr_3pl;

    reco_track_daughter_Bragg_fwd_ka_pl0[track_i][daughter_i] = Bragg_fwd_ka[0];
    reco_track_daughter_Bragg_fwd_pr_pl0[track_i][daughter_i] = Bragg_fwd_pr[0];
    reco_track_daughter_Bragg_fwd_pi_pl0[track_i][daughter_i] = Bragg_fwd_pi[0];
    reco_track_daughter_Bragg_fwd_mu_pl0[track_i][daughter_i] = Bragg_fwd_mu[0];
    reco_track_daughter_Bragg_fwd_ka_pl1[track_i][daughter_i] = Bragg_fwd_ka[1];
    reco_track_daughter_Bragg_fwd_pr_pl1[track_i][daughter_i] = Bragg_fwd_pr[1];
    reco_track_daughter_Bragg_fwd_pi_pl1[track_i][daughter_i] = Bragg_fwd_pi[1];
    reco_track_daughter_Bragg_fwd_mu_pl1[track_i][daughter_i] = Bragg_fwd_mu[1];
    reco_track_daughter_Bragg_fwd_ka_pl2[track_i][daughter_i] = Bragg_fwd_ka[2];
    reco_track_daughter_Bragg_fwd_pr_pl2[track_i][daughter_i] = Bragg_fwd_pr[2];
    reco_track_daughter_Bragg_fwd_pi_pl2[track_i][daughter_i] = Bragg_fwd_pi[2];
    reco_track_daughter_Bragg_fwd_mu_pl2[track_i][daughter_i] = Bragg_fwd_mu[2];

    reco_track_daughter_MIP_pl0[track_i][daughter_i] = No_Bragg[0];
    reco_track_daughter_MIP_pl1[track_i][daughter_i] = No_Bragg[1];
    reco_track_daughter_MIP_pl2[track_i][daughter_i] = No_Bragg[2];



  }

/*
            Bragg_fwd_pr[0] = -9999;
            Bragg_fwd_ka[0] = -9999;
            Bragg_fwd_pi[0] = -9999;
            Bragg_fwd_mu[0] = -9999;
            No_Bragg[0] = -9999;

            Bragg_fwd_pr[1] = -9999;
            Bragg_fwd_ka[1] = -9999;
            Bragg_fwd_pi[1] = -9999;
            Bragg_fwd_mu[1] = -9999;
            No_Bragg[1] = -9999;

            Bragg_fwd_pr[2] = -9999;
            Bragg_fwd_ka[2] = -9999;
            Bragg_fwd_pi[2] = -9999;
            Bragg_fwd_mu[2] = -9999;
            No_Bragg[2] = -9999;

*/


}

  simb::MCParticle const* CCKaonAnalyzerRebuild::truthMatchTrack(std::vector<art::Ptr<recob::Hit>>& hits_from_track, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit){

    simb::MCParticle const* matched_mcparticle = NULL;

    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

    for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
      particle_vec.clear(); match_vec.clear();
      particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

      //cout << "track's particle_vec.size(): " << particle_vec.size() << endl;                                                                                                                               

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
        trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
        //trkide[ particle_vec[i_p]->TrackId() ] ++;                                                                                                                                                          
        tote += match_vec[i_p]->energy;
        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
          maxe = trkide[ particle_vec[i_p]->TrackId() ];
          matched_mcparticle = particle_vec[i_p];
        }
      }
    }
    if(matched_mcparticle) return matched_mcparticle;
    else return NULL;
  }

  simb::MCParticle const* CCKaonAnalyzerRebuild::truthMatchHit(art::Ptr<recob::Hit>& hit, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit){

    simb::MCParticle const* matched_mcparticle = NULL;

    double maxe=-1;
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(), particle_vec, match_vec);

    //cout << particle_vec.size() << endl;
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      //cout<< "match_vec[i_p]->energy: " << match_vec[i_p]->energy << endl;

      if(match_vec[i_p]->energy>maxe){
	maxe = match_vec[i_p]->energy;
	matched_mcparticle = particle_vec[i_p];
      }
    }
    if(matched_mcparticle) return matched_mcparticle;
    else return NULL;
  }

void CCKaonAnalyzerRebuild::fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{

  //std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
  //art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);

  //cout << "hits_from_track.size(): " << hits_from_track.size() << endl;

  simb::MCParticle const* matched_mcparticle = NULL;
  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

    //cout << "track's particle_vec.size(): " << particle_vec.size() << endl;

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      //trkide[ particle_vec[i_p]->TrackId() ] ++;
      tote += match_vec[i_p]->energy;
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
        maxe = trkide[ particle_vec[i_p]->TrackId() ];
        matched_mcparticle = particle_vec[i_p];
      }
    }
  }
 

  if(matched_mcparticle){
    //const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;

    if (daughter_i<0) {
      cout << "this case track_i is: " << track_i << endl;
      cout << "matched_mcparticle->PdgCode(): " << matched_mcparticle->PdgCode() << endl;
      reco_track_true_pdg[track_i] = matched_mcparticle->PdgCode();
      reco_track_true_origin[track_i] = 1;//int(mc_truth->Origin());
      reco_track_true_primary[track_i] = matched_mcparticle->Process()=="primary";
      reco_track_true_end_inTPC[track_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_end_in5cmTPC[track_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_end_inCCInclusiveTPC[track_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_true_length[track_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
    }
    else {
      reco_track_daughter_true_pdg[track_i][daughter_i] = matched_mcparticle->PdgCode();
      reco_track_daughter_true_origin[track_i][daughter_i] = 1;//int(mc_truth->Origin());
      reco_track_daughter_true_primary[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
      reco_track_daughter_true_end_inTPC[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_in5cmTPC[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_inCCInclusiveTPC[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_length[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      //for (auto const& pPart : ptList) {
      //  if (pPart->TrackId()==matched_mcparticle2->Mother()) {
      //    reco_track_daughter_true_mother[ntracks][ndaughters] = pPart->PdgCode();
      //    break;
      //  }
      //}
    }
  }else{
      cout << "no matched_mcparticle found in fillTrue" << endl;
  }

}


void CCKaonAnalyzerRebuild::fillTrackMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{
  std::map<int,double> trkide;
  std::map<int,double> trkidhit;
  std::map<int,int> trkidpdg;    
  std::map<int,int> trkidmrgid;
  std::map<int,int> trkidmother;

  map<int, int>::iterator it_trkidpdg;
  map<int, int>::iterator it_trkidmrgid;
  map<int, int>::iterator it_trkidmrgid2;
  map<int, double>::iterator it_trkide;
  map<int, double>::iterator it_trkidhit;

  std::map<double, int, std::greater<int>> epdg;
  std::map<double, int, std::greater<int>> hitpdg;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_track.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();

    particles_per_hit.get(hits_from_track[i_h].key(), particle_vec, match_vec);

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      trkidhit[ particle_vec[i_p]->TrackId() ]++;
      trkidpdg[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->PdgCode();
      trkidmrgid[ particle_vec[i_p]->TrackId() ] = -9;
      trkidmother[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->Mother();
    }
  }


  //get merged id for two generations

  int currentMergedId = 1; 

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
        
    if(trkidmrgid[it_trkidpdg->first] != -9) continue;  //not filled yet
    trkidmrgid[it_trkidpdg->first] = currentMergedId; //fill current id
    int currentMotherTrackId = trkidmother[it_trkidpdg->first];
    
    while (currentMotherTrackId > 0) {//have mother

      if(trkidpdg.find(currentMotherTrackId)==trkidpdg.end()) break;
      if(trkidpdg[currentMotherTrackId]!=it_trkidpdg->second) break;
      
	trkidmrgid[currentMotherTrackId] = currentMergedId;
	currentMotherTrackId = trkidmother[currentMotherTrackId];
      }
    ++currentMergedId; 
  }                                                                                             
 
  
  for (it_trkidmrgid = trkidmrgid.begin(); it_trkidmrgid != trkidmrgid.end(); it_trkidmrgid++){
    for (it_trkidmrgid2 = trkidmrgid.begin(); it_trkidmrgid2 != trkidmrgid.end(); it_trkidmrgid2++){
      
      if(it_trkidmrgid->first == it_trkidmrgid2->first) continue; //skip if they have the same trackID
      
      if(it_trkidmrgid->second == it_trkidmrgid2->second){ //check if they have same MergedId 
	//trkide.find( it_trkidmrgid->first );
	trkide[ it_trkidmrgid->first ] += trkide[ it_trkidmrgid2->first ];
	trkidhit[ it_trkidmrgid->first ] += trkidhit[ it_trkidmrgid2->first ];
	
	trkidmrgid[ it_trkidmrgid2->first ] = -1; //removed if it's been merged 
	trkide[ it_trkidmrgid2->first ] = -1;
	trkidhit[ it_trkidmrgid2->first ] = -1;
      }      
    }	
  }
  

  it_trkide = trkide.begin();
  it_trkidhit = trkidhit.begin();
  it_trkidmrgid = trkidmrgid.begin();

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){

    if(it_trkidmrgid->second ==-1) continue;
    epdg[ it_trkide->second ] = it_trkidpdg->second;
    hitpdg[ it_trkidhit->second ] = it_trkidpdg->second;

    it_trkide++;//? how to initialise
    it_trkidhit++;
    it_trkidmrgid++;
  }


  Int_t merge_i = 0;  
 
  /* 
cout << "TRKID AND MRGID" << endl;

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
    if(trkidpdg.find( trkidmother[ it_trkidpdg->first] ) == trkidpdg.end()) continue;
    cout << "trkid: " << it_trkidpdg->first << ", mrgid: " << trkidmrgid[ it_trkidpdg->first] << ", PDG: " << it_trkidpdg->second << endl;
    cout << "trkid: " << it_trkidpdg->first << ", mother: " << trkidmother[ it_trkidpdg->first] << ", PDG: " << trkidpdg[ trkidmother[ it_trkidpdg->first] ] << endl;
  }

  cout << "ENERGY AND ITS PDG" << endl;
  */

  for (auto const& x : epdg){

    if(daughter_i<0){
      reco_track_match_e[track_i][merge_i] = x.first;
      //cout << "after fill reco_track_match_e[track_i][merge_i]: " << reco_track_match_e[track_i][merge_i] << endl;
      reco_track_match_epdg[track_i][merge_i] = x.second;

    }else{
      reco_track_daughter_match_e[track_i][daughter_i][merge_i] = x.first;
      reco_track_daughter_match_epdg[track_i][daughter_i][merge_i] = x.second;
    }
    merge_i++;   
  }
  
  merge_i=0;
  //cout << "HITS AND ITS PDG" << endl;
  for (auto const& y : hitpdg){
    if(daughter_i<0){ 
      reco_track_match_hit[track_i][merge_i] = y.first;     
      reco_track_match_hitpdg[track_i][merge_i] = y.second;
    }else{
      reco_track_daughter_match_e[track_i][daughter_i][merge_i] = y.first;
      reco_track_daughter_match_epdg[track_i][daughter_i][merge_i] = y.second;
    }
    merge_i++;
  }
  
}


void CCKaonAnalyzerRebuild::fillTrueMatching_sh(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
                                         art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                         int track_i,
                                         int daughter_i)
{

  simb::MCParticle const* matched_mcparticle = NULL;
  std::unordered_map<int,double> shwride;
  double maxe=-1, tote=0;
  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hits_from_shower[i_h].key(), particle_vec, match_vec);

    //cout << "shower's particle_vec.size(): " << particle_vec.size() << endl;

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      shwride[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      tote += match_vec[i_p]->energy;

      //cout << "particle_vec[i_p]->TrackId(): " << particle_vec[i_p]->TrackId() << endl;
      //cout << "match_vec[i_p]->energy: " << match_vec[i_p]->energy << endl;
      //cout << "shwride[ particle_vec[i_p]->TrackId() ]: " << shwride[ particle_vec[i_p]->TrackId() ] << endl;
      if( shwride[ particle_vec[i_p]->TrackId() ] > maxe ){
        maxe = shwride[ particle_vec[i_p]->TrackId() ];
        matched_mcparticle = particle_vec[i_p];
      }
    }
  }
  //cout << "MAXE: " << maxe << endl;


  if(matched_mcparticle){
    //const art::Ptr<simb::MCTruth> mc_truth=TrackIDToMCTruth(evt,"largeant",matched_mcparticle->TrackId());
    TLorentzVector mcstart, mcend;
    unsigned int pstarti, pendi;

    //cout << "found matched_mcparticle" << endl;
    //cout << "and its pdg: " << matched_mcparticle->PdgCode() << endl;

    if(track_i>=0 && daughter_i>=0){

      //cout << "filling reco_track_daughter_true_pdg_sh[track_i][daughter_i] with track_i: " << track_i << ", daughter_i: " << daughter_i << ", pdg" << matched_mcparticle->PdgCode() << endl;
      //cout << "before filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[track_i][daughter_i] << endl;
      //cout << "after filling reco_track_daughter_true_pdg_sh[track_i][daughter_i]: " << reco_track_daughter_true_pdg_sh[track_i][daughter_i] << endl;

      reco_track_daughter_true_pdg_sh[track_i][daughter_i] = matched_mcparticle->PdgCode();

      reco_track_daughter_true_origin_sh[track_i][daughter_i] = 1;//int(mc_truth->Origin());
      reco_track_daughter_true_primary_sh[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
      reco_track_daughter_true_end_inTPC_sh[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_in5cmTPC_sh[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_end_inCCInclusiveTPC_sh[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
      reco_track_daughter_true_length_sh[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      //}
    }
  }
}


void CCKaonAnalyzerRebuild::fillShowerMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i)
{
  std::map<int,double> trkide;
  std::map<int,double> trkidhit;
  std::map<int,int> trkidpdg;
  std::map<int,int> trkidmrgid;
  std::map<int,int> trkidmother;

  map<int, int>::iterator it_trkidpdg;
  map<int, int>::iterator it_trkidmrgid;
  map<int, int>::iterator it_trkidmrgid2;
  map<int, double>::iterator it_trkide;
  map<int, double>::iterator it_trkidhit;

  std::map<double, int, std::greater<int>> epdg;
  std::map<double, int, std::greater<int>> hitpdg;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_shower.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();

    particles_per_hit.get(hits_from_shower[i_h].key(), particle_vec, match_vec);

    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p) {
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
      trkidhit[ particle_vec[i_p]->TrackId() ]++;
      trkidpdg[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->PdgCode();
      trkidmrgid[ particle_vec[i_p]->TrackId() ] = -9;
      trkidmother[ particle_vec[i_p]->TrackId() ] = particle_vec[i_p]->Mother();
    }
  }
 
  //get merged id for two generations

  int currentMergedId = 1; 

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
        
    if(trkidmrgid[it_trkidpdg->first] != -9) continue;  //not filled yet
    trkidmrgid[it_trkidpdg->first] = currentMergedId; //fill current id
    int currentMotherTrackId = trkidmother[it_trkidpdg->first];
    
    while (currentMotherTrackId > 0) {//have mother

      if(trkidpdg.find(currentMotherTrackId)==trkidpdg.end()) break;
      if(trkidpdg[currentMotherTrackId]!=it_trkidpdg->second) break;
      
	trkidmrgid[currentMotherTrackId] = currentMergedId;
	currentMotherTrackId = trkidmother[currentMotherTrackId];
      }
    ++currentMergedId; 
  }                                                                                             
 
  
  for (it_trkidmrgid = trkidmrgid.begin(); it_trkidmrgid != trkidmrgid.end(); it_trkidmrgid++){
    for (it_trkidmrgid2 = trkidmrgid.begin(); it_trkidmrgid2 != trkidmrgid.end(); it_trkidmrgid2++){
      
      if(it_trkidmrgid->first == it_trkidmrgid2->first) continue; //skip if they have the same trackID
      
      if(it_trkidmrgid->second == it_trkidmrgid2->second){ //check if they have same MergedId 
	trkide[ it_trkidmrgid->first ] += trkide[ it_trkidmrgid2->first ];
	trkidhit[ it_trkidmrgid->first ] += trkidhit[ it_trkidmrgid2->first ];
	
	trkidmrgid[ it_trkidmrgid2->first ] = -1; //removed if it's been merged 
	trkide[ it_trkidmrgid2->first ] = -1;
	trkidhit[ it_trkidmrgid2->first ] = -1;
      }      
    }	
  }
  

  it_trkide = trkide.begin();
  it_trkidhit = trkidhit.begin();
  it_trkidmrgid = trkidmrgid.begin();

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){

    if(it_trkidmrgid->second ==-1) continue;
    epdg[ it_trkide->second ] = it_trkidpdg->second;
    hitpdg[ it_trkidhit->second ] = it_trkidpdg->second;

    it_trkide++;
    it_trkidhit++;
    it_trkidmrgid++;
  }

  Int_t merge_i = 0; 

  /*
  cout << "TRKID AND MRGID" << endl;

  for (it_trkidpdg = trkidpdg.begin(); it_trkidpdg != trkidpdg.end(); it_trkidpdg++){
    cout << "shwid: " << it_trkidpdg->first << ", mrgid: " << trkidmrgid[ it_trkidpdg->first] << ", PDG: " << it_trkidpdg->second << endl;
    if(trkidpdg.find( trkidmother[ it_trkidpdg->first] ) == trkidpdg.end()) continue;
    cout << "shwid: " << it_trkidpdg->first << ", mother: " << trkidmother[ it_trkidpdg->first] << ", PDG: " << trkidpdg[ trkidmother[ it_trkidpdg->first] ] << endl;
  }

  cout << "ENERGY AND ITS PDG" << endl;
  */

  for (auto const& x : epdg){

    if(track_i>=0 && daughter_i>=0){
      //cout << "track_i: " << track_i << ", daughter_i: " << daughter_i << ", merge_i: " << merge_i << endl;
      //cout << "before fill reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i]: " << reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] << endl;
      reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] = x.first;
      //cout << "after fill reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i]: " << reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] << endl;

      reco_track_daughter_shower_match_epdg[track_i][daughter_i][merge_i] = x.second;
    }
    merge_i++;   
  }
  
  merge_i=0;

  //cout << "HITS AND ITS PDG" << endl;

  for (auto const& y : hitpdg){
    if(track_i>=0 && daughter_i>=0){ 
      reco_track_daughter_shower_match_e[track_i][daughter_i][merge_i] = y.first;
      reco_track_daughter_shower_match_epdg[track_i][daughter_i][merge_i] = y.second;
    }
    merge_i++;
  }
  
}



///////////////////////////////////////////////////////////////////////////////////////////

//========================================================================



  //DEFINE_ART_MODULE(CCKaonAnalyzerRebuild)
}
//}
