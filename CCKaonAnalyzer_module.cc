#include "CCKaonAnalyzer_module.h"

//#ifdef __MAKECINT__
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+; 
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+; 
#endif

using namespace std;
namespace Kaon_Analyzer{

	
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
CCKaonAnalyzer::CCKaonAnalyzer(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitsModuleLabel                 (pset.get< std::string >("HitsModuleLabel","gaushit")),
  fLArG4ModuleLabel                (pset.get< std::string >("LArG4ModuleLabel","largeant")),
  fGenieGenModuleLabel             (pset.get< std::string >("GenieGenModuleLabel","generator")),  
  fTrackModuleLabel                (pset.get< std::string >("TrackModuleLabel","pandora")),
  fShowerModuleLabel               (pset.get< std::string >("ShowerModuleLabel","pandora")),
  fCalorimetryModuleLabel          (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),
  fPIDLabel                        (pset.get< std::string >("PIDLabel","pandoracalipidSCE")),
  fHitTruthAssns                   (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")), 
  fHitTrackAssns                   (pset.get< std::string >("HitTrackAssn","pandora")), 
  fHitShowerAssns                  (pset.get< std::string >("HitShowerAssn","pandora")), 
  m_pfp_producer                   (pset.get< std::string >("pfp_producer","pandora")),
  fPFParticleLabel                 (pset.get< std::string >("PFParticleLabel", "pandora")),
  fSpacePointproducer              (pset.get< std::string >("SpacePointproducer", "pandora")),
  isMC                             (pset.get< bool >("IsMC",true)),
  fRebuiltTrackModuleLabel         (pset.get< std::string >("RebuiltTrackModuleLabel","CCKaonProducer")),
  fRebuiltCalorimetryModuleLabel   (pset.get< std::string >("RebuiltCalorimetryModuleLabel","pandoraRebuildTrackCaliSCE")),
  fRebuiltPIDLabel                 (pset.get< std::string >("RebuiltPIDLabel","pandoraRebuildTrackCaliPidSCE")),
  fRebuiltHitTrackAssns            (pset.get< std::string >("RebuiltHitTrackAssn","CCKaonProducer"))

{

  m_is_verbose = pset.get<bool>("Verbose",false);
  m_pandoraLabel = pset.get<std::string>("PandoraLabel"); 


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
CCKaonAnalyzer::~CCKaonAnalyzer()
{
  //destructor
}

//========================================================================
void CCKaonAnalyzer::endSubRun(const art::SubRun &subrun)
{
  //if (!m_isData)
  //{
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    m_pot = subrun.getByLabel("generator", potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    // -- std::cout << "[CCKaonAnalyzer::endSubRun] Storing POT info!" << std::endl;
  //}

  m_run = subrun.run();
  m_subrun = subrun.subRun();
  fSubrunTree->Fill();
}

//========================================================================
void CCKaonAnalyzer::beginJob()
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

  fEventTree->Branch("cheat_peak_pdg", &cheat_peak_pdg, "cheat_peak_pdg[20][10]/F");
  fEventTree->Branch("cheat_peak_theta", &cheat_peak_theta, "cheat_peak_theta[20][10]/F");
  fEventTree->Branch("cheat_pip_trkln", &cheat_pip_trkln, "cheat_pip_trkln/F");
  fEventTree->Branch("cheat_mup_trkln", &cheat_mup_trkln, "cheat_mup_trkln/F");

  fEventTree->Branch("cheat_peak_phi", &cheat_peak_phi, "cheat_peak_phi[20][10]/F");
  fEventTree->Branch("best_peak_theta", &best_peak_theta, "best_peak_theta[20][10]/F");
  fEventTree->Branch("best_peak_phi", &best_peak_phi, "best_peak_phi[20][10]/F");
  fEventTree->Branch("best_peak_trkln", &best_peak_trkln, "best_peak_trkln[20][10]/F");

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

  fEventTree->Branch("true_nkaons_anti", &true_nkaons_anti, "true_nkaons_anti/I");
  fEventTree->Branch("true_kaon_length_anti", &true_kaon_length_anti, "true_kaon_length_anti/F");
  fEventTree->Branch("true_kaon_p_anti", &true_kaon_p_anti, "true_kaon_p_anti/F");
  fEventTree->Branch("true_kaon_ke_anti", &true_kaon_ke_anti, "true_kaon_ke_anti/F");
  fEventTree->Branch("true_kaon_theta_anti", &true_kaon_theta_anti, "true_kaon_theta_anti/F");
  fEventTree->Branch("true_kaon_costheta_anti", &true_kaon_costheta_anti, "true_kaon_costheta_anti/F");
  fEventTree->Branch("true_kaon_phi_anti", &true_kaon_phi_anti, "true_kaon_phi_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_angle_anti", &true_kaon_ccmuon_angle_anti, "true_kaon_ccmuon_angle_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_cosangle_anti", &true_kaon_ccmuon_cosangle_anti, "true_kaon_ccmuon_cosangle_anti/F");
  fEventTree->Branch("true_kaon_end_process_anti", &true_kaon_end_process_anti, "true_kaon_end_process_anti/I");
  fEventTree->Branch("true_kaon_end_ke_anti", &true_kaon_end_ke_anti, "true_kaon_end_ke_anti/F");
  fEventTree->Branch("true_kaon_end_x_anti", &true_kaon_end_x_anti, "true_kaon_end_x_anti/F");
  fEventTree->Branch("true_kaon_end_y_anti", &true_kaon_end_y_anti, "true_kaon_end_y_anti/F");
  fEventTree->Branch("true_kaon_end_z_anti", &true_kaon_end_z_anti, "true_kaon_end_z_anti/F");
  fEventTree->Branch("true_kaon_end_inTPC_anti", &true_kaon_end_inTPC_anti, "true_kaon_end_inTPC_anti/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC_anti", &true_kaon_end_in5cmTPC_anti, "true_kaon_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC_anti", &true_kaon_end_inCCInclusiveTPC_anti, "true_kaon_end_inCCInclusiveTPC_anti/O");

  fEventTree->Branch("true_kaon_ndaughters_anti", &true_kaon_ndaughters_anti, "true_kaon_ndaughters_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_decay_anti", &true_kaon_ndaughters_decay_anti, "true_kaon_ndaughters_decay_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_inelastic_anti", &true_kaon_ndaughters_inelastic_anti, "true_kaon_ndaughters_inelastic_anti/I");
  fEventTree->Branch("true_kaon_ndecmup_anti", &true_kaon_ndecmup_anti, "true_kaon_ndecmup_anti/I");
  fEventTree->Branch("true_kaon_ndecpip_anti", &true_kaon_ndecpip_anti, "true_kaon_ndecpip_anti/I");
  fEventTree->Branch("true_kaon_ninekap_anti", &true_kaon_ninekap_anti, "true_kaon_ninekap_anti/I");
  fEventTree->Branch("true_kaon_ninepip_anti", &true_kaon_ninepip_anti, "true_kaon_ninepip_anti/I");
  fEventTree->Branch("true_kaon_ninepro_anti", &true_kaon_ninepro_anti, "true_kaon_ninepro_anti/I");
  fEventTree->Branch("true_kaon_daughter_length_anti", &true_kaon_daughter_length_anti, "true_kaon_daughter_length_anti/F");
  fEventTree->Branch("true_kaon_daughter_p_anti", &true_kaon_daughter_p_anti, "true_kaon_daughter_p_anti/F");
  fEventTree->Branch("true_kaon_daughter_ke_anti", &true_kaon_daughter_ke_anti, "true_kaon_daughter_ke_anti/F");
  fEventTree->Branch("true_kaon_daughter_theta_anti", &true_kaon_daughter_theta_anti, "true_kaon_daughter_theta_anti/F");
  fEventTree->Branch("true_kaon_daughter_costheta_anti", &true_kaon_daughter_costheta_anti, "true_kaon_daughter_costheta_anti/F");
  fEventTree->Branch("true_kaon_daughter_angle_anti", &true_kaon_daughter_angle_anti, "true_kaon_daughter_angle_anti/F");
  fEventTree->Branch("true_kaon_daughter_cosangle_anti", &true_kaon_daughter_cosangle_anti, "true_kaon_daughter_cosangle_anti/F");
  fEventTree->Branch("true_kaon_daughter_pdg_anti", &true_kaon_daughter_pdg_anti, "true_kaon_daughter_pdg_anti/I");
  fEventTree->Branch("true_kaon_daughter_end_x_anti", &true_kaon_daughter_end_x_anti, "true_kaon_daughter_end_x_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_y_anti", &true_kaon_daughter_end_y_anti, "true_kaon_daughter_end_y_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_z_anti", &true_kaon_daughter_end_z_anti, "true_kaon_daughter_end_z_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC_anti", &true_kaon_daughter_end_inTPC_anti, "true_kaon_daughter_endInTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC_anti", &true_kaon_daughter_end_in5cmTPC_anti, "true_kaon_daughter_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC_anti", &true_kaon_daughter_end_inCCInclusiveTPC_anti, "true_kaon_daughter_end_inCCInclusiveTPC_anti/O");


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

  fEventTree->Branch("reco_track_start_x", &reco_track_start_x, "reco_track_start_x[20]/F");
  fEventTree->Branch("reco_track_start_y", &reco_track_start_y, "reco_track_start_y[20]/F");
  fEventTree->Branch("reco_track_start_z", &reco_track_start_z, "reco_track_start_z[20]/F");
  fEventTree->Branch("reco_track_end_x", &reco_track_end_x, "reco_track_end_x[20]/F");
  fEventTree->Branch("reco_track_end_y", &reco_track_end_y, "reco_track_end_y[20]/F");
  fEventTree->Branch("reco_track_end_z", &reco_track_end_z, "reco_track_end_z[20]/F");


  fEventTree->Branch("reco_track_distance", &reco_track_distance, "reco_track_distance[20]/F");

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



  fEventTree->Branch("reco_track_chi2ka_3pl", &reco_track_chi2ka_3pl, "reco_track_chi2ka_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pr_3pl", &reco_track_chi2pr_3pl, "reco_track_chi2pr_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pi_3pl", &reco_track_chi2pi_3pl, "reco_track_chi2pi_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2mu_3pl", &reco_track_chi2mu_3pl, "reco_track_chi2mu_3pl[20]/F");
  fEventTree->Branch("reco_track_mean_dedx_3pl", &reco_track_mean_dedx_3pl, "reco_track_mean_dedx_3pl[20]/F");
  fEventTree->Branch("reco_track_mean_dedx_pl0", &reco_track_mean_dedx_pl0, "reco_track_mean_dedx_pl0[20]/F");
  fEventTree->Branch("reco_track_mean_dedx_pl1", &reco_track_mean_dedx_pl1, "reco_track_mean_dedx_pl1[20]/F");
  fEventTree->Branch("reco_track_mean_dedx_pl2", &reco_track_mean_dedx_pl2, "reco_track_mean_dedx_pl2[20]/F");
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

  fEventTree->Branch("reco_track_daughter_match_e_rebuild", &reco_track_daughter_match_e_rebuild, "reco_track_daughter_match_e_rebuild[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_hit_rebuild", &reco_track_daughter_match_hit_rebuild, "reco_track_daughter_match_hit_rebuild[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_epdg_rebuild", &reco_track_daughter_match_epdg_rebuild, "reco_track_daughter_match_epdg_rebuild[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_match_hitpdg_rebuild", &reco_track_daughter_match_hitpdg_rebuild, "reco_track_daughter_match_hitpdg_rebuild[20][20][10]/I");

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
  fEventTree->Branch("reco_track_daughter_vtx_distance_rebuild", &reco_track_daughter_vtx_distance_rebuild, "reco_track_daughter_vtx_distance_rebuild[20][20]/F");
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
  fEventTree->Branch("reco_track_daughter_mean_dedx_3pl", &reco_track_daughter_mean_dedx_3pl, "reco_track_daughter_mean_dedx_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl0", &reco_track_daughter_mean_dedx_pl0, "reco_track_daughter_mean_dedx_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl1", &reco_track_daughter_mean_dedx_pl1, "reco_track_daughter_mean_dedx_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl2", &reco_track_daughter_mean_dedx_pl2, "reco_track_daughter_mean_dedx_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_3pl", &reco_track_daughter_llrpid_3pl, "reco_track_daughter_llrpid_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_k_3pl", &reco_track_daughter_llrpid_k_3pl, "reco_track_daughter_llrpid_k_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC", &reco_track_daughter_vtx_inTPC, "reco_track_daughter_vtx_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC", &reco_track_daughter_vtx_in5cmTPC, "reco_track_daughter_vtx_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC", &reco_track_daughter_vtx_inCCInclusiveTPC, "reco_track_daughter_vtx_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC", &reco_track_daughter_end_inTPC, "reco_track_daughter_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC", &reco_track_daughter_end_in5cmTPC, "reco_track_daughter_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC", &reco_track_daughter_end_inCCInclusiveTPC, "reco_track_daughter_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_pdg_rebuild", &reco_track_daughter_true_pdg_rebuild, "reco_track_daughter_true_pdg_rebuild[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin_rebuild", &reco_track_daughter_true_origin_rebuild, "reco_track_daughter_true_origin_rebuild[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary_rebuild", &reco_track_daughter_true_primary_rebuild, "reco_track_daughter_true_primary_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC_rebuild", &reco_track_daughter_true_end_inTPC_rebuild, "reco_track_daughter_true_end_inTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC_rebuild", &reco_track_daughter_true_end_in5cmTPC_rebuild, "reco_track_daughter_true_end_in5cmTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC_rebuild", &reco_track_daughter_true_end_inCCInclusiveTPC_rebuild, "reco_track_daughter_true_end_inCCInclusiveTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length_rebuild", &reco_track_daughter_true_length_rebuild, "reco_track_daughter_true_length_rebuild[20][20]/F");

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


  fEventTree->Branch("reco_track_daughter_start_x_rebuild", &reco_track_daughter_start_x_rebuild, "reco_track_daughter_start_x_rebuild[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_y_rebuild", &reco_track_daughter_start_y_rebuild, "reco_track_daughter_start_y_rebuild[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_start_z_rebuild", &reco_track_daughter_start_z_rebuild, "reco_track_daughter_start_z_rebuild[20][20]/F"); 
  fEventTree->Branch("reco_track_daughter_end_x_rebuild", &reco_track_daughter_end_x_rebuild, "reco_track_daughter_end_x_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_y_rebuild", &reco_track_daughter_end_y_rebuild, "reco_track_daughter_end_y_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_z_rebuild", &reco_track_daughter_end_z_rebuild, "reco_track_daughter_end_z_rebuild[20][20]/F");


  fEventTree->Branch("reco_track_ndaughters_rebuild", &reco_track_ndaughters_rebuild, "reco_track_ndaughters_rebuild[20]/I");
  fEventTree->Branch("reco_track_daughter_distance_rebuild", &reco_track_daughter_distance_rebuild, "reco_track_daughter_distance_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_length_rebuild", &reco_track_daughter_length_rebuild, "reco_track_daughter_length_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta_rebuild", &reco_track_daughter_theta_rebuild, "reco_track_daughter_theta_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi_rebuild", &reco_track_daughter_phi_rebuild, "reco_track_daughter_phi_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl0_rebuild", &reco_track_daughter_chi2ka_pl0_rebuild, "reco_track_daughter_chi2ka_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl0_rebuild", &reco_track_daughter_chi2pr_pl0_rebuild, "reco_track_daughter_chi2pr_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl0_rebuild", &reco_track_daughter_chi2pi_pl0_rebuild, "reco_track_daughter_chi2pi_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl0_rebuild", &reco_track_daughter_chi2mu_pl0_rebuild, "reco_track_daughter_chi2mu_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl1_rebuild", &reco_track_daughter_chi2ka_pl1_rebuild, "reco_track_daughter_chi2ka_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl1_rebuild", &reco_track_daughter_chi2pr_pl1_rebuild, "reco_track_daughter_chi2pr_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl1_rebuild", &reco_track_daughter_chi2pi_pl1_rebuild, "reco_track_daughter_chi2pi_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl1_rebuild", &reco_track_daughter_chi2mu_pl1_rebuild, "reco_track_daughter_chi2mu_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl2_rebuild", &reco_track_daughter_chi2ka_pl2_rebuild, "reco_track_daughter_chi2ka_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl2_rebuild", &reco_track_daughter_chi2pr_pl2_rebuild, "reco_track_daughter_chi2pr_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl2_rebuild", &reco_track_daughter_chi2pi_pl2_rebuild, "reco_track_daughter_chi2pi_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl2_rebuild", &reco_track_daughter_chi2mu_pl2_rebuild, "reco_track_daughter_chi2mu_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_3pl_rebuild", &reco_track_daughter_chi2ka_3pl_rebuild, "reco_track_daughter_chi2ka_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_3pl_rebuild", &reco_track_daughter_chi2pr_3pl_rebuild, "reco_track_daughter_chi2pr_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_3pl_rebuild", &reco_track_daughter_chi2pi_3pl_rebuild, "reco_track_daughter_chi2pi_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_3pl_rebuild", &reco_track_daughter_chi2mu_3pl_rebuild, "reco_track_daughter_chi2mu_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_3pl_rebuild", &reco_track_daughter_mean_dedx_3pl_rebuild, "reco_track_daughter_mean_dedx_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl0_rebuild", &reco_track_daughter_mean_dedx_pl0_rebuild, "reco_track_daughter_mean_dedx_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl1_rebuild", &reco_track_daughter_mean_dedx_pl1_rebuild, "reco_track_daughter_mean_dedx_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_mean_dedx_pl2_rebuild", &reco_track_daughter_mean_dedx_pl2_rebuild, "reco_track_daughter_mean_dedx_pl2_rebuild[20][20]/F");

  fEventTree->Branch("reco_track_daughter_likepr_3pl_rebuild", &reco_track_daughter_likepr_3pl_rebuild, "reco_track_daughter_likepr_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_3pl_rebuild", &reco_track_daughter_llrpid_3pl_rebuild, "reco_track_daughter_llrpid_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_k_3pl_rebuild", &reco_track_daughter_llrpid_k_3pl_rebuild, "reco_track_daughter_llrpid_k_3pl_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC_rebuild", &reco_track_daughter_vtx_inTPC_rebuild, "reco_track_daughter_vtx_inTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC_rebuild", &reco_track_daughter_vtx_in5cmTPC_rebuild, "reco_track_daughter_vtx_in5cmTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC_rebuild", &reco_track_daughter_vtx_inCCInclusiveTPC_rebuild, "reco_track_daughter_vtx_inCCInclusiveTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC_rebuild", &reco_track_daughter_end_inTPC_rebuild, "reco_track_daughter_end_inTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC_rebuild", &reco_track_daughter_end_in5cmTPC_rebuild, "reco_track_daughter_end_in5cmTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC_rebuild", &reco_track_daughter_end_inCCInclusiveTPC_rebuild, "reco_track_daughter_end_inCCInclusiveTPC_rebuild[20][20]/O");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl0_rebuild", &reco_track_daughter_Bragg_fwd_ka_pl0_rebuild, "reco_track_daughter_Bragg_fwd_ka_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl0_rebuild", &reco_track_daughter_Bragg_fwd_pr_pl0_rebuild, "reco_track_daughter_Bragg_fwd_pr_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl0_rebuild", &reco_track_daughter_Bragg_fwd_pi_pl0_rebuild, "reco_track_daughter_Bragg_fwd_pi_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl0_rebuild", &reco_track_daughter_Bragg_fwd_mu_pl0_rebuild, "reco_track_daughter_Bragg_fwd_mu_pl0_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl1_rebuild", &reco_track_daughter_Bragg_fwd_ka_pl1_rebuild, "reco_track_daughter_Bragg_fwd_ka_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl1_rebuild", &reco_track_daughter_Bragg_fwd_pr_pl1_rebuild, "reco_track_daughter_Bragg_fwd_pr_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl1_rebuild", &reco_track_daughter_Bragg_fwd_pi_pl1_rebuild, "reco_track_daughter_Bragg_fwd_pi_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl1_rebuild", &reco_track_daughter_Bragg_fwd_mu_pl1_rebuild, "reco_track_daughter_Bragg_fwd_mu_pl1_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl2_rebuild", &reco_track_daughter_Bragg_fwd_ka_pl2_rebuild, "reco_track_daughter_Bragg_fwd_ka_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl2_rebuild", &reco_track_daughter_Bragg_fwd_pr_pl2_rebuild, "reco_track_daughter_Bragg_fwd_pr_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl2_rebuild", &reco_track_daughter_Bragg_fwd_pi_pl2_rebuild, "reco_track_daughter_Bragg_fwd_pi_pl2_rebuild[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl2_rebuild", &reco_track_daughter_Bragg_fwd_mu_pl2_rebuild, "reco_track_daughter_Bragg_fwd_mu_pl2_rebuild[20][20]/F");



  fEventTree->Branch("k_can_trkid", &k_can_trkid,"k_can_trkid/I");
  fEventTree->Branch("mu_can_trkid", &mu_can_trkid,"mu_can_trkid/I");
  fEventTree->Branch("k_mu_can_dis", &k_mu_can_dis,"k_mu_can_dis/F");
  fEventTree->Branch("k_mu_open_angle", &k_mu_open_angle,"k_mu_open_angle/F");
  fEventTree->Branch("k_vtx_dis", &k_vtx_dis,"k_vtx_dis/F");

  fSubrunTree = tfs->make<TTree>("subruns", "SubRun Tree");
  fSubrunTree->Branch("run", &m_run, "run/i");
  fSubrunTree->Branch("subRun", &m_subrun, "subRun/i");
  //if (!m_isData)
  fSubrunTree->Branch("pot", &m_pot, "pot/F");
    
}

void CCKaonAnalyzer::analyze( const art::Event& evt){
  reset();  
  run=evt.run();
  subrun=evt.subRun();
  event=evt.id().event(); 
  cout << "Run " << run;
  cout << " Subrun " << subrun;
  cout << " Event " << event << endl;

  art::Handle< std::vector<simb::MCTruth> > mctruths;
  std::vector< art::Ptr<simb::MCParticle> > ptList;


  //-----------scan--------------
  
    n_prip=0;
    n_prip_k_dau=0;
    Int_t prim_k_id = -999;
    prip = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
    prip_k_dau = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
    IsKaon = false;
    IsSingleKaon = false;
    IsAssociatedKaon = false;
    IsMuBR = false;
    IsPiBR = false;
    flag_prim_k = true;
    

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
    


    // print true information
    for (auto const& pPart : ptList) {

      /*      
      cout << "trackId " << pPart->TrackId();
      cout << ", mother " << pPart->Mother();
      cout << ", pdg " << pPart->PdgCode();
      cout << ", process " << pPart->Process() << " - " << pPart->EndProcess();
      cout << ", ke " << pPart->E()-pPart->Mass() << " - " << pPart->EndE()-pPart->Mass();
      cout << ", P " << pPart->Momentum().Vect().Mag();
      cout << ", Eng " << pPart->E() << ", EndE " << pPart->EndE();
      cout << ", X " << pPart->Vx() << " - " << pPart->EndX();
      cout << ", Y " << pPart->Vy() << " - " << pPart->EndY();
      cout << ", Z " << pPart->Vz() << " - " << pPart->EndZ();
      cout << ", EndPointx " << pPart->EndPosition()[0];
      cout << ", EndPointy " << pPart->EndPosition()[1];
      cout << ", EndPointz " << pPart->EndPosition()[2];
      cout << ", delta Z " << pPart->EndZ()-pPart->Vz();
      cout << ", ndaughters " << pPart->NumberDaughters();
      cout << endl;
      */
      
      //cout << n_prip << endl;
      if(n_prip<20 && pPart->Mother()==0){//primary interaction

	prip.at(n_prip) = pPart->PdgCode();

	if( (flag_prim_k ==true)  && (prip.at(n_prip) == 321)){//kaon
	  prim_k_id = pPart->TrackId();
	  flag_prim_k = false;
	}
	++n_prip;

      }

      else continue;
      //++n_prip;  

    }
    n_prip=0;



    for (auto const& pPart : ptList) {  

      if(n_prip_k_dau<20 && pPart->Mother()==prim_k_id){
	prip_k_dau.at(n_prip_k_dau) = pPart->PdgCode();
	++n_prip_k_dau;
      }

      else continue;
      //++n_prip_k_dau;      

    }
    n_prip_k_dau=0;
    flag_prim_k = true;


    if(IsKaon == false){

      //if(std::find(prip.begin(),prip.end(),321)!=prip.end() && std::find(prip.begin(),prip.end(),-321)==prip.end() && std::find(prip.begin(),prip.end(),13)!=prip.end()){
    if(std::find(prip.begin(),prip.end(),321)!=prip.end() && std::find(prip.begin(),prip.end(),13)!=prip.end()){
      IsKaon = true;
      cout << "IS K+" << endl;

      if( std::find(prip.begin(),prip.end(),3212)!=prip.end() || std::find(prip.begin(),prip.end(),3222)!=prip.end() || std::find(prip.begin(),prip.end(),3112)!=prip.end() || std::find(prip.begin(),prip.end(),3122)!=prip.end()){
	IsAssociatedKaon = true;
	
	if(true_nu_pdg==14 && true_nu_ccnc==0 && true_nu_vtx_inCCInclusiveTPC==1){
	  
	  cout << "IS ASSOCIATED KAON" << endl;
	  cout << "Run " << run;
	  cout << " Subrun " << subrun;
	  cout << " Event " << event << endl;
	  cout << "Primary K+ ID: " << prim_k_id << endl;
	  
	}
      }
      else{
	IsSingleKaon = true;
	if(true_nu_pdg==14 && true_nu_ccnc==0 && true_nu_vtx_inCCInclusiveTPC==1){
	//if(true_nu_pdg==14 && true_nu_ccnc==0 && true_kaon_ke>=0 && true_nu_vtx_inCCInclusiveTPC==1){
	  cout << "IS SINGLE KAON!!!!" << endl;
	  cout << "Run " << run;
	  cout << " Subrun " << subrun;
	  cout << " Event " << event << endl;
	  cout << "Primary K+ ID: " << prim_k_id << endl;
	  
	  }
      }

      if( std::find(prip_k_dau.begin(),prip_k_dau.end(),-13)!=prip_k_dau.end()){
	IsMuBR = true;
	cout << "IS KAON -> MU" << endl;
      }
      else if( std::find(prip_k_dau.begin(),prip_k_dau.end(),211)!=prip_k_dau.end() ){
	IsPiBR = true;
	cout << "IS KAON -> PI" << endl;
      }

    }

    }

  //-----------------------------


  if (isMC) {

    // get MC event weights
    art::InputTag eventweight_tag("eventweight");
    art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
    
    if (evt.getByLabel(eventweight_tag, eventweights_handle)) {
      std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
      art::fill_ptr_vector(eventweights, eventweights_handle);
      std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
      //std::cout << "Weight: " <<std::endl;
      for (auto const& x : evtwgt_map) {
        //std::cout << "Weight - Size fun:\t" << x.first << " " << x.second.size() << std::endl;
        std::vector<float> wg;
      //}
      const std::vector<double> &weights = evtwgt_map.at(x.first);
      //event_weight = weights.front();
      //std::cout << "Event Weight:  " << event_weight << std::endl;
      for (unsigned int i=0;i<weights.size();i++) {
        //cout << i << "\t" << weights[i] << endl;
        wg.push_back(weights[i]);
      }

      event_weight[x.first] = wg;
      evtwgt_funcname.push_back(x.first);
      evtwgt_weight.push_back(wg);
      evtwgt_nweight.push_back(x.second.size());

      }
    }
    else {
      //std::cout << "***** Failed obtaining eventweight ******" << std::endl;
    }

    // get MCTruth
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //std::cout << "Number of MCTruths objects in event " << mctruths->size() << std::endl;
      return;
      }
      simb::MCTruth mctruth = mctruths->at(0);

    // get MCParticles
    /*
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
    if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
      art::fill_ptr_vector(ptList, mcParticleHandle); 
    }
    */

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

    // print true information
    /*
    cout << mctruths->at(0) << endl;
    for (auto const& pPart : ptList) {
      if (pPart->Process()!="Primary" && pPart->E()-pPart->Mass()<0.001) continue;
      cout << "trackId " << pPart->TrackId();
      cout << ", mother " << pPart->Mother();
      cout << ", pdg " << pPart->PdgCode();
      cout << ", process " << pPart->Process() << " - " << pPart->EndProcess();
      cout << ", ke " << pPart->E()-pPart->Mass() << " - " << pPart->EndE()-pPart->Mass();
      cout << ", X " << pPart->Vx() << " - " << pPart->EndX();
      cout << ", Z " << pPart->Vz() << " - " << pPart->EndZ();
      cout << ", delta Z " << pPart->EndZ()-pPart->Vz();
      cout << ", ndoughters " << pPart->NumberDaughters();
      cout << endl;
    }
    */

    // find the the highest momentum k+
    double true_kaon_maxp = -1;
    int true_kaon_maxp_id = -1;
    true_nkaons = 0;
    true_nprimary = 0;
    true_nhyperons = 0;
    //double true_hyperon_maxp = -1;
    //int true_hyperon_maxp_id = -1;

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

      // -- cout << "True kaon end process " << true_kaon_end_process << endl;

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

  
    // find the the highest momentum k-
    double true_kaon_maxp_anti = -1;
    int true_kaon_maxp_id_anti = -1;
    true_nkaons_anti = 0;

    for (auto const& pPart : ptList) {
      if (pPart->Process()=="primary" && pPart->PdgCode()==-321) {
        if(int(mctruths->at(0).Origin())==1) { // make sure origin is neutrino beam
          if (true_kaon_maxp_anti < pPart->P()) {
            true_kaon_maxp_anti = pPart->P();
            true_kaon_maxp_id_anti = pPart->TrackId();
            true_nkaons_anti++;
          }
        }  
      }
    }
  
    // true signal interaction: nu CC with at least one k- // selection with true info
    if (true_nu_pdg==14 && true_nu_ccnc==0 && true_nkaons_anti>0) {

      vector<int> decay_daughters_anti;
      vector<int> decay_muplus_id_anti;
      vector<int> decay_piplus_id_anti;
      vector<int> inelastic_daughters_anti;
      vector<int> inelastic_kaplus_id_anti;
      TVector3 true_kaon_pvector_anti;

      for (auto const& pPart : ptList) {

        // primary kaon- (w/ highest momentum)
        if (pPart->TrackId()==true_kaon_maxp_id_anti) {
          TLorentzVector mcstart_anti, mcend_anti;
          unsigned int pstarti_anti, pendi_anti;
          true_kaon_length_anti = length(*pPart, mcstart_anti, mcend_anti, pstarti_anti, pendi_anti);
          true_kaon_p_anti  = pPart->P();
          true_kaon_ke_anti = pPart->E()-pPart->Mass();
          true_kaon_theta_anti = pPart->Momentum().Theta();
          true_kaon_costheta_anti = pPart->Momentum().CosTheta();
          true_kaon_phi_anti = pPart->Momentum().Phi();
          true_kaon_ccmuon_angle_anti = pPart->Momentum().Angle(true_lepton_pvector);
          true_kaon_ccmuon_cosangle_anti = TMath::Cos(pPart->Momentum().Angle(true_lepton_pvector));
          true_kaon_pvector_anti = pPart->Momentum().Vect();
          true_kaon_end_ke_anti = pPart->EndE()-pPart->Mass();
          true_kaon_end_x_anti = pPart->EndX();
          true_kaon_end_y_anti = pPart->EndY();
          true_kaon_end_z_anti = pPart->EndZ();
          true_kaon_end_inTPC_anti = isInsideVolume("TPC",pPart->EndPosition().Vect());
          true_kaon_end_in5cmTPC_anti = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
          true_kaon_end_inCCInclusiveTPC_anti = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
        }

        // find kaon- daughters
        else if (pPart->Mother()==true_kaon_maxp_id_anti) {

          // decay daughters
          if (pPart->Process()=="Decay") {
            decay_daughters_anti.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==13) {
              decay_muplus_id_anti.push_back(pPart->TrackId());
            }
            else if (pPart->PdgCode()==-211) {
              decay_piplus_id_anti.push_back(pPart->TrackId());
            }
            //cout << "Decay " << pPart->PdgCode() << endl;
          }

          // inelastic interaction daughters
          if (pPart->Process()=="kaon+Inelastic") {
            inelastic_daughters_anti.push_back(pPart->PdgCode());
            if (pPart->PdgCode()==-321) {
              inelastic_kaplus_id_anti.push_back(pPart->TrackId());
            }
            //cout << "Inelastic " << pPart->PdgCode() << endl;
          }

        }

      }//MC particles loop
    
      true_kaon_ndaughters_anti = decay_daughters_anti.size() + inelastic_daughters_anti.size();
      true_kaon_ndaughters_decay_anti = decay_daughters_anti.size();
      true_kaon_ndaughters_inelastic_anti = inelastic_daughters_anti.size();

      int n_decay_muplus_anti     = decay_muplus_id_anti.size();
      int n_decay_piplus_anti     = decay_piplus_id_anti.size();
      int n_decay_numu_anti       = count(decay_daughters_anti.begin(),decay_daughters_anti.end(),-14);
      int n_decay_pi0_anti        = count(decay_daughters_anti.begin(),decay_daughters_anti.end(),111);
      int n_inelastic_kaplus_anti = inelastic_kaplus_id_anti.size();
      int n_inelastic_piplus_anti = count(inelastic_daughters_anti.begin(),inelastic_daughters_anti.end(),-211);
      int n_inelastic_proton_anti = count(inelastic_daughters_anti.begin(),inelastic_daughters_anti.end(),-2212);

      true_kaon_ndecmup_anti = n_decay_muplus_anti;
      true_kaon_ndecpip_anti = n_decay_piplus_anti;
      true_kaon_ninekap_anti = n_inelastic_kaplus_anti;
      true_kaon_ninepip_anti = n_inelastic_piplus_anti;
      true_kaon_ninepro_anti = n_inelastic_proton_anti;

      // define kaon end process
      int daughter_id_anti = -1;
      if (decay_daughters_anti.size()) {
        //kaon+ -> muon+ numu
        if (decay_daughters_anti.size()==2 && n_decay_muplus_anti==1 && n_decay_numu_anti==1) {
          true_kaon_end_process_anti = 0;
          daughter_id_anti = decay_muplus_id_anti.front();
        }
        //kaon+ -> pion+ + pi0
        else if (decay_daughters_anti.size()==2 && n_decay_piplus_anti==1 && n_decay_pi0_anti==1) {
          true_kaon_end_process_anti = 1;
          daughter_id_anti = decay_piplus_id_anti.front();
        }
        //kaon+ -> muon+ nopion+
        else if (n_decay_muplus_anti==1 && n_decay_piplus_anti==0) {
          true_kaon_end_process_anti = 2;
          daughter_id_anti = decay_muplus_id_anti.front();
        }
        //kaon+ -> pion+ and nomuon+
        else if (n_decay_muplus_anti==0 && n_decay_piplus_anti==1) {
          true_kaon_end_process_anti = 3;
          daughter_id_anti = decay_piplus_id_anti.front();
        }
        //other decays
        else {
          true_kaon_end_process_anti = 4;
        }
      }

      if (inelastic_daughters_anti.size()) {
        //inelasitic with kaon-
        if (n_inelastic_kaplus_anti==1) {
          true_kaon_end_process_anti = 5;
          daughter_id_anti = inelastic_kaplus_id_anti.front();
        }
        //other inelastic
        else {
          true_kaon_end_process_anti = 6;
        }
      }

      //this shouldn't happen but it does! rarely
      if (decay_daughters_anti.size() && inelastic_daughters_anti.size()) {
        // -- cout << "Primary kaon+ decay and inelastic at the same time!!!" << endl;
        true_kaon_end_process_anti = 7;
      }

      // -- cout << "True kaon end process " << true_kaon_end_process << endl;

      //check kaon daughter if it exists
      //only if there is one pion+, one muon+ or one kaon+
      if (daughter_id_anti!=-1) {
        for (auto const& pPart : ptList) {
      
          // kaon- daughter
          if (pPart->TrackId()==daughter_id_anti) {
            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            true_kaon_daughter_length_anti = length(*pPart, mcstart, mcend, pstarti, pendi);
            true_kaon_daughter_p_anti = pPart->P();
            true_kaon_daughter_ke_anti = pPart->E()-pPart->Mass();
            true_kaon_daughter_theta_anti = pPart->Momentum().Theta();
            true_kaon_daughter_costheta_anti = pPart->Momentum().CosTheta();
            true_kaon_daughter_angle_anti = pPart->Momentum().Angle(true_kaon_pvector_anti);
            true_kaon_daughter_cosangle_anti = TMath::Cos(pPart->Momentum().Angle(true_kaon_pvector_anti));
            true_kaon_daughter_pdg_anti = pPart->PdgCode();
            true_kaon_daughter_end_x_anti = pPart->EndX();
            true_kaon_daughter_end_y_anti = pPart->EndY();
            true_kaon_daughter_end_z_anti = pPart->EndZ();
            true_kaon_daughter_end_inTPC_anti = isInsideVolume("TPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_in5cmTPC_anti = isInsideVolume("5cmTPC",pPart->EndPosition().Vect());
            true_kaon_daughter_end_inCCInclusiveTPC_anti = isInsideVolume("CCInclusiveTPC",pPart->EndPosition().Vect());
            break;
          }
      
        }//MC particles loop
      }//kaon daughter exists

    }//is true signal

    //this->filter(evt);
  
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


  /*
  if (!reco_nu_cc_filter) {
    fEventTree->Fill();
    return;
  }
  */

  // Collect all recontructed particles
  art::Handle<std::vector<recob::PFParticle>> pfparticles;
  evt.getByLabel(m_pfp_producer, pfparticles);

  if (pfparticles->size()==0) {
    std::cout << "No PFParticles found" << std::endl;
    fEventTree->Fill();
    return;
  }

  // Get PFParticle associations
  art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
  if(!pfp_muon_assn.isValid()){
    // -- cout << "PFParticle-T0 associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
  if(!pfparticleTrackAssn.isValid()){
    // -- cout << "PFParticle-Track associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
  if(!pfparticleVertexAssn.isValid()){
    // -- cout << "PFParticle-Vertex associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

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

  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

  reco_nu_vtx_x = pfparticleVertexAssn.at(pfnu.key()).front()->position().X();
  reco_nu_vtx_y = pfparticleVertexAssn.at(pfnu.key()).front()->position().Y();
  reco_nu_vtx_z = pfparticleVertexAssn.at(pfnu.key()).front()->position().Z();
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

      //reco_nu_daughters_id.push_back(track.key());
      reco_nu_daughters_id.push_back(track->ID());

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
  //art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.ID()).front();

  reco_ccmu_vtx_x = pfparticleVertexAssn.at(pfmuon.key()).front()->position().X();
  reco_ccmu_vtx_y = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Y();
  reco_ccmu_vtx_z = pfparticleVertexAssn.at(pfmuon.key()).front()->position().Z();
  //--cout << "CC muon start (x,y,z) = " << reco_ccmu_vtx_x << ", " << reco_ccmu_vtx_y << ", " << reco_ccmu_vtx_z << endl;
  reco_ccmu_vtx_inTPC = isInsideVolume("TPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_in5cmTPC = isInsideVolume("5cmTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);
  reco_ccmu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",reco_ccmu_vtx_x,reco_ccmu_vtx_y,reco_ccmu_vtx_z);

  // get hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitlist, hitListHandle);
  }
  //--cout << "Found " << hitlist.size() << " hits" << endl;

  // get tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }

  art::Handle< std::vector<recob::Track> > rebuilttrackListHandle;
  std::vector<art::Ptr<recob::Track> > rebuilttracklist;
  if(evt.getByLabel(fRebuiltTrackModuleLabel,rebuilttrackListHandle)) {
    art::fill_ptr_vector(rebuilttracklist, rebuilttrackListHandle);
  }
  //--cout << "Found " << tracklist.size() << " tracks" << endl;

  // get showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
    art::fill_ptr_vector(showerlist, showerListHandle);
  }

  // get track associations
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  if(!fmcal.isValid()) cout << "fmcal is invalid" << endl;
  if(!fmcal.isValid()){
    //--cout << "Track-Calorimetry associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }
  art::FindManyP<anab::Calorimetry> fmcal_rebuilt(rebuilttrackListHandle, evt, fRebuiltCalorimetryModuleLabel);

  art::FindManyP<recob::Hit> hits_from_rebuilttracks(rebuilttrackListHandle, evt, fRebuiltHitTrackAssns);
  art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
  art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
  if(!hits_from_tracks.isValid()){
    //--cout << "Track-Hit associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  art::FindManyP<anab::ParticleID> rebuilttrackPIDAssn(rebuilttrackListHandle, evt, fRebuiltPIDLabel);
  art::FindManyP<anab::ParticleID> trackPIDAssn(trackListHandle, evt, fPIDLabel);
  if(!trackPIDAssn.isValid()){
    //--cout << "Track PID associations are not valid" << endl;
    fEventTree->Fill();
    return;
  }

  // find track multiplicity around vertex
  int NTracks=tracklist.size();
  int RebuiltNTracks=rebuilttracklist.size();
  int ntracks = 0;

  int NShowers=showerlist.size();


  //selection cuts
  std::vector<int> kaon_can_trkID;
  std::vector<int> muon_can_trkID;

  std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>> assocMCPartt;
  if (isMC) {
    assocMCPartt = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> (new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, evt, fHitTruthAssns));
  }

 
 // loop over nu daughters
  cout << "reco_nu_ndaughters: " << reco_nu_ndaughters << endl;

  for (int i=0; i<NTracks; i++) {


    art::Ptr<recob::Track> ptrack(trackListHandle,i);
    const recob::Track& track = *ptrack;

    // skip cc muon track
    if (track.ID()==trkmuon->ID()){
      cout << "THIS IS CC MUON" << endl;
      continue;
    }

    cout << "track.key(): " << ptrack.key() << ", track.ID(): "<< track.ID() << endl;
    cout << "trkmuon.key(): " << trkmuon.key() << endl;

    // take all primary tracks    
    bool skip = true;
    for (int k=0; k<reco_nu_ndaughters; k++) {
	cout << "reco_nu_daughters_id[k]: " << reco_nu_daughters_id[k] << endl;

      if (int(track.ID())==reco_nu_daughters_id[k]) {

      //if (int(ptrack.key())!=reco_nu_daughters_id[k]) {
	skip=false;
	break;
      }
    }
    if (skip) continue;

    cout << "THIS TRACK IS NOT SKIPPED" << endl;

    // check track start and end
    TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
    TVector3 end(track.End().X(),track.End().Y(),track.End().Z());

    reco_track_start_x[ntracks] = track.Vertex().X();
    reco_track_start_y[ntracks] = track.Vertex().Y();
    reco_track_start_z[ntracks] = track.Vertex().Z();

    reco_track_end_x[ntracks] = track.End().X();
    reco_track_end_y[ntracks] = track.End().Y();
    reco_track_end_z[ntracks] = track.End().Z();

    double st_vtx=TMath::Sqrt((reco_nu_vtx_x-pos.X())*(reco_nu_vtx_x-pos.X()) +
                              (reco_nu_vtx_y-pos.Y())*(reco_nu_vtx_y-pos.Y()) +
                              (reco_nu_vtx_z-pos.Z())*(reco_nu_vtx_z-pos.Z()));
    reco_track_distance[ntracks] = st_vtx;//distance between vtx and track start position

    reco_track_vtx_inTPC[ntracks] = isInsideVolume("TPC",pos);
    reco_track_vtx_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",pos);
    reco_track_vtx_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",pos);

    reco_track_end_inTPC[ntracks] = isInsideVolume("TPC",end);
    reco_track_end_in5cmTPC[ntracks] = isInsideVolume("5cmTPC",end);
    reco_track_end_inCCInclusiveTPC[ntracks] = isInsideVolume("CCInclusiveTPC",end);

    // track length and angles
    double trklen=track.Length();
    reco_track_length[ntracks] = trklen;
    reco_track_theta[ntracks] = track.Theta();
    reco_track_phi[ntracks] = track.Phi();

    reco_track_P_vtx[ntracks] = track.VertexMomentum();
    reco_track_P_str[ntracks] = track.StartMomentum();
    reco_track_P_end[ntracks] = track.EndMomentum();

    double end_dis=TMath::Sqrt((reco_nu_vtx_x-end.X())*(reco_nu_vtx_x-end.X()) +
                               (reco_nu_vtx_y-end.Y())*(reco_nu_vtx_y-end.Y()) +
                               (reco_nu_vtx_z-end.Z())*(reco_nu_vtx_z-end.Z()));
    reco_track_dir[ntracks] = (st_vtx<end_dis);

    //fillCalorimetry(fmcal.at(ptrack.key()),ntracks);
    //fillCalorimetry(fmcal.at(ptrack.key()),ntracks,-1,false);//no daughter track id
    fillCalorimetry(fmcal.at(ptrack.key()), ptrack, ntracks, -1, false);//no daughter track id
  
    // check PID
    if (trackPIDAssn.isValid()){
      double delta_z = end.Z()-pos.Z();
      double delta_y = end.Y()-pos.Y();
      double angle_y = TMath::ATan2(delta_z, delta_y);
      //fillPID(trackPIDAssn.at(ptrack.key()), angle_y, ntracks);
      fillPID(trackPIDAssn.at(ptrack.key()), angle_y, ntracks, -1, false);
    }

    // find true matched particle
    if (isMC) {
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
      //assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> (new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, evt, fHitTruthAssns));
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      fillTrueMatching(hits_from_track, particles_per_hit, ntracks, -1, true, false);
      mergeChecker(hits_from_track, particles_per_hit, ntracks, -1, true, false);

    }//isMC

    int ndaughters = 0;
    for (int j=0; j<NTracks; j++) {
      art::Ptr<recob::Track> ptrack_dau(trackListHandle,j);
      const recob::Track& track_dau = *ptrack_dau;
    
      // skip all primary tracks
      if (track_dau.ID()==trkmuon->ID()){
	continue;
      }
      
      bool skip = false;
      for (int k=0; k<reco_nu_ndaughters; k++) {
        //if (int(ptrack_dau.key())==reco_nu_daughters_id[k]) {
        if (int(track_dau.ID())==reco_nu_daughters_id[k]) {
          skip=true;
          break;
        }
      }
      if (skip) continue;
     

      TVector3 pos2(track_dau.Vertex().X(),track_dau.Vertex().Y(),track_dau.Vertex().Z());

      double track_dau_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
                                         (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
                                         (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

      // check distance to vertex
      if (track_dau_distance<10) { //7cm, 

        reco_track_daughter_distance[ntracks][ndaughters] = track_dau_distance;

        TVector3 end2(track_dau.End().X(),track_dau.End().Y(),track_dau.End().Z());

	reco_track_daughter_start_x[ntracks][ndaughters] = track_dau.Vertex().X();
	reco_track_daughter_start_y[ntracks][ndaughters] = track_dau.Vertex().Y();
	reco_track_daughter_start_z[ntracks][ndaughters] = track_dau.Vertex().Z();
	reco_track_daughter_end_x[ntracks][ndaughters] = track_dau.End().X();
	reco_track_daughter_end_y[ntracks][ndaughters] = track_dau.End().Y();
	reco_track_daughter_end_z[ntracks][ndaughters] = track_dau.End().Z();

        reco_track_daughter_vtx_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",pos2);
        reco_track_daughter_vtx_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",pos2);
        reco_track_daughter_vtx_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",pos2);
      
        reco_track_daughter_end_inTPC[ntracks][ndaughters] = isInsideVolume("TPC",end2);
        reco_track_daughter_end_in5cmTPC[ntracks][ndaughters] = isInsideVolume("5cmTPC",end2);
        reco_track_daughter_end_inCCInclusiveTPC[ntracks][ndaughters] = isInsideVolume("CCInclusiveTPC",end2);

        // track length and angles
        reco_track_daughter_length[ntracks][ndaughters] = track_dau.Length();
        reco_track_daughter_theta[ntracks][ndaughters] = track_dau.Theta();
        reco_track_daughter_phi[ntracks][ndaughters] = track_dau.Phi();

        //fillCalorimetry_old(fmcal.at(ptrack_dau_old.key()),ntracks,ndaughters_old);
        fillCalorimetry(fmcal.at(ptrack_dau.key()), ptrack_dau, ntracks, ndaughters, false);

        // check PID
        if (rebuilttrackPIDAssn.isValid()) {
          double delta_z = end2.Z()-pos2.Z();
          double delta_y = end2.Y()-pos2.Y();
          double angle_y = TMath::ATan2(delta_z, delta_y);
          //fillPID_old(trackPIDAssn.at(ptrack_dau_old.key()), angle_y, ntracks, ndaughters_old);
          fillPID(trackPIDAssn.at(ptrack_dau.key()), angle_y, ntracks, ndaughters, false);
        }      

	if(isMC){
          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
          std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack_dau.key());
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  fillTrueMatching(hits_from_track, particles_per_hit, ntracks, ndaughters, true, false);
	  mergeChecker(hits_from_track, particles_per_hit, ntracks, ndaughters, true, false); 
	}
      
        ndaughters++;

      } // remove condition on gap between K+ end and daughter start
    
      reco_track_ndaughters[ntracks] = ndaughters;
 
    }


    // check if there is a track at the end
    int ndaughters_rebuild = 0;
    for (int j=0; j<RebuiltNTracks; j++) {
    
      art::Ptr<recob::Track> ptrack_dau(rebuilttrackListHandle,j);
      const recob::Track& track_dau = *ptrack_dau;
    
      // skip all primary tracks
      if (track_dau.ID()==trkmuon->ID()){
	continue;
      }
      //if (track_dau.ID()==trkmuon->ID()) continue;
      
      bool skip = false;
      for (int k=0; k<reco_nu_ndaughters; k++) {
        //if (int(ptrack_dau.key())==reco_nu_daughters_id[k]) {
        if (int(track_dau.ID())==reco_nu_daughters_id[k]) {
          skip=true;
          break;
        }
      }
      if (skip) continue;
      
      TVector3 pos2(track_dau.Vertex().X(),track_dau.Vertex().Y(),track_dau.Vertex().Z());

      double track_dau_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
					    (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
					    (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

      // check distance to vertex
      if (track_dau_distance<10) { //7cm, 
        reco_track_daughter_distance_rebuild[ntracks][ndaughters_rebuild] = track_dau_distance;

        TVector3 end2(track_dau.End().X(),track_dau.End().Y(),track_dau.End().Z());

	reco_track_daughter_start_x_rebuild[ntracks][ndaughters_rebuild] = track_dau.Vertex().X();
	reco_track_daughter_start_y_rebuild[ntracks][ndaughters_rebuild] = track_dau.Vertex().Y();
	reco_track_daughter_start_z_rebuild[ntracks][ndaughters_rebuild] = track_dau.Vertex().Z();
	reco_track_daughter_end_x_rebuild[ntracks][ndaughters_rebuild] = track_dau.End().X();
	reco_track_daughter_end_y_rebuild[ntracks][ndaughters_rebuild] = track_dau.End().Y();
	reco_track_daughter_end_z_rebuild[ntracks][ndaughters_rebuild] = track_dau.End().Z();

        reco_track_daughter_vtx_inTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("TPC",pos2);
        reco_track_daughter_vtx_in5cmTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("5cmTPC",pos2);
        reco_track_daughter_vtx_inCCInclusiveTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("CCInclusiveTPC",pos2);
      
        reco_track_daughter_end_inTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("TPC",end2);
        reco_track_daughter_end_in5cmTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("5cmTPC",end2);
        reco_track_daughter_end_inCCInclusiveTPC_rebuild[ntracks][ndaughters_rebuild] = isInsideVolume("CCInclusiveTPC",end2);

        // track length and angles
        reco_track_daughter_length_rebuild[ntracks][ndaughters_rebuild] = track_dau.Length();
        reco_track_daughter_theta_rebuild[ntracks][ndaughters_rebuild] = track_dau.Theta();
        reco_track_daughter_phi_rebuild[ntracks][ndaughters_rebuild] = track_dau.Phi();

        double start_dis=TMath::Sqrt((reco_nu_vtx_x-pos2.X())*(reco_nu_vtx_x-pos2.X()) +
                                     (reco_nu_vtx_y-pos2.Y())*(reco_nu_vtx_y-pos2.Y()) +
                                     (reco_nu_vtx_z-pos2.Z())*(reco_nu_vtx_z-pos2.Z()));    

        reco_track_daughter_vtx_distance_rebuild[ntracks][ndaughters_rebuild] = start_dis;

	//cout << "track length is " << track_dau.Length() << endl;
        //fillCalorimetry(fmcal_rebuilt.at(ptrack_dau.key()),ntracks,ndaughters_rebuild,true);
        fillCalorimetry(fmcal_rebuilt.at(ptrack_dau.key()), ptrack_dau, ntracks, ndaughters_rebuild, true);

        // check PID
        if (rebuilttrackPIDAssn.isValid()) {
          double delta_z = end2.Z()-pos2.Z();
          double delta_y = end2.Y()-pos2.Y();
          double angle_y = TMath::ATan2(delta_z, delta_y);
          //fillPID(rebuilttrackPIDAssn.at(ptrack_dau.key()), angle_y, ntracks, ndaughters);
          fillPID(rebuilttrackPIDAssn.at(ptrack_dau.key()), angle_y, ntracks, ndaughters_rebuild, true);
        }      

        // find true matched particle
        if (isMC) {
          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
          std::vector<art::Ptr<recob::Hit>> hits_from_rebuilttrack = hits_from_rebuilttracks.at(ptrack_dau.key());
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  fillTrueMatching(hits_from_rebuilttrack, particles_per_hit, ntracks, ndaughters, true, true);
	  mergeChecker(hits_from_rebuilttrack, particles_per_hit, ntracks, ndaughters, true, true);
	}//isMC

        ndaughters_rebuild++;

      } // remove condition on gap between K+ end and daughter start
    
      reco_track_ndaughters_rebuild[ntracks] = ndaughters_rebuild;
          	
    }//NTrack loop

    int ndaughters_sh = 0;
    for (int j_s=0; j_s<NShowers; j_s++) {
      art::Ptr<recob::Shower> pshower_dau(showerListHandle,j_s);
      const recob::Shower& shower_dau = *pshower_dau;
      
      // skip all primary showers

      bool skip2 = false;
      for (int k_sh=0; k_sh<reco_nu_ndaughters; k_sh++) {
	if (int(pshower_dau.key())==reco_nu_daughters_id[k_sh]) {
	  skip2=true;
	  break;
	}
      }
      if (skip2) continue;
      TVector3 pos2_sh(shower_dau.ShowerStart().X(),shower_dau.ShowerStart().Y(),shower_dau.ShowerStart().Z());

      double shower_dau_distance=TMath::Sqrt((end.X()-pos2_sh.X())*(end.X()-pos2_sh.X()) +
					     (end.Y()-pos2_sh.Y())*(end.Y()-pos2_sh.Y()) +
					     (end.Z()-pos2_sh.Z())*(end.Z()-pos2_sh.Z()));

      if (shower_dau_distance<40) {
	if(ndaughters_sh>20) break;
        if (isMC) {
          art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, evt, fHitTruthAssns);
          std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_rebuilttracks.at(pshower_dau.key());
	  fillTrueMatching(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh, false, false);
	  mergeChecker(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh, false, false);
	}
	ndaughters_sh++;
      }
    }

    ntracks++;

  } // loop over K  reco trks	

  reco_ntracks = ntracks;
  
  if(flg_reco_k == true){
    cout << "Event has reconstrcted K+ track" << endl;
    //if(true_kaon_end_process == 0){
    //cout << "K+ -> Mu+ decay" << endl;
    if(flg_reco_mu == true) cout << "Event has reconstrcted mu+ track" << endl;
    else cout << "Event DOES NOT have reconstrcted mu+ track" << endl;
    //}
    //if(true_kaon_end_process == 1){
    //cout << "K+ -> Pi+ decay" << endl;
    if(flg_reco_pi == true) cout << "Event has reconstrcted pi+ track" << endl; 
    else cout << "Event DOES NOT have reconstrcted pi+ track" << endl;
    //}
  }
  else cout << "Event DO NOT have reconstrcted K+ track" << endl;
  
  flg_reco_k  = false;
  flg_reco_mu  = false;
  flg_reco_pi  = false;
  

  fEventTree->Fill();


} // end of analyze function
 

 /////////////////////////////////////////// Reset Function ///////////////////////////////
 
 void CCKaonAnalyzer::reset(){

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
   
   true_nkaons_anti=-9999;
   true_kaon_length_anti=-9999;
   true_kaon_p_anti=-9999;
   true_kaon_ke_anti=-9999;
   true_kaon_theta_anti=-9999;
   true_kaon_costheta_anti=-9999;
   true_kaon_phi_anti=-9999;
   true_kaon_ccmuon_angle_anti=-9999;
   true_kaon_ccmuon_cosangle_anti=-9999;
   true_kaon_end_process_anti=-9999;
   true_kaon_end_ke_anti=-9999;
   true_kaon_end_x_anti=-9999;
   true_kaon_end_y_anti=-9999;
   true_kaon_end_z_anti=-9999;
   true_kaon_end_inTPC_anti=false;
   true_kaon_end_in5cmTPC_anti=false;
   true_kaon_end_inCCInclusiveTPC_anti=false;
   
   true_kaon_ndaughters_anti=-9999;
   true_kaon_ndaughters_decay_anti=-9999;
   true_kaon_ndaughters_inelastic_anti=-9999;
   true_kaon_ndecmup_anti=-9999;
   true_kaon_ndecpip_anti=-9999;
   true_kaon_ninekap_anti=-9999;
   true_kaon_ninepip_anti=-9999;
   true_kaon_ninepro_anti=-9999;
   true_kaon_daughter_length_anti=-9999;
   true_kaon_daughter_p_anti=-9999;
   true_kaon_daughter_ke_anti=-9999;
   true_kaon_daughter_theta_anti=-9999;
   true_kaon_daughter_costheta_anti=-9999;
   true_kaon_daughter_angle_anti=-9999;
   true_kaon_daughter_cosangle_anti=-9999;
   true_kaon_daughter_pdg_anti=-9999;
   true_kaon_daughter_end_x_anti=-9999;
   true_kaon_daughter_end_y_anti=-9999;
   true_kaon_daughter_end_z_anti=-9999;
   true_kaon_daughter_end_inTPC_anti=false;
   true_kaon_daughter_end_in5cmTPC_anti=false;
   true_kaon_daughter_end_inCCInclusiveTPC_anti=false;
   
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
   
   reco_ntracks=0;
   reco_nshowers=0;
   
   for(int k=0; k<kMaxTracks; k++){


     for(int l=0; l<kMaxParticles; l++){ 
       cheat_peak_pdg[k][l] = -9999;
       cheat_peak_theta[k][l] = -9999;
       cheat_peak_phi[k][l] = -9999;
       best_peak_theta[k][l] = -9999;
       best_peak_phi[k][l] = -9999;
       best_peak_trkln[k][l] = -9999;
     }  

     reco_track_start_x[k]=-9;
     reco_track_start_y[k]=-9;
     reco_track_start_z[k]=-9;
     reco_track_end_x[k]=-9;
     reco_track_end_y[k]=-9;
     reco_track_end_z[k]=-9;
     
     reco_track_distance[k]=-9;
     
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
     
     
     reco_track_chi2ka_3pl[k]=-99;
     reco_track_chi2pr_3pl[k]=-99;
     reco_track_chi2pi_3pl[k]=-99;
     reco_track_chi2mu_3pl[k]=-99;
     reco_track_likepr_3pl[k]=-99;
     reco_track_mean_dedx_3pl[k] = -9;
     reco_track_mean_dedx_pl0[k] = -9;
     reco_track_mean_dedx_pl1[k] = -9;
     reco_track_mean_dedx_pl2[k] = -9;
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
     reco_track_ndaughters_rebuild[k]=0;
     reco_shower_ndaughters[k]=0;
     
     for(int l=0; l<kMaxMerge; l++){ 
       reco_track_match_e[k][l] = -999;
       reco_track_match_hit[k][l] = -999;
       reco_track_match_epdg[k][l] = -999;
       reco_track_match_hitpdg[k][l] = -999;
     }
     
      
     for(int m=0; m<kMaxShowers; m++){ 
       
       reco_track_daughter_start_x[k][m]=-9;
       reco_track_daughter_start_y[k][m]=-9;
       reco_track_daughter_start_z[k][m]=-9;
       reco_track_daughter_end_x[k][m]=-9;
       reco_track_daughter_end_y[k][m]=-9;
       reco_track_daughter_end_z[k][m]=-9;
       
       reco_track_daughter_distance_sh[k][m]=-9;
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

       reco_track_daughter_start_x_rebuild[k][m]=-9;
       reco_track_daughter_start_y_rebuild[k][m]=-9;
       reco_track_daughter_start_z_rebuild[k][m]=-9;
       reco_track_daughter_end_x_rebuild[k][m]=-9;
       reco_track_daughter_end_y_rebuild[k][m]=-9;
       reco_track_daughter_end_z_rebuild[k][m]=-9;       

       reco_track_daughter_distance_rebuild[k][m]=-9;
       reco_track_daughter_vtx_distance_rebuild[k][m]=-9;
       reco_track_daughter_length_rebuild[k][m]=-99;
       reco_track_daughter_theta_rebuild[k][m]=-99;
       reco_track_daughter_phi_rebuild[k][m]=-99;
       reco_track_daughter_chi2ka_pl0_rebuild[k][m]=-9;
       reco_track_daughter_chi2pr_pl0_rebuild[k][m]=-9;
       reco_track_daughter_chi2pi_pl0_rebuild[k][m]=-9;
       reco_track_daughter_chi2mu_pl0_rebuild[k][m]=-9;
       reco_track_daughter_chi2ka_pl1_rebuild[k][m]=-9;
       reco_track_daughter_chi2pr_pl1_rebuild[k][m]=-9;
       reco_track_daughter_chi2pi_pl1_rebuild[k][m]=-9;
       reco_track_daughter_chi2mu_pl1_rebuild[k][m]=-9;
       reco_track_daughter_chi2ka_pl2_rebuild[k][m]=-9;
       reco_track_daughter_chi2pr_pl2_rebuild[k][m]=-9;
       reco_track_daughter_chi2pi_pl2_rebuild[k][m]=-9;
       reco_track_daughter_chi2mu_pl2_rebuild[k][m]=-9;
       reco_track_daughter_chi2ka_3pl_rebuild[k][m]=-99;
       reco_track_daughter_chi2pr_3pl_rebuild[k][m]=-99;
       reco_track_daughter_chi2pi_3pl_rebuild[k][m]=-99;
       reco_track_daughter_chi2mu_3pl_rebuild[k][m]=-99;
       reco_track_daughter_likepr_3pl_rebuild[k][m]=-99;
       reco_track_daughter_mean_dedx_3pl_rebuild[k][m] = -9;
       reco_track_daughter_mean_dedx_pl0_rebuild[k][m] = -9;
       reco_track_daughter_mean_dedx_pl1_rebuild[k][m] = -9;
       reco_track_daughter_mean_dedx_pl2_rebuild[k][m] = -9;
       reco_track_daughter_llrpid_3pl_rebuild[k][m]=-9;
       reco_track_daughter_llrpid_k_3pl_rebuild[k][m]=-9;
       reco_track_daughter_vtx_inTPC_rebuild[k][m]=false;
       reco_track_daughter_vtx_in5cmTPC_rebuild[k][m]=false;
       reco_track_daughter_vtx_inCCInclusiveTPC_rebuild[k][m]=false;
       reco_track_daughter_end_inTPC_rebuild[k][m]=false;
       reco_track_daughter_end_in5cmTPC_rebuild[k][m]=false;
       reco_track_daughter_end_inCCInclusiveTPC_rebuild[k][m]=false;
       reco_track_daughter_true_pdg_rebuild[k][m]=-999;

       reco_track_daughter_Bragg_fwd_ka_pl0_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl0_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl0_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl0_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl1_rebuild[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pr_pl1_rebuild[k][m] = -999;  
       reco_track_daughter_Bragg_fwd_pi_pl1_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl1_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_ka_pl2_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pr_pl2_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_pi_pl2_rebuild[k][m] = -999; 
       reco_track_daughter_Bragg_fwd_mu_pl2_rebuild[k][m] = -999; 


       reco_track_daughter_distance[k][m]=-9;
       reco_track_daughter_vtx_distance[k][m]=-9;
       reco_angle_track_daughter[k][m]=-9;
       
       for(int l=0; l<kMaxMerge; l++){ 
	 reco_track_daughter_match_e[k][m][l] = -999;
	 reco_track_daughter_match_hit[k][m][l] = -999;
	 reco_track_daughter_match_epdg[k][m][l] = -999;
	 reco_track_daughter_match_hitpdg[k][m][l] = -999;
	 reco_track_daughter_match_e_rebuild[k][m][l] = -999;
	 reco_track_daughter_match_hit_rebuild[k][m][l] = -999;
	 reco_track_daughter_match_epdg_rebuild[k][m][l] = -999;
	 reco_track_daughter_match_hitpdg_rebuild[k][m][l] = -999;
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
       reco_track_daughter_mean_dedx_3pl[k][m] = -9;
       reco_track_daughter_mean_dedx_pl0[k][m] = -9;
       reco_track_daughter_mean_dedx_pl1[k][m] = -9;
       reco_track_daughter_mean_dedx_pl2[k][m] = -9;
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
       
       
     }
   }
   
   k_can_trkid=-999;
   mu_can_trkid=-999;
   k_mu_can_dis=-999;
   k_mu_open_angle=-999;
   k_vtx_dis=-999;
 }
 
 /////////////////////////////////////////////////////////////////////////////////////////
 
bool CCKaonAnalyzer::isInsideVolume(string volume, double x, double y, double z)
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

double CCKaonAnalyzer::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
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

  //void CCKaonAnalyzer::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, int track_i, int daughter_i, bool wTrackRebuilder)
  void CCKaonAnalyzer::fillCalorimetry(const std::vector<art::Ptr<anab::Calorimetry>> &calos, art::Ptr<recob::Track>& ptrack, int track_i, int daughter_i, bool wTrackRebuilder)
{

  double llr_pid_total = 0;
  double llr_pid_total_k = 0;

  for (auto const &calo : calos) {

    //int planenum = calo->PlaneID().Plane;

    auto const &plane = calo->PlaneID().Plane;
    auto const &dedx_values = calo->dEdx();

    /*
    if(planenum==0) dv0 = calo->dEdx();
    if(planenum==0) rv0 = calo->ResidualRange();

    if(planenum==1) dv1 = calo->dEdx();
    if(planenum==1) rv1 = calo->ResidualRange();

    if(planenum==2) dv2 = calo->dEdx();
    if(planenum==2) rv2 = calo->ResidualRange();
    */

    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();

    std::vector<std::vector<float>> par_values;
    par_values.push_back(rr);
    par_values.push_back(pitch);

    if (calo->ResidualRange().size() == 0) continue;

    float calo_energy = 0;
    for (size_t i = 0; i < dedx_values.size(); i++) {
      calo_energy += dedx_values[i] * pitch[i];
    }
    
    llr_pid_total += llr_pid_calculator.LLR_many_hits_one_plane(dedx_values, par_values, plane);
    llr_pid_total_k += llr_pid_calculator_k.LLR_many_hits_one_plane(dedx_values, par_values, plane);

  }

  double llr_pid_score = atan(llr_pid_total / 100.) * 2 / 3.14159266;
  double llr_pid_score_k = atan(llr_pid_total_k / 100.) * 2 / 3.14159266;


  //Mean dE/dX
  std::vector<std::pair<int,double>> MeandEdXs = MeandEdX(calos);
  double thisThreePlaneMeandEdX = ThreePlaneMeandEdX(ptrack,MeandEdXs);
  if( thisThreePlaneMeandEdX != thisThreePlaneMeandEdX ) std::cout << "NAN for three plane dedx!" << std::endl;
  //ThisPrimaryDaughter.MeandEdX_ThreePlane = thisThreePlaneMeandEdX;

  //add single plane mean dEdX


  if(wTrackRebuilder == true){

    if (daughter_i<0) {
    }
    else {

      for(size_t i_plane=0;i_plane<MeandEdXs.size();i_plane++){
	if( MeandEdXs.at(i_plane).first == 0 ) reco_track_daughter_mean_dedx_pl0_rebuild[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 1 ) reco_track_daughter_mean_dedx_pl1_rebuild[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 2 ) reco_track_daughter_mean_dedx_pl2_rebuild[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
      }

      reco_track_daughter_mean_dedx_3pl_rebuild[track_i][daughter_i] = ThreePlaneMeandEdX(ptrack,MeandEdXs);
      reco_track_daughter_llrpid_3pl_rebuild[track_i][daughter_i] = llr_pid_score;
      reco_track_daughter_llrpid_k_3pl_rebuild[track_i][daughter_i] = llr_pid_score_k;

    }

  }else{

    if (daughter_i<0) {

      for(size_t i_plane=0;i_plane<MeandEdXs.size();i_plane++){
	if( MeandEdXs.at(i_plane).first == 0 ) reco_track_mean_dedx_pl0[track_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 1 ) reco_track_mean_dedx_pl1[track_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 2 ) reco_track_mean_dedx_pl2[track_i] = MeandEdXs.at(i_plane).second;
      }

      reco_track_mean_dedx_3pl[track_i] = ThreePlaneMeandEdX(ptrack,MeandEdXs);
      reco_track_llrpid_3pl[track_i] = llr_pid_score;
      reco_track_total_llrpid_3pl[track_i] = llr_pid_total;
      reco_track_llrpid_k_3pl[track_i] = llr_pid_score_k;

    }else{

      for(size_t i_plane=0;i_plane<MeandEdXs.size();i_plane++){
	if( MeandEdXs.at(i_plane).first == 0 ) reco_track_daughter_mean_dedx_pl0[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 1 ) reco_track_daughter_mean_dedx_pl1[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
	if( MeandEdXs.at(i_plane).first == 2 ) reco_track_daughter_mean_dedx_pl2[track_i][daughter_i] = MeandEdXs.at(i_plane).second;
      }

      reco_track_daughter_mean_dedx_3pl[track_i][daughter_i] = ThreePlaneMeandEdX(ptrack,MeandEdXs);
      reco_track_daughter_llrpid_3pl[track_i][daughter_i] = llr_pid_score;
      reco_track_daughter_llrpid_k_3pl[track_i][daughter_i] = llr_pid_score_k;
    }

  }
  
}


void CCKaonAnalyzer::fillPID(const std::vector<art::Ptr<anab::ParticleID>> &trackPID, double angle_y, int track_i, int daughter_i, bool wTrackRebuilder)
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
  //Float_t No_Bragg[3] = {-1,-1,-1};
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
              
            for(int pl=0; pl<3; pl++){

            if (UBPID::uB_getSinglePlane(AlgScore.fPlaneMask)==pl) {
            if (TMath::Abs(AlgScore.fAssumedPdg)==2212) Bragg_fwd_pr[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==321) Bragg_fwd_ka[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==13) Bragg_fwd_mu[pl] = AlgScore.fValue;
            if (TMath::Abs(AlgScore.fAssumedPdg)==211) Bragg_fwd_pi[pl] = AlgScore.fValue;
            //if (TMath::Abs(AlgScore.fAssumedPdg)==0) No_Bragg[pl] = AlgScore.fValue;
            }

            }
          //}
      }
    }

  }

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

  if(wTrackRebuilder==true){
    if (daughter_i<0) {
    }
    else {

      reco_track_daughter_chi2ka_pl0_rebuild[track_i][daughter_i] = chi2ka[0];
      reco_track_daughter_chi2pr_pl0_rebuild[track_i][daughter_i] = chi2pr[0];
      reco_track_daughter_chi2pi_pl0_rebuild[track_i][daughter_i] = chi2pi[0];
      reco_track_daughter_chi2mu_pl0_rebuild[track_i][daughter_i] = chi2mu[0];
      reco_track_daughter_chi2ka_pl1_rebuild[track_i][daughter_i] = chi2ka[1];
      reco_track_daughter_chi2pr_pl1_rebuild[track_i][daughter_i] = chi2pr[1];
      reco_track_daughter_chi2pi_pl1_rebuild[track_i][daughter_i] = chi2pi[1];
      reco_track_daughter_chi2mu_pl1_rebuild[track_i][daughter_i] = chi2mu[1];
      reco_track_daughter_chi2ka_pl2_rebuild[track_i][daughter_i] = chi2ka[2];
      reco_track_daughter_chi2pr_pl2_rebuild[track_i][daughter_i] = chi2pr[2];
      reco_track_daughter_chi2pi_pl2_rebuild[track_i][daughter_i] = chi2pi[2];
      reco_track_daughter_chi2mu_pl2_rebuild[track_i][daughter_i] = chi2mu[2];
      reco_track_daughter_chi2ka_3pl_rebuild[track_i][daughter_i] = chi2ka_3pl;
      reco_track_daughter_chi2pr_3pl_rebuild[track_i][daughter_i] = chi2pr_3pl;
      reco_track_daughter_chi2pi_3pl_rebuild[track_i][daughter_i] = chi2pi_3pl;
      reco_track_daughter_chi2mu_3pl_rebuild[track_i][daughter_i] = chi2mu_3pl;
      reco_track_daughter_likepr_3pl_rebuild[track_i][daughter_i] = likepr_3pl;
      
      reco_track_daughter_Bragg_fwd_ka_pl0_rebuild[track_i][daughter_i] = Bragg_fwd_ka[0];
      reco_track_daughter_Bragg_fwd_pr_pl0_rebuild[track_i][daughter_i] = Bragg_fwd_pr[0];
      reco_track_daughter_Bragg_fwd_pi_pl0_rebuild[track_i][daughter_i] = Bragg_fwd_pi[0];
      reco_track_daughter_Bragg_fwd_mu_pl0_rebuild[track_i][daughter_i] = Bragg_fwd_mu[0];
      reco_track_daughter_Bragg_fwd_ka_pl1_rebuild[track_i][daughter_i] = Bragg_fwd_ka[1];
      reco_track_daughter_Bragg_fwd_pr_pl1_rebuild[track_i][daughter_i] = Bragg_fwd_pr[1];
      reco_track_daughter_Bragg_fwd_pi_pl1_rebuild[track_i][daughter_i] = Bragg_fwd_pi[1];
      reco_track_daughter_Bragg_fwd_mu_pl1_rebuild[track_i][daughter_i] = Bragg_fwd_mu[1];
      reco_track_daughter_Bragg_fwd_ka_pl2_rebuild[track_i][daughter_i] = Bragg_fwd_ka[2];
      reco_track_daughter_Bragg_fwd_pr_pl2_rebuild[track_i][daughter_i] = Bragg_fwd_pr[2];
      reco_track_daughter_Bragg_fwd_pi_pl2_rebuild[track_i][daughter_i] = Bragg_fwd_pi[2];
      reco_track_daughter_Bragg_fwd_mu_pl2_rebuild[track_i][daughter_i] = Bragg_fwd_mu[2];
      
    }
  }else{
    
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
      
    }
    
  }
  
}



void CCKaonAnalyzer::fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_recoobj,
                                      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
                                      int track_i,
                                      int daughter_i,
				      bool isTrack,
				      bool wTrackRebuilder)
{

  simb::MCParticle const* matched_mcparticle = NULL;
  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  
  for(size_t i_h=0; i_h<hits_from_recoobj.size(); i_h++) {
    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hits_from_recoobj[i_h].key(), particle_vec, match_vec);
    
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
    
    if(isTrack == true){

      if(wTrackRebuilder == true){
	if (daughter_i<0) {
	} else {
	  reco_track_daughter_true_pdg_rebuild[track_i][daughter_i] = matched_mcparticle->PdgCode();
	  reco_track_daughter_true_origin_rebuild[track_i][daughter_i] = 1;//int(mc_truth->Origin());
	  reco_track_daughter_true_primary_rebuild[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
	  reco_track_daughter_true_end_inTPC_rebuild[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_end_in5cmTPC_rebuild[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_end_inCCInclusiveTPC_rebuild[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_length_rebuild[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
	  
	}
      } else {
	
	if (daughter_i<0) {
	  reco_track_true_pdg[track_i] = matched_mcparticle->PdgCode();
	  reco_track_true_origin[track_i] = 1;//int(mc_truth->Origin());
	  reco_track_true_primary[track_i] = matched_mcparticle->Process()=="primary";
	  reco_track_true_end_inTPC[track_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_true_end_in5cmTPC[track_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_true_end_inCCInclusiveTPC[track_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_true_length[track_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
	} else {
	  reco_track_daughter_true_pdg[track_i][daughter_i] = matched_mcparticle->PdgCode();
	  reco_track_daughter_true_origin[track_i][daughter_i] = 1;//int(mc_truth->Origin());
	  reco_track_daughter_true_primary[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
	  reco_track_daughter_true_end_inTPC[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_end_in5cmTPC[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_end_inCCInclusiveTPC[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
	  reco_track_daughter_true_length[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
	}
      }
    } else {
      if(track_i>=0 && daughter_i>=0){	
	reco_track_daughter_true_pdg_sh[track_i][daughter_i] = matched_mcparticle->PdgCode();
	reco_track_daughter_true_origin_sh[track_i][daughter_i] = 1;//int(mc_truth->Origin());
	reco_track_daughter_true_primary_sh[track_i][daughter_i] = matched_mcparticle->Process()=="primary";
	reco_track_daughter_true_end_inTPC_sh[track_i][daughter_i] = isInsideVolume("TPC", matched_mcparticle->EndPosition().Vect());
	reco_track_daughter_true_end_in5cmTPC_sh[track_i][daughter_i] = isInsideVolume("5cmTPC", matched_mcparticle->EndPosition().Vect());
	reco_track_daughter_true_end_inCCInclusiveTPC_sh[track_i][daughter_i] = isInsideVolume("CCInclusiveTPC", matched_mcparticle->EndPosition().Vect());
	reco_track_daughter_true_length_sh[track_i][daughter_i] = length(*matched_mcparticle, mcstart, mcend, pstarti, pendi);
      }
    }
    
  } else {
    cout << "no matched_mcparticle found in fillTrue" << endl;
  }
  
}


void CCKaonAnalyzer::mergeChecker(std::vector<art::Ptr<recob::Hit>>& hits_from_recoobj,
				  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
				  int track_i,
				  int daughter_i,
				  bool isTrack,
				  bool wTrackRebuilder)
{

  std::map<int, double> trkide;
  std::map<int, double> trkidhit;
  std::map<int, int> trkidpdg;
  std::map<int, int> trkidmrgid;
  std::map<int, int> trkidmother;

  std::map<double, int, std::greater<int>> epdg;
  std::map<double, int, std::greater<int>> hitpdg;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;


  // Process each hit in the given collection
  for (auto const& hit : hits_from_recoobj) {
    particle_vec.clear();
    match_vec.clear();
    particles_per_hit.get(hit.key(), particle_vec, match_vec);

    // Accumulate information for each particle associated with the hit
    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
      int trackID = particle_vec[i_p]->TrackId();
      trkide[trackID] += match_vec[i_p]->energy;
      trkidhit[trackID]++;
      trkidpdg[trackID] = particle_vec[i_p]->PdgCode();
      trkidmrgid[trackID] = -9;
      trkidmother[trackID] = particle_vec[i_p]->Mother();
    }
  }

  // Generate merged IDs for related tracks
  int currentMergedID = 1;
  for (auto& pdg_entry : trkidpdg) {
    if (trkidmrgid[pdg_entry.first] != -9) continue;
    trkidmrgid[pdg_entry.first] = currentMergedID;
    int currentMotherTrackID = trkidmother[pdg_entry.first];

    //while (currentMotherTrackId > 0 && trkidpdg.find(currentMotherTrackId) != trkidpdg.end() && trkidpdg[currentMotherTrackId] == pdg_entry.second) {
    while (currentMotherTrackID > 0) {
      if( trkidpdg.find(currentMotherTrackID) != trkidpdg.end() &&
	  trkidpdg[currentMotherTrackID] == pdg_entry.second ){
	trkidmrgid[currentMotherTrackID] = currentMergedID;
	currentMotherTrackID = trkidmother[currentMotherTrackID];
      }
    }
    ++currentMergedID;
  }
  
  // Merge tracks with the same merged ID
  for (auto& mrgid1 : trkidmrgid) {
    for (auto& mrgid2 : trkidmrgid) {
      if (mrgid1.first == mrgid2.first || mrgid1.second != mrgid2.second) continue;
      trkide[mrgid1.first] += trkide[mrgid2.first];
      trkidhit[mrgid1.first] += trkidhit[mrgid2.first];
      trkidmrgid[mrgid2.first] = -1;  // Mark as merged
      //trkide[mrgid2.first] = -1;
      //trkidhit[mrgid2.first] = -1;
    }
  }

  // Prepare data for output
  for (auto const& ide : trkide) {
    if (trkidmrgid[ide.first] == -1) continue;
    epdg[ide.second] = trkidpdg[ide.first];
    hitpdg[trkidhit[ide.first]] = trkidpdg[ide.first];
  }

  // Output data

  auto assignEnergy = [&](double energy, int pdgCode, int index) {
    if (isTrack) {
      if (wTrackRebuilder) {
	if (daughter_i < 0) {
	} else {
	  reco_track_daughter_match_e_rebuild[track_i][daughter_i][index] = energy;
	  reco_track_daughter_match_epdg_rebuild[track_i][daughter_i][index] = pdgCode;
	}
      } else {
	if (daughter_i < 0) {
	  reco_track_match_e[track_i][index] = energy;
	  reco_track_match_epdg[track_i][index] = pdgCode;
	} else {
	  reco_track_daughter_match_e[track_i][daughter_i][index] = energy;
	  reco_track_daughter_match_epdg[track_i][daughter_i][index] = pdgCode;
	}
      }
    } else {
      reco_track_daughter_shower_match_e[track_i][daughter_i][index] = energy;
      reco_track_daughter_shower_match_epdg[track_i][daughter_i][index] = pdgCode;
    }
  };
  
  auto assignHits = [&](double hits, int pdgCode, int index) {
    if (isTrack) {
      if (wTrackRebuilder) {
	if (daughter_i < 0) {
	} else {
	  reco_track_daughter_match_hit_rebuild[track_i][daughter_i][index] = hits;
	  reco_track_daughter_match_hitpdg_rebuild[track_i][daughter_i][index] = pdgCode; 
	}
      } else {
	if (daughter_i < 0) {
	  reco_track_match_hit[track_i][index] = hits;
	  reco_track_match_hitpdg[track_i][index] = pdgCode;
	} else {
	  reco_track_daughter_match_hit[track_i][daughter_i][index] = hits;
	  reco_track_daughter_match_hitpdg[track_i][daughter_i][index] = pdgCode;
	}
      }
    } else {
      reco_track_daughter_shower_match_hit[track_i][daughter_i][index] = hits;
      reco_track_daughter_shower_match_hitpdg[track_i][daughter_i][index] = pdgCode;
    }
  };

  int merge_i = 0;
  for (auto const& x : epdg) {
    assignEnergy(x.first, x.second, merge_i++);
  }

  merge_i = 0; // Reset index for next loop
  for (auto const& y : hitpdg) {
    assignHits(y.first, y.second, merge_i++);
  }

}




///////////////////////////////////////////////////////////////////////////////////////////



  // DEFINE_ART_MODULE(CCKaonAnalyzer)
}
//}
