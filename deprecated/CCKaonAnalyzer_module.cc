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
  isMC                             (pset.get< bool >("IsMC",true)),
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
  
  //begin by resetting everything
  reset();
  
  /////////////////////////
  // EVENT ID INFORMATION /
  /////////////////////////

  fEventID = evt.id().event();
  run = evt.run();
  subrun = evt.subRun();
  event = evt.event();


  //////////////////////////////
  // GET EVENT GENERATOR INFO //
  //////////////////////////////
  
  if(isMC) {
    
    art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mcTrVect;
    if(e.getByLabel(fGenieGenModuleLabel,mctruthListHandle)) art::fill_ptr_vector(mcTrVect,mctruthListHandle); 

    if(!mcTrVect.size()) return;
    
    fNMCTruths = mcTrVect.size();
    art::Ptr<simb::MCTruth> MCtruth = mcTrVect.at(0);
    
    // true neutrino information
    simb::MCNeutrino Nu = MCtruth->GetNeutrino();
    mode = Nu.Mode();
    ccnc = Nu.CCNC();
    
    if(ccnc == 0) fCCNC = "CC";
    else fCCNC = "NC";
    
    if(mode == 0) fMode = "QEL";
    else if(mode == 1) fMode = "RES";
    else if(mode == 2) fMode = "DIS";
    else if(mode == 3) fMode = "COH";
    else if(mode == 5) fMode = "ElectronScattering";
    else if(mode == 10) fMode = "MEC";
    else if(mode == 11) fMode = "Diffractive";
    else fMode = "Other";

    for(int k_particles=0;k_particles<MCtruth->NParticles();k_particles++){
      
      simb::MCParticle Part = MCtruth->GetParticle(k_particles);
      
      // Get list of particles from true PV, if lepton or neutrino set PV
      if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1){ 
	fTruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());
	fInActiveTPC=inActiveTPC(fTruePrimaryVertex);
      }//if lepton
      
      // Record info about the neutrino
      if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
	SimParticle P;
	P.SetKinematics(Part.Momentum(),Part.Mass());
	P.SetPositions(Part.Position(),Part.EndPosition());
	P.PDG = Part.PdgCode();
	P.Origin = 0;
	fNeutrino.push_back( P );
	
	//true_nu_vtx_inTPC = isInsideVolume("TPC",mctruth.GetNeutrino().Nu().Position().Vect());
	//true_nu_vtx_in5cmTPC = isInsideVolume("5cmTPC",mctruth.GetNeutrino().Nu().Position().Vect());
	//true_nu_vtx_inCCInclusiveTPC = isInsideVolume("CCInclusiveTPC",mctruth.GetNeutrino().Nu().Position().Vect());

      }      
    }
  }


  ////////////////////////////////////////////
  // Get Geant Information
  ///////////////////////////////////////////
  
  if(isMC) {
    
    //get list of geant4 particles
    art::Handle<std::vector<simb::MCParticle>>g4particleHandle;
    std::vector< art::Ptr<simb::MCParticle>>g4particleVect;
    g4particleVect.clear();
    
    if(e.getByLabel(fLArG4ModuleLabel,g4particleHandle)) art::fill_ptr_vector(g4particleVect,g4particleHandle);
    else
      std::cout << "Geant Particles Missing!" << std::endl;
    
    
    //id numbers of particles produced at primary vertex
    //KaonPlus_IDs.clear();
    primary_IDs.clear();
    
    //id numbers of primary particle daughters
    daughter_IDs.clear();
    KaonPlus_daughter_IDs.clear();
    KaonMinus_daughter_IDs.clear();
    
    //map between particle ID numbers (g4p->TrackId()) and pointers to simb::MCParticle
    partByID.clear();
    
    
    for(const art::Ptr<simb::MCParticle> &g4p : g4particleVect){
      
      //if mother == 0 particle is produced at primary vertex
      if(g4p->Mother() == 0){ 
	
	//store vector of id's of primary particles
	primary_IDs.push_back(g4p->TrackId());
	
	if(isKaon(g4p->PdgCode())){ //if the particle is a kaon, store daughter ids
	  
	  fIsKaon = true;
	  
	  if(g4p->PdgCode() == 321){//K+
	    
	    fIsKaonPlus = true;
	    
	    if(g4p->EndProcess() == "Decay"){
	      for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
		KaonPlus_daughter_IDs.push_back( g4p->Daughter(i_d) );
	      }
	    }
	    else if(g4p->EndProcess() == "kaon+Inelastic"){
	      for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
		KaonPlus_Inelastic_daughter_IDs.push_back( g4p->Daughter(i_d) );
	      }
	    }	 
	    
	    // record decay vertex of K+
	    fDecayVertex.SetXYZ( g4p->EndPosition().X() , g4p->EndPosition().Y() , g4p->EndPosition().Z() );
	    
	  } else if(g4p->PdgCode() == -321){//K-
	    
	    fIsKaonMinus = true;
	    
	    if(g4p->EndProcess() == "Decay"){
	      for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
		KaonMinus_daughter_IDs.push_back( g4p->Daughter(i_d) );
	      }	      
	    }
	    
	  } else {//Other Ks
	    if(g4p->EndProcess() == "Decay"){ 
	      for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
		KaonOthers_daughter_IDs.push_back( g4p->Daughter(i_d) );
	      }
	    }
	  }
	  
	} else if(isHyperon(g4p->PdgCode())) { // hyperon from associated production
	  
	  fIsHyperon = true;
	  
	  if(g4p->EndProcess() == "Decay"){ 
	    for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
	      Hyperon_daughter_IDs.push_back( g4p->Daughter(i_d) );
	    }
	  }
	  
	}
	
      }
      // fill the map of MCParticle and TrackID
      part_and_id = std::make_pair(g4p->TrackId() , g4p);
      partByID.insert( part_and_id );
      
    }
    
    
    for(size_t i_p=0;i_p<primary_IDs.size();i_p++){
      
      //geant does not always keep everything it simulates, make sure particle is in list of IDs (crashes otherwise!)
      if(partByID.find(primary_IDs[i_p]) == partByID.end()) continue;
      art::Ptr<simb::MCParticle> part = partByID[primary_IDs.at(i_p)];
      
      if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
      
      SimParticle P = MakeSimParticle(*part);
      P.Origin = getOrigin(part->TrackId());
      
      if( isKaon(part->PdgCode()) ) { 
	//K+ from primary vertex
	if( isKaonPlus(part->PdgCode()) ) fKaonPlus.push_back( P );
	
	//K- from primary vertex
	else if( isKaonMinus(part->PdgCode()) ) fKaonMinus.push_back( P );
	
	//Other K from primary vertex
	else fKaonOthers.push_back( P );
      }
      
      //lepton produced at primary vertex
      if( isLepton(part->PdgCode()) || isNeutrino(part->PdgCode()) ) fLepton.push_back( P );
      
      //hyperon produced at primary vertex
      if( isHyperon(part->PdgCode()) ) fPrimaryHyperon.push_back( P );
      
      //nucleon produced at primary vertex
      if( isNucleon(part->PdgCode()) ) fPrimaryNucleon.push_back( P );
      
      //pion produced at primary vertex
      if( isPion(part->PdgCode()) ) fPrimaryPion.push_back( P );
      
    }
    
    
    //check all decay products are actually produced at decay vertex - sometimes you get some electrons thrown in
    if(fKaonPlus.size() == 1){
      
      //if using GENIE as evgen, hyperon events get labelled as QEL - change this to HYP
      if( fMode == "QEL" || fMode == "RES" || fMode == "DIS" ) fMode = "KAON";
      //if( fMode == "RES" || fMode == "DIS" ) fMode = "KAON";
      
      std::vector<int> KaonPlus_daughter_IDs_tmp;
      
      for(size_t i_d=0;i_d<KaonPlus_daughter_IDs.size();i_d++){
	
	//geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
	if(partByID.find(KaonPlus_daughter_IDs[i_d]) == partByID.end()) continue;
	
	art::Ptr<simb::MCParticle> part = partByID[KaonPlus_daughter_IDs[i_d]];
	
	if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
	
	if(part->Position().X() == fKaonPlus.at(0).EndX && part->Position().Y() == fKaonPlus.at(0).EndY && part->Position().Z() == fKaonPlus.at(0).EndZ){
	  
	  KaonPlus_daughter_IDs_tmp.push_back(KaonPlus_daughter_IDs.at(i_d));
	  
	}
      }
      
      //refined id list of K+
      KaonPlus_daughter_IDs = KaonPlus_daughter_IDs_tmp;
      
    }
    
    if(fKaonMinus.size() == 1){
      
      std::vector<int> KaonMinus_daughter_IDs_tmp;
      
      for(size_t i_d=0;i_d<KaonMinus_daughter_IDs.size();i_d++){
	
	//geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
	if(partByID.find(KaonMinus_daughter_IDs[i_d]) == partByID.end()) continue;
	
	art::Ptr<simb::MCParticle> part = partByID[KaonMinus_daughter_IDs[i_d]];
	
	if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
	
	if(part->Position().X() == fKaonMinus.at(0).EndX && part->Position().Y() == fKaonMinus.at(0).EndY && part->Position().Z() == fKaonMinus.at(0).EndZ){
	  
	  KaonMinus_daughter_IDs_tmp.push_back(KaonMinus_daughter_IDs.at(i_d));
	  
	}
      }
      
      KaonMinus_daughter_IDs = KaonMinus_daughter_IDs_tmp;
      
    }
    
    
    if(fHyperon.size() == 1){
      
      std::vector<int> Hyperon_daughter_IDs_tmp;
      
      for(size_t i_d=0;i_d<Hyperon_daughter_IDs.size();i_d++){
	
	//geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
	if(partByID.find(Hyperon_daughter_IDs[i_d]) == partByID.end()) continue;
	
	art::Ptr<simb::MCParticle> part = partByID[Hyperon_daughter_IDs[i_d]];
	
	if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
	
	if(part->Position().X() == fPrimaryHyperon.at(0).EndX && part->Position().Y() == fPrimaryHyperon.at(0).EndY && part->Position().Z() == fPrimaryHyperon.at(0).EndZ){
	  
	  Hyperon_daughter_IDs_tmp.push_back(Hyperon_daughter_IDs.at(i_d));
	  
	}
      }
      
      Hyperon_daughter_IDs = Hyperon_daughter_IDs_tmp;
      
    }
    
    
    //now go through list of kaon daughters, get info about the decay
    for(size_t i_d=0;i_d<KaonPlus_daughter_IDs.size();i_d++){
      
      //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(KaonPlus_daughter_IDs[i_d]) == partByID.end()) continue;
      
      art::Ptr<simb::MCParticle> part = partByID[KaonPlus_daughter_IDs[i_d]];
      
      if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
      
      SimParticle Decay = MakeSimParticle(*part);
      Decay.Origin = getOrigin(part->TrackId());
      fDecay.push_back( Decay );
      
    }
    
    
    
    //we want to check K+ -> Nu_mu + mu+, K+ -> Pi+ + Pi0
    
    if(fIsKaonPlus){
      
      int countMuP = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == -13; });
      int countEP = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == -11; });
      int countNuMu = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == 14; });
      int countNuE = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == 12; });
      int countPiN = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == 111; });
      int countPiP = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == 211; });
      int countPiM = std::count_if(fDecay.begin(), fDecay.end(), [](const SimParticle& p) { return p.PDG == -211; });
      
      if(fDecay.size() == 2){
	
	fDecayOpeningAngle = (180/3.1416)*TMath::ACos( (fDecay.at(0).Px*fDecay.at(1).Px + fDecay.at(0).Py*fDecay.at(1).Py + fDecay.at(0).Pz*fDecay.at(1).Pz)/(fDecay.at(0).ModMomentum*fDecay.at(1).ModMomentum));
	
	if( countNuMu==1 && countMuP==1 ) fIsKaonPlus_NuMuP = true;
	else if( countPiN==1 && countPiP==1 ) fIsKaonPlus_PiPPiN = true;
	
      }else if(fDecay.size() == 3){
	
	if( countPiP==2 && countPiM==1 ) fIsKaonPlus_2PiPPiM = true;
	else if( countPiN==1 && countEP==1 && countNuE==1 ) fIsKaonPlus_ENuE = true;
	else if( countPiN==2 && countPiP==1 ) fIsKaonPlus_2PiNPiP = true;
	else fIsKaonPlus_Others = true;
	
      }
    }
    
    if(!KaonPlus_Inelastic_daughter_IDs.empty()){
      if(fIsKaonPlus) fIsInelastic_KaonPlus = true;
      else fIsInelastic_Others = true;
    }
    
    
    // Add Hyperon daughters
    for(size_t i_d=0;i_d<Hyperon_daughter_IDs.size();i_d++){
      
      //geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Hyperon_daughter_IDs[i_d]) == partByID.end()) continue;
      
      art::Ptr<simb::MCParticle> part = partByID[Hyperon_daughter_IDs[i_d]];
      
      if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these
      
      SimParticle HyperonDecay = MakeSimParticle(*part);
      HyperonDecay.Origin = getOrigin(part->TrackId());
      fHyperonDecay.push_back( KaonDecay );
      
    }
    
  }// if isMC
  
  
  //FV cut to daughter tracks?
  if(fNeutrino.size() == 1 && fInActiveTPC && fIsKaonPlus && fNeutrino.at(0).PDG == 14 && ( fMode == "QEL" || fMode == "KAON")){ 
    
    //add kinematic thresholds?
    if(fIsKaonPlus_NuMuP = true) fIsSignal_NuMuP = true;
    else if(fIsKaonPlus_PiPPiN = true) fIsSignal_PiPPiN = true;
    else{
      fIsSignal_NuMuP = false;
      fIsSignal_PiPPiN = false;
    }
  }
  
  
  
  ///////////////////////////////////////////////////////////////////////////
  //Get Reconstructed Info
  //////////////////////////////////////////////////////////////////////////
  
  
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
  
  
  //setup handles
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  art::Handle< std::vector<recob::Hit> > hitHandle;
  art::Handle< std::vector<recob::Track> > trackHandle;
  art::Handle< std::vector<recob::Track> > trackRebuiltHandle;
  art::Handle< std::vector<recob::Shower> > showerHandle;
  art::Handle< std::vector<anab::ParticleID> > pidHandle;
  art::Handle< std::vector<anab::Calorimetry> > caloHandle;
  art::Handle< larpandoraobj::PFParticleMetadata > particleMetadataHandle;
  art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
  
  std::vector< art::Ptr<recob::PFParticle> >pfparticleVect;
  std::vector<art::Ptr<recob::Hit> > hitVect;
  std::vector<art::Ptr<recob::Track> > trackVect;
  std::vector<art::Ptr<recob::Track> > trackRebuiltVect;
  std::vector<art::Ptr<recob::Shower> > showerVect;
  std::vector < art::Ptr<anab::ParticleID> > pidVect;
  std::vector<art::Ptr<anab::Calorimetry> > caloVect;
  std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > metadataVect;
  std::vector<art::Ptr<recob::SpacePoint> > spacepointVect;
  
  if(evt.getByLabel(fPFParticleLabel,pfparticleHandle)){
    art::fill_ptr_vector(pfparticleVect,pfparticleHandle);
  }
  
  if(!pfparticleVect.size()) return;
  
  
  if(evt.getByLabel(fTrackModuleLabel,trackHandle)) {
    art::fill_ptr_vector(trackVect, trackHandle);
  }
  else std::cout << "Track handle not setup" << std::endl;
  
  if(evt.getByLabel(fRebuiltTrackModuleLabel,rebuilttrackHandle)) {
    art::fill_ptr_vector(rebuilttrackVect, rebuilttrackHandle);
  }
  else std::cout << "Rebuild Track handle not setup" << std::endl;
  
  if(evt.getByLabel(fShowerModuleLabel,showerHandle)) {
    art::fill_ptr_vector(showerVect, showerHandle);
  }
  else std::cout << "Shower handle not setup" << std::endl;
  
  if(evt.getByLabel(fHitsModuleLabel,hitHandle)){
    art::fill_ptr_vector(hitVect, hitHandle);
  }
  else std::cout << "Hit handle not setup" << std::endl;
  
  if(evt.getByLabel(fSpacePointproducer,spacepointHandle)){
    art::fill_ptr_vector(spacepointVect, spacepointHandle);
  }
  else std::cout << "SpacePoint handle not setup" << std::endl;
  
  
  //vertices, tracks and showers assoc with PFPs
  art::FindManyP<anab::T0> pfp_muon_assn(pfparticlesHandle, evt, "NuCCproducer");
  art::FindManyP<recob::Vertex> vertexAssoc(pfparticleVect,evt,fVertexLabel);
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect,evt,fTrackModuleLabel);
  art::FindManyP<recob::Track> trackRebuiltAssoc(pfparticleVect,evt,fTrackRebuiltModuleLabel);
  art::FindManyP<recob::Shower> showerAssoc(pfparticleVect,evt,fShowerModuleLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssoc(pfparticleVect,evt,fMetadataLabel);
  
  art::FindManyP<anab::ParticleID> PIDAssoc(trackVect,evt,fPIDLabel);
  art::FindManyP<anab::ParticleID> PIDRebuiltAssoc(trackRebuiltVect,evt,fPIDRebuiltLabel);
  art::FindManyP<anab::Calorimetry> caloTrackAssoc(trackVect,evt,fCaloLabel);
  art::FindManyP<anab::Calorimetry> caloTrackRebuiltAssoc(trackRebuiltVect,evt,fCaloRebuiltLabel);
  
  //get hits assoc with tracks
  art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,evt,fTrackHitAssnLabel);
  //get hits assoc with showers
  art::FindManyP<recob::Hit> showerHitAssoc(showerHandle,evt,fShowerHitAssnLabel);
  
  //backtracker
  art::FindMany< simb::MCParticle , anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle,evt,fHitTruthAssnLabel);
  
  //spacepoints assoc with PFPs
  art::FindManyP<recob::SpacePoint> pfpSpacePointAssoc(pfparticleHandle,e,fPFPSpacePointAssnLabel);
  
  size_t neutrinoID = 99999;
  
  
  //go through the list of pandora PFP's
  lar_pandora::PFParticleVector pfneutrinos(0);
  lar_pandora::PFParticleVector pfmuons(0);
  
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
    
    std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());
    
    //get reconstructed neutrino
    if( pfp->IsPrimary() && std::abs(pfp->PdgCode())==14 ){
      
      neutrinoID = pfp->Self();
      fNPrimaryDaughters = pfp->NumDaughters();     
      pfneutrinos.push_back(pfp);
      
      //get the reconstructed primary vertex
      for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){
	
	//correct for space charge
	geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };                
	geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
	
	//w SC correction - forward
	fRecoPrimaryVertex.SetXYZ( vtx->position().X() + sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z());
	//WE NEED TO STORE RECO NEUTRINO INFORMATION!!
	
      }
    }
  }
  
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  
  // Find CC muon and daughters
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){ 
    
    // look at particles with neutrino parent and one associated track
    if ( pfp->Parent()==pfnu->Self() && trackAssoc.at(pfp.key()).size()==1) { 
      
      art::Ptr<recob::Track> track = trackAssoc.at(pfp.key()).front();
      // CC muon has a T0 associated 
      if (pfp_muon_assn.at(pfp.key()).size()==1) 
	pfmuons.push_back(pfparticle);
      
    }
  }
  
  if(pfmuons.size()==1) pfmuon = pfmuons.front();
  art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();

  
  
  //go through rest of particles, get lots of useful info!
  //we need primary tracks and their daughter tracks
  
  std::vector<TVector3> TrackStarts;
  
  for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){ // usual reconstruction
    
    RecoParticle ThisPrimary;
    RecoParticle ThisPrimaryDaughter;
    RecoParticle ThisPrimaryDaughterRebuilt;
    
    // decalre versatile recoparticle and fill at the end?
    
    //get data from every PFP, not just neutrino daughters
    //if(pfp->Parent() != neutrinoID) continue;
    //if(pfp->Parent() != neutrinoID) {
    
    std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());
    std::vector< art::Ptr<recob::Vertex> > pfpVertex = vertexAssoc.at(pfp.key());
    
    std::vector< art::Ptr<recob::Shower> > pfpShowers = showerAssoc.at(pfp.key());
    std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> >pfpMeta = metadataAssoc.at(pfp.key());
    
    std::vector< art::Ptr<recob::SpacePoint> > pfpSpacePoints = pfpSpacePointAssoc.at(pfp.key());
    
    if(pfp->Parent() == neutrinoID) ThisPrimary.PDG = pfp->PdgCode();
    else if(pfp->Parent() != neutrinoID) ThisPrimaryDaughter.PDG = pfp->PdgCode();
    
    for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){
      
      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());
      
      if (!pfParticlePropertiesMap.empty()){
	for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){
	  
	  if(it->first == "TrackScore"){
	    ThisPrimaryDaughter.TrackShowerScore = it->second;
	  }
	}
      }
    }
    
    
    
    //////////////////////////////////////////////
    //            Get track info                //
    //////////////////////////////////////////////


    //if pfp is a track like object
    if(pfpTracks.size() == 1 && !pfpVertex.empty()){
      
      for(const art::Ptr<recob::Track> &trk : pfpTracks){
	
	//sets track length/position related variables in ThisPrimaryDaughter
	SetTrackVariables(ThisPrimaryDaughter , trk);
	TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
	
	if(isMC) ThisPrimaryDaughter.HasTruth = true;
	
	if(isMC){
	  
	  //get hits assoc with track
	  std::vector< art::Ptr< recob::Hit> > hits = trackHitAssoc.at(trk.key());
	  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle, evt, fHitTruthAssns);
	  
	  simb::MCParticle const* matchedParticle = NULL;
	  matchedParticle = fillTrueMatching(hits, particles_per_hit);
	  
	  
	  if(matchedParticle != NULL){
	    
	    SimParticle P = MakeSimParticle(*matchedParticle);
	    P.Origin = getOrigin(matchedParticle->TrackId());
	    
	    ThisPrimaryDaughter.HasTruth = true;
	    ThisPrimaryDaughter.TrackTruthPurity = maxe/tote;
	    
	    ThisPrimaryDaughter.TrackTruePDG = P.PDG;
	    ThisPrimaryDaughter.TrackTrueE = P.E;
	    ThisPrimaryDaughter.TrackTruePx = P.Px;
	    ThisPrimaryDaughter.TrackTruePy = P.Py;
	    ThisPrimaryDaughter.TrackTruePz = P.Pz;
	    
	    ThisPrimaryDaughter.TrackTrueModMomentum = P.ModMomentum;
	    ThisPrimaryDaughter.TrackTrueKE = P.KE;
	    
	    ThisPrimaryDaughter.TrackTrueLength = P.Travel;
	    
	    ThisPrimaryDaughter.TrackTrueOrigin = P.Origin;
	    
	  }
	  else ThisPrimaryDaughter.HasTruth = false;
	  
	}//if isMC
	
	
	
	///////////////////////
	// get PID for track //
	///////////////////////
	
	// Setup Calo Assn
	std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());
	std::vector<art::Ptr<anab::Calorimetry>> PIDFromTrack = PIDAssoc.at(trk.key());
	
	fillCalorimetry(caloFromTrack, trk, ThisPrimaryDaughter);//fill values
	fillPID(PIDFromTrack, trk, ThisPrimaryDaughter);//fill values
	
      }//end of pfptracks
      
      
      ////////////////////////////
      // Get Vertex information //
      ////////////////////////////
      
      for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){
	
	geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };                
	geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
	
	//w SC correction - forward
	TVector3 pos( vtx->position().X() + sce_corr.X() , vtx->position().Y() - sce_corr.Y() , vtx->position().Z() - sce_corr.Z() );
	
	ThisPrimaryDaughter.SetVertex( pos );
	ThisPrimaryDaughter.Displacement = (pos - fRecoPrimaryVertex).Mag();
	
      }
    }// end of pfpTrack loop
    
    if(ThisPrimaryDaughter.PDG == -13 && pfpTracks.size() == 1){
      fTrackPrimaryDaughters_NuMuP.push_back( ThisPrimaryDaughter );
      fNPrimaryTrackDaughters_NuMuP++;
    }
    else if(ThisPrimaryDaughter.PDG == 211 && pfpTracks.size() == 1){
      fTrackPrimaryDaughters_PiPPiN.push_back( ThisPrimaryDaughter );
      fNPrimaryTrackDaughters_PiPPiN++;
    }
    
    
  }//end of pfp loop


  //set indices in particle vectors
  for(size_t i_tr=0;i_tr<fTrackPrimaryDaughters.size();i_tr++) fTrackPrimaryDaughters[i_tr].Index = i_tr;

  //store truth matching info for muon, decay proton and pion
  StoreTrackTruth();
  
  if(fPrint) PrintInfo();
  
  FinishEvent();
  
}



  ///////////////////////////////////////////
  // Print some useful info from the event //
  ///////////////////////////////////////////

  void CCKaonAnalyzer::printInfo(){

    if( fKaonPlus.size() != 1 ) return;

    if( !(fKaonPlus.at(0).PDG == fPrintPdg || fKaonPlus == -1) || !fInActiveTPC) return;

    std::cout << std::endl;
    std::cout << "EventID: " << fEventID-1 << std::endl;

    std::cout << "Truth info" << std::endl;

    std::cout << "Kaon: " << std::endl;
    for(size_t i=0;i<fKaon.size();i++) fKaon.at(i).Print();

    std::cout << "Lepton: " << std::endl;
    for(size_t i=0;i<fLepton.size();i++) fLepton.at(i).Print();

    std::cout << "Kaon decay products: " << std::endl;
    for(size_t i=0;i<fDecay.size();i++) fDecay.at(i).Print();

    if(fDecay.size() == 2){
      if(fIsKaonPlus_NuMuP) std::cout << "This is K+ -> Nu + Mu+ decay" << std::endl;
      if(fIsKaonPlus_PiPPiN) std::cout << "This is K+ -> Pi+ + Pi0 decay" << std::endl;
      std::cout << "Decay Opening Angle:  " << fDecayOpeningAngle << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Reco info" << std::endl;
    if(fIsKaonPlus_NuMuP) {

      std::cout << "Num tracklike daughters of NuMuP: " << fNPrimaryTrackDaughters_NuMuP << std::endl;
      std::cout << "List tracklike of daughters of NuMuP: " << std::endl;
      
      for(size_t i=0;i<fTrackPrimaryDaughters_NuMuP.size();i++){
	fTrackPrimaryDaughters_NuMuP.at(i).Print();
      }

    }
    else if(fIsKaonPlus_PiPPiN) {

      std::cout << "Num tracklike daughters of PiPPiN: " << fNPrimaryTrackDaughters_PiPPiN << std::endl;
      std::cout << "List tracklike of daughters of PiPPiN: " << std::endl;

      for(size_t i=0;i<fTrackPrimaryDaughters_PiPPiN.size();i++){
	fTrackPrimaryDaughters_PiPPiN.at(i).Print();
      }

    }

    std::cout << std::endl;

  }


  /////////////////////////////////////////////
  // Check origin of particle by its TrackId //
  /////////////////////////////////////////////

  int CCKaonAnalyzer::getOrigin(int idnum){

    //search list of primaries
    for(size_t i_p=0;i_p<primary_IDs.size();i_p++){
      if(primary_IDs.at(i_p) == idnum){ 
	return 1;
      }
    }

    //search list of K+ decay products
    for(size_t i_d=0;i_d<KaonPlus_daughter_IDs.size();i_d++){
      if(KaonPlus_daughter_IDs.at(i_d) == idnum){ 
	return 2;
      }
    }

    //search list of K- decay products
    for(size_t i_d=0;i_d<KaonMinus_daughter_IDs.size();i_d++){
      if(KaonMinus_daughter_IDs.at(i_d) == idnum){ 
	return 3;
      }
    }

    //search list of hyperon decay products
    for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){
      if(Hyperon_daughter_IDs.at(i_d) == idnum){ 
	return 4;
      }
    }


    return 5;

  }



  ///////////////////////////////////////////////////////////////
  // Finished processing event - update Metadata and fill tree //
  ///////////////////////////////////////////////////////////////

  void CCKaonAnalyzer::FinishEvent(){

    if(fDebug) std::cout << "Finishing Event" << std::endl;

    if(!fIsData){

      if(fCCNC == "CC") fNChargedCurrent++;
      else fNNeutralCurrent++;

      if( fNeutrino.size() != 1 ) std::cout << "Number of simulated neutrinos in this event != 1 !!" << std::endl;
      else if(fNeutrino.at(0).PDG == 12) fNnue++;
      else if(fNeutrino.at(0).PDG == 14) fNnuMu++;
      else if(fNeutrino.at(0).PDG == -12) fNnueBar++;
      else if(fNeutrino.at(0).PDG == -14) fNnuMuBar++;


      //genie uses QEL for hyperon events, NuWro uses HYP
      if(fNeutrino.size() == 1 && fInActiveTPC && fIsKaonPlus && fNeutrino.at(0).PDG == -14 && ( fMode == "QEL" || fMode == "KAON")){ 

	if(fIsKaonPlus_NuMuP = true) fIsSignal_NuMuP = true;
	else if(fIsKaonPlus_PiPPiN = true) fIsSignal_PiPPiN = true;
	else
	  fIsKaonPlus_PiPPiN = false, fIsSignal_PiPPiN = false;

      } 

      //asscociated hyperon/kaon production tagger - hyperon and kaon in the final state
      if( fKaonPlus.size() && fPrimaryHyperon.size() ) fIsAssociatedKaonPlus = true;
      if( fKaonMinus.size() && fPrimaryHyperon.size() ) fIsAssociatedKaonMinus = true;
      if( fKaonOthers.size() && fPrimaryHyperon.size() ) fIsAssociatedKaonOthers = true;

      //if is a signal event, check if decay products were reconstructed
      if( fIsSignal_NuMuP && fTrueDecayMuonIndex >= 0 ) fGoodReco_NuMuP = true;
      if( fIsSignal_PiPPiN && fTrueDecayPionIndex >= 0 ) fGoodReco_PiPPiN = true;

    }

    //Perform selection
    fSelectedEvent = PerformSelection();

    //store info for this event
    fOutputTree->Fill();


    //update metadata

    fNEvents++; //total events in sample
    if(fIsKaonPlus && fInActiveTPC) fNKaonPlus++; //total K+ events in active vol
    if(fSelectedEvent) fNSelectedEvents++; //total events passing selection
    if(fSelectedEvent && fInActiveTPC && fIsKaonPlus) fNSelectedKaonPlus++; //total hyperons in active vol passing preselection

    //signal events
    if(fIsSignal_NuMuP) fNSignal_NuMuP++;
    if(fIsSignal_PiPPiN) fNSignal_PiPPiN++;
    if(fIsSignal_NuMuP && fSelectedEvent) fNSelectedSignal_NuMuP++;
    if(fIsSignal_PiPPiN && fSelectedEvent) fNSelectedSignal_PiPPiN++;

    //signal events
    if(fGoodReco) fNGoodReco++;
    if(fGoodReco && fSelectedEvent) fNSelectedGoodReco++;

    
  }


  void CCKaonAnalyzer::StoreTrackTruth(){

    fTrueMuonIndex=-1;
    fTrueKaonIndex=-1;
    fTrueDecayMuonIndex=-1;
    fTrueDecayPionIndex=-1;

    //can be multiple tracks corresponding to the muon/proton/pion
    //first get all the tracks truth matching to primary muon, decay proton, decay pion

    std::vector<int> CCMuons;
    std::vector<int> Kaons;
    std::vector<int> Protons;
    std::vector<int> DecayMuons;
    std::vector<int> DecayPions;

    //add primary track

    for(size_t i_tr=0;i_tr<fTrackPrimaryDaughters.size();i_tr++){

      //if track does not have a truth matching, skip
      if(!fTrackPrimaryDaughters.at(i_tr).HasTruth) continue;

      //if muon produced at primary vertex
      if( abs(fTrackPrimaryDaughters.at(i_tr).TrackTruePDG) == 13 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 1 ) CCMuons.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

      //if proton produced at primary vertex
      if( fTrackPrimaryDaughters.at(i_tr).TrackTruePDG == 2212 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 1 ) Protons.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

      //if mu+ produced from kaon decay
      if( fTrackPrimaryDaughters.at(i_tr).TrackTruePDG == -13 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 2 ) DecayMuons.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

      //if pi+ produced from kaon decay
      if( fTrackPrimaryDaughters.at(i_tr).TrackTruePDG == -211 && fTrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 2 ) DecayPions.push_back( fTrackPrimaryDaughters.at(i_tr).Index );

    }

    //if there are no muons, protons or pions, exit here
    if( !CCMuons.size() && !DecayMuons.size() && !DecayPions.size() ) return;

    //set muon information
    //if multiple muons found, choose the one closest to reco'd primary vertex
    double min_dist=10000;
    for(size_t i_m=0;i_m<CCMuons.size();i_m++){

      TVector3 MuonStart(fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartX,fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartY,fTrackPrimaryDaughters.at(Muons.at(i_m)).TrackStartZ);

      double d = (MuonStart - fRecoPrimaryVertex).Mag();
      if(d < min_dist) { fTrueMuonIndex = Muons.at(i_m); min_dist = d;  }

    }


    //if there are no decay products, exit here
    if( !DecayMuons.size() && !DecayPions.size() ) return;

    if( DecayMuons.size() ) {

      double min_dist=10000; 
      for(size_t i_mu=0;i_mu<DecayMuons.size();i_mu++){ 

	TVector3 DecayMuonStart(fTrackPrimaryDaughters.at(DecayMuons.at(i_m)).TrackStartX,fTrackPrimaryDaughters.at(DecayMuons.at(i_m)).TrackStartY,fTrackPrimaryDaughters.at(DecayMuons.at(i_m)).TrackStartZ);
	
	double d = (DecayMuonStart - fRecoPrimaryKaonEnd).Mag();
	if(d < min_dist) { fTrueMuonIndex = DecayMuons.at(i_m); min_dist = d; }
      
      }

    }
    if( DecayPions.size() ) {
      
      double min_dist=10000;
      for(size_t i_pi=0;i_pi<DecayPions.size();i_pi++){

	TVector3 DecayPionStart(fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartX,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartY,fTrackPrimaryDaughters.at(DecayPions.at(i_pi)).TrackStartZ);
	double d = (DecayPionStart - fRecoPrimaryKaonEnd).Mag();
	if(d < min_dist) { fTruePionIndex = DecayPions.at(i_m); min_dist = d; }

      }

    }

  }



} // end of analyze function
 

 /////////////////////////////////////////// Reset Function ///////////////////////////////
 
 void CCKaonAnalyzer::reset(){

   //begin by resetting everything

   /////////////
   // General //
   /////////////

   fNMCTruths=0;

   fInActiveTPC=false;
   fIsHyperon=false;
   fIsSigmaZero=false;
   fIsLambda=false;
   fIsLambdaCharged=false;
   fIsAssociatedHyperon=false;
   fIsSignal=false;
   fGoodReco=false;

   fWeight=1.0;

   //default mode for real data
   fMode = "NONE";

   /////////////
   // G4 Info //
   /////////////

   //neutrino that interacted
   fNeutrino.clear();

   //lepton produced in interaction
   fLepton.clear();

   //hyperon produced
   fHyperon.clear();

   //nucleons, pions and kaons produced at primary vtx
   fPrimaryNucleon.clear();
   fPrimaryPion.clear();
   fPrimaryKaon.clear();

   //vertex information
   fTruePrimaryVertex.SetXYZ(-1000,-1000,-1000);

   //hyperon decay products
   fDecayVertex.SetXYZ(-1000,-1000,-1000);
   fDecay.clear();

   //sigma zero decay products
   fSigmaZeroDecayPhoton.clear();
   fSigmaZeroDecayLambda.clear();

   fDecayOpeningAngle=-1; //opening angle between hyperon decay products   
   fLeptonPionAngle=-1; //openining angle between lepton and pion
   fLeptonNucleonAngle=-1; //opening angle between lepton and nucleon

   fKaonDecay.clear();

   ///////////////
   // Reco Info //
   ///////////////

   fSelectedEvent = false; //true if event passes some selection criteria
   fNPrimaryDaughters = 0; //number of primary daughters

   fNPrimaryTrackDaughters=0; //num of track like primary daughters
   fNPrimaryShowerDaughters=0; //num of shower like primary daughters

   fRecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

   fTrackPrimaryDaughters.clear();
   fShowerPrimaryDaughters.clear();
   fMuonIndex=-1;

   //Truth Matched info - indices of muon, proton and pion from decay
   //in fTrackPrimaryDaughters (if they exist)
   fTrueMuonIndex=-1;
   fTrueDecayProtonIndex=-1;
   fTrueDecayPionIndex=-1;

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
