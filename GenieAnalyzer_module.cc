////////////////////////////////////////////////////////////////////////
// Class:       GenieAnalyzer
// Plugin Type: analyzer (art v3_01_02)
// File:        GenieAnalyzer_module.cc
//
// Generated at Wed Sep  4 23:47:09 2019 by Guillermo Fiorentini using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1.h"
#include "TH2.h"
#include "nusimdata/SimulationBase/MCTruth.h"

class GenieAnalyzer;


class GenieAnalyzer : public art::EDAnalyzer {
public:
  explicit GenieAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GenieAnalyzer(GenieAnalyzer const&) = delete;
  GenieAnalyzer(GenieAnalyzer&&) = delete;
  GenieAnalyzer& operator=(GenieAnalyzer const&) = delete;
  GenieAnalyzer& operator=(GenieAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  double fMaxEnergy;
  double fKaonMass;

  TH1F* hEnu;
  TH1F* hMode;
  TH1F* hKaons;
  TH1F* hKaonPlusEnergy;
  TH2F* hEnuMode;
  TH2F* hEnuKaons;
  TH2F* hModeKaons;

};


GenieAnalyzer::GenieAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fMaxEnergy = p.get<float>("MaxNuEnergy", 5.0);
  fKaonMass = 493.677/1000.;//GeV
}

void GenieAnalyzer::analyze(art::Event const& e)
{

  art::Handle< std::vector<simb::MCTruth> > mctruths;
  e.getByLabel("generator", mctruths);

  if (mctruths->size()!=1) {
    std::cout << "Number of MCTruths objects in event: " << mctruths->size() << std::endl;
    return;
  }


  simb::MCTruth mctruth = mctruths->at(0);

  double nu_energy = mctruth.GetNeutrino().Nu().E();
  int nu_pdg = mctruth.GetNeutrino().Nu().PdgCode();
  int ccnc = mctruth.GetNeutrino().CCNC();//0=CC,1=NC
  int mode = mctruth.GetNeutrino().Mode();//0-13

  if( nu_energy >= fMaxEnergy ) nu_energy = fMaxEnergy - 0.001;

  mode = mode+ccnc*15.;//numu CC: 0-14, NC: 15-29
  if ( nu_pdg == -14 ) mode += 1.*30;//anumu CC: 30-44, NC: 45-59
  else if ( nu_pdg == 12 ) mode += 2.*30;//nue CC: 60-74, NC: 75-89
  else if ( nu_pdg == -12 ) mode += 3.*30;//anue CC: 90-104, NC: 105-119

  hEnu->Fill(nu_energy);
  hMode->Fill(mode);
  hEnuMode->Fill(nu_energy,mode);

  //count final state particles
  //and get kaon energy
  int nparticles = mctruth.NParticles();
  int nproton  = 0;
  int nneutron = 0;
  int npip = 0;
  int npim = 0;
  int npi0 = 0;
  int nkp  = 0;
  int nkm  = 0;
  int nk0  = 0;
  int nk0b = 0;
  int nk0s = 0;
  int nk0l = 0;
  std::vector<double> kp_ke;

  for (int i=0; i<nparticles; i++) {

    simb::MCParticle particle = mctruth.GetParticle(i);

    //count number of FSI particles
    if (particle.StatusCode() == 1) {
      if (particle.PdgCode() == 2212) nproton++;
      else if (particle.PdgCode() == 2112) nneutron++;
      else if (particle.PdgCode() == 211) npip++;
      else if (particle.PdgCode() == -211) npim++;
      else if (particle.PdgCode() == 111) npi0++;
      else if (particle.PdgCode() == 130) nk0l++;
      else if (particle.PdgCode() == 310) nk0s++;
      else if (particle.PdgCode() == 311) nk0++;
      else if (particle.PdgCode() == -311) nk0b++;
      else if (particle.PdgCode() == 321) kp_ke.push_back(particle.EndE()-fKaonMass);
      else if (particle.PdgCode() == -321) nkm++;
    }

  }

  nkp = kp_ke.size();

  //kaons
  if (nkp>0){

    //hKaons->Fill(kaon_mode);

    //numu CC with 1 kaon
    hKaonPlusEnergy->Fill(kp_ke[0]);

  }

}

void GenieAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  hEnu = tfs->make<TH1F>("henu", ";Neutrino Energy (GeV);Events",100,0,fMaxEnergy);
  hMode = tfs->make<TH1F>("hmode", ";Interaction Mode;Events",120,0,120);
  hKaons = tfs->make<TH1F>("hkaons", ";Kaon ID;Events",11,0,11);
  hKaonPlusEnergy = tfs->make<TH1F>("hkaonplusenergy", ";Kaon+ Kinetic Energy;Events",100,0,1);
  hEnuMode = tfs->make<TH2F>("henumode", ";Neutrino Energy (GeV);Interaction Mode",100,0,fMaxEnergy,120,0,120);
  hEnuKaons = tfs->make<TH2F>("henukaons", ";Neutrino Energy (GeV);Kaon ID",100,0,fMaxEnergy,11,0,11);
  hModeKaons = tfs->make<TH2F>("hmodekaons", ";Interaction Mode;Kaon ID",120,0,120,11,0,11);

}

void GenieAnalyzer::endJob()
{
  std::cout << "Filled " << hEnu->GetEntries() << " neutrinos" << std::endl;
  std::cout << "CCQES events " << 100.*hMode->GetBinContent(1)/hEnu->GetEntries() << " %" << std::endl;
  std::cout << "CCRES events " << 100.*hMode->GetBinContent(2)/hEnu->GetEntries() << " %" << std::endl;
  std::cout << "CCDIS events " << 100.*hMode->GetBinContent(3)/hEnu->GetEntries() << " %" << std::endl;
  std::cout << "CCCOH events " << 100.*hMode->GetBinContent(4)/hEnu->GetEntries() << " %" << std::endl;
  std::cout << "Events with kaons " << 100.*hEnuKaons->GetEntries()/hEnu->GetEntries() << " %" << std::endl;
  std::cout << "Events with 1 kaon+ " << 100.*(hKaons->GetBinContent(1))/hEnu->GetEntries() << " %" << std::endl;
}

DEFINE_ART_MODULE(GenieAnalyzer)
