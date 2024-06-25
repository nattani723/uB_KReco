#include "AngularDistributionDrawer.h"

namespace kaon_reconstruction
{



  AngularDistributionDrawer::AngularDistributionDrawer(const ParticleDirectionFinder& particle_direction_finder) :
  m_theta_bin_size(particle_direction_finder.get_theta_bin_size()),
    m_phi_bin_size(particle_direction_finder.get_phi_bin_size()),
    m_num_bin_theta(static_cast<int>(M_PI / particle_direction_finder.get_theta_bin_size())),
    m_num_bin_phi(static_cast<int>(2*M_PI / particle_direction_finder.get_phi_bin_size()))
      {
      }


  //------------------------------------------------------------------------------------------------------------------------------------------   

  void AngularDistributionDrawer::runDrawer(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map, TCanvas* &c)
  {

    this->fill_angular_distribution_map_cheated_pdg(sp_list_roi, k_end, angular_distribution_map_cheated_pdg, h_angular_distribution_cheated_pdg, spacepointToHitMap, hit_pdg_map);

    this->draw_hist_angular_distribution_map_cheated_pdg(h_angular_distribution_cheated_pdg, "test.root", c);

  } 

  //------------------------------------------------------------------------------------------------------------------------------------------      

  void AngularDistributionDrawer::fill_angular_distribution_map_cheated_pdg(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map) const
  {
    for(auto it_sp = sp_list_roi.begin(); it_sp != sp_list_roi.end(); ++it_sp){
      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      double theta = displacement_vector.Theta();
      double phi = displacement_vector.Phi();
      int theta_factor = (int)(std::floor(theta / m_theta_bin_size));
      int phi_factor = (int)(std::floor(phi / m_phi_bin_size));

      // retrieve PDG info
      art::Ptr<recob::Hit> phit = spacepointToHitMap.at(*it_sp);
      int pdg = hit_pdg_map[*phit];

      angular_distribution_map_cheated_pdg[pdg][theta_factor][phi_factor] += TMath::Sin(theta);
    }

    if(angular_distribution_map_cheated_pdg.empty()) return;

    for(const auto& angular_distribution_map : angular_distribution_map_cheated_pdg) {
      TH2D* h = new TH2D("", "; #Theta [rad]; #Phi [rad]; Hits [a.u.]", m_num_bin_theta, 0, M_PI, m_num_bin_phi, -M_PI, M_PI);
      for(const auto& theta_entry : angular_distribution_map.second) {
	for(const auto& phi_entry : theta_entry.second) {
	  h->Fill(m_theta_bin_size * theta_entry.first, m_phi_bin_size * phi_entry.first, phi_entry.second);
	}
      }
      h_angular_distribution_cheated_pdg[angular_distribution_map.first] = h;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void AngularDistributionDrawer::draw_hist_angular_distribution_map_cheated_pdg(std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, TString outfile_name, TCanvas* &c) const
  {
    if(h_angular_distribution_cheated_pdg.empty()) return;

    TFile outfile(outfile_name, "update");
    auto h_stack = new THStack("h_stack", "");

    c->Clear();
    c->SetFillStyle(1001);

    TLegend * leg = new TLegend(0.1, 0.0, 0.9, 1.0, "");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.2);  
    leg->SetNColumns(2);

    const double Single_PadSplit = 0.85;
    TPad *p_plot = new TPad("p_plot","p_plot",0,0,1,Single_PadSplit);
    TPad *p_legend = new TPad("p_legend","p_legend",0,Single_PadSplit,1,1);
    p_legend->SetBottomMargin(0);
    p_legend->SetTopMargin(0.1);
    p_plot->SetTopMargin(0.01);

    std::set<int> added_pids;

    for(const auto& it_pdg : h_angular_distribution_cheated_pdg) {
      int pdg_code = it_pdg.first;
      TH2D* histogram = it_pdg.second;
      histogram->SetFillStyle(1001);

      int fill_color = kGray; // Default color

      switch (pdg_code) {
      case 321:
	fill_color = kBlue;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "K^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case -13:
	fill_color = kCyan;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "#mu^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case 211:
	fill_color = kMagenta;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "#pi^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case -11:
      case 11:
	fill_color = kGreen + 2;
	if (added_pids.find(11) == added_pids.end() && added_pids.find(-11) == added_pids.end()) {
	  leg->AddEntry(histogram, "Shower", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case 2212:
      case 2112:
	fill_color = kRed;
	if (added_pids.find(2212) == added_pids.end() && added_pids.find(2112) == added_pids.end() ) {
	  leg->AddEntry(histogram, "Nucleon", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      default:
	if (added_pids.find(999) == added_pids.end()) {
	  leg->AddEntry(histogram, "Others", "f");
	  added_pids.insert(999);
	}
	break; // Use default color
      }

      histogram->SetFillColor(fill_color);
      histogram->SetFillStyle(1001);
      h_stack->Add(histogram);
    }

    p_legend->Draw();
    p_legend->cd();
    leg->Draw();
    c->cd();
    p_plot->Draw();
    p_plot->cd();
    h_stack->Draw("hist lego3 0");

    h_stack->GetXaxis()->SetTitle("#Phi [rad]");
    h_stack->GetYaxis()->SetTitle("#Theta [rad]");
    h_stack->GetXaxis()->SetTitleOffset(1.5);
    h_stack->GetYaxis()->SetTitleOffset(1.5);

    //h_stack->GetZaxis()->SetTitle("Event [a.u.]");
    c->Modified();
    c->Write();

  }

}//namespace kaon_reconstruction 
