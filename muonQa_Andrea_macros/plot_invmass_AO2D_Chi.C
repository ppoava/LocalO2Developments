TFile *fAnalysisResults;

double scaleME = 1;
double massMin = 2.3;
double massMax = 4.1;

// Global canvas and legend definitions to compare different data sets and geometries
TCanvas *c_MuonKine_MuonCuts = new TCanvas("c_MuonKine_MuonCuts", "c_MuonKine_MuonCuts", 1200, 800);
TCanvas c("c", "c", 1200, 800);

TLegend leg_MuonKine_MuonCuts;
std::vector<TPaveText*> all_texts;

TH1* GetTH1(TFile* f, TString histname)
{
  std::cout << "Reading " << histname << " from TFile" << std::endl;
  TH1* hist = (TH1*)f->Get(histname);
  if (hist == nullptr) { std::cout << ">> error retrieving histogram" << std::endl; }
  else { std::cout << ">> histogram sucessfully read from TFile" << std::endl; }
  return hist;

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  //std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  std::cout << "Reading " << histname << " from TFile" << std::endl;
  TH2* hist = (TH2*)f->Get(histname);
  if (hist == nullptr) { std::cout << ">> error retrieving histogram" << std::endl; }
  else { std::cout << ">> histogram sucessfully read from TFile" << std::endl; }
  return hist;

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  //std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}

void NormalizeInvmassME(TH1* h, TH1* hME)
{
  int binMin = h->GetXaxis()->FindBin(2.0);
  int binMax = h->GetXaxis()->FindBin(2.5);
  double integral = h->Integral(binMin, binMax);
  double integralME = hME->Integral(binMin, binMax);
  if (integralME > 0) hME->Scale(integral / integralME);
}


void RemoveEntries(TH1F* h, float min, float max){
  for (int iBin=h->GetXaxis()->FindBin(min); iBin<h->GetXaxis()->FindBin(max); ++iBin) {
    h->SetBinContent(iBin,0);
    h->SetBinError(iBin,0);
  }
}

int color_matching_type[]={
  TColor::GetColor(95,193,199),
  TColor::GetColor(230,229,0),
  TColor::GetColor(29,71,99)
};

void FitJPsi(TH1* mass_all, TString label, Style_t lineStyle, Width_t lineWidth)
{
  auto bg_mass_all = (TH1F*)mass_all->Clone();
  RemoveEntries(bg_mass_all,2.7,3.3);
  RemoveEntries(bg_mass_all,3.4,3.7);

  auto fit_bg_mass_all = new TF1(TString::Format("fit_bg_mass_all_%s", mass_all->GetName()),"[0]*exp([1]*x) + [2]*exp([3]*x)",massMin,massMax);
  fit_bg_mass_all->SetLineWidth(4);
  bg_mass_all->Fit(fit_bg_mass_all,"RNQ");

  auto fit_sig_mass_all = new TF1(TString::Format("fit_cb_mass_global_all_%s", mass_all->GetName()),"[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])",massMin,massMax);
  fit_sig_mass_all->SetParameters(100, 0.6, -2.13903e+06, 0.05, 3.09);
  fit_sig_mass_all->SetNpx(1000);

  auto sig_mass_all = (TH1F*)mass_all->Clone();
  sig_mass_all->Add(fit_bg_mass_all,-1);
  sig_mass_all->Fit(fit_sig_mass_all,"RNQ");

  auto fit_total_mass_all = new TF1(TString::Format("fit_total_mass_global_all_%s", mass_all->GetName()),"[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]*exp([6]*x) + [7]*exp([8]*x)",massMin,massMax);
  fit_total_mass_all->SetLineWidth(4);
  fit_total_mass_all->SetParameters(fit_sig_mass_all->GetParameter(0),
         fit_sig_mass_all->GetParameter(1),
         fit_sig_mass_all->GetParameter(2),
         fit_sig_mass_all->GetParameter(3),
         fit_sig_mass_all->GetParameter(4),
         fit_bg_mass_all->GetParameter(0),
         fit_bg_mass_all->GetParameter(1),
         fit_bg_mass_all->GetParameter(2),
         fit_bg_mass_all->GetParameter(3));

  fit_total_mass_all->SetLineColor(color_matching_type[2]);
  fit_total_mass_all->SetLineStyle(lineStyle);
  fit_total_mass_all->SetLineWidth(lineWidth);
  fit_bg_mass_all->SetLineColor(color_matching_type[0]);
  fit_bg_mass_all->SetLineStyle(lineStyle);
  fit_bg_mass_all->SetLineWidth(lineWidth);

  mass_all->Fit(fit_total_mass_all,"RNQ");
  mass_all->Fit(fit_total_mass_all,"RNQ");
  TFitResultPtr r = mass_all->Fit(fit_total_mass_all,"RNS");

  fit_bg_mass_all->SetParameters(fit_total_mass_all->GetParameter(5),
            fit_total_mass_all->GetParameter(6),
            fit_total_mass_all->GetParameter(7),
            fit_total_mass_all->GetParameter(8));

  fit_total_mass_all->SetNpx(1000);
  fit_bg_mass_all->SetNpx(1000);

  fit_total_mass_all->Draw("SAME");
  fit_bg_mass_all->Draw("SAME"); 

  double peak = mass_all->GetMaximum();
  mass_all->SetMaximum(peak * 1.5);

  // Simple way to move text for multiple entries
  TPaveText* text;
  if (lineStyle == 1) {
    text = new TPaveText(3.4, peak, 4.05, peak * 1.3);
    text->SetName(Form("fitlabel_%s", label.Data()));
    text->AddText(Form("%s:", label.Data()));
  } else if (lineStyle == 3) {
    text = new TPaveText(3.4, peak * 0.3, 4.05, peak * 0.6);
    text->SetName(Form("fitlabel_%s", label.Data()));
    text->AddText(Form("%s:", label.Data()));
  } else {
    std::cout << "ERROR: linestyle doesn't match function call. Aborting" << std::endl; 
    return;
  }

  text->SetBorderSize(0);
  text->SetFillColor(kWhite);
  //text->AddText(std::format("N(J/#psi) = {:0.2f}", fit_sig_mass_all->Integral(massMin, massMax) / mass_all->GetXaxis()->GetBinWidth(1), fit_sig_mass_all->IntegralError(massMin, massMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())).c_str());
  text->AddText(std::format("#mu = {:0.3f} #pm {:0.3f}", fit_sig_mass_all->GetParameter(4), fit_sig_mass_all->GetParError(4)).c_str());
  text->AddText(std::format("#sigma = {:0.3f} #pm {:0.3f}", fit_sig_mass_all->GetParameter(3), fit_sig_mass_all->GetParError(3)).c_str());
  text->Draw("SAME"); 

  all_texts.push_back(text);
}

void skeleton_plot_invmass(std::vector<TString> vfNames, std::vector<TString> vLabels, std::vector<Style_t> vLineStyles, std::vector<Width_t> vLineWidths)
{
  for (Int_t i = 0; i < vfNames.size(); i++) {
    TString fName = vfNames[i];
    TString label = vLabels[i];
    Style_t lineStyle = vLineStyles[i];
    Width_t lineWidth = vLineWidths[i];

    fAnalysisResults = new TFile(fName);

    int rebin = 1;

    gStyle->SetOptStat(0);
    //gStyle->SetOptStat(1111111);
    gStyle->SetOptFit(1111);

    c_MuonKine_MuonCuts->cd();
    TH1* invmassMuonTracks1 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_MuonCuts");
    TH1* invmassMuonTracks1ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_MuonCuts");
    NormalizeInvmassME(invmassMuonTracks1, invmassMuonTracks1ME);
    invmassMuonTracks1->SetLineColor(kRed);
    invmassMuonTracks1->SetLineStyle(lineStyle);
    invmassMuonTracks1->SetLineWidth(lineWidth);
    invmassMuonTracks1->Rebin(rebin);
    invmassMuonTracks1->Scale(1.0 / invmassMuonTracks1->Integral());
    invmassMuonTracks1->GetXaxis()->SetRangeUser(massMin, massMax);
    if (i == 0) { invmassMuonTracks1->Draw("HIST E"); 
    } else { invmassMuonTracks1->Draw("SAME HIST E"); }
    //invmassMuonTracks1ME->Draw("same");
    FitJPsi(invmassMuonTracks1, label, lineStyle, lineWidth);
    leg_MuonKine_MuonCuts.AddEntry(invmassMuonTracks1, Form("%s", label.Data()), "l");
    leg_MuonKine_MuonCuts.Draw();

    c.cd();
    TH1* invmassMuonTracks2 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts");
    TH1* invmassMuonTracks2ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts");
    NormalizeInvmassME(invmassMuonTracks2, invmassMuonTracks2ME);
    invmassMuonTracks2->SetLineColor(kBlue - 2);
    invmassMuonTracks2->SetLineStyle(lineStyle);
    invmassMuonTracks2->SetLineWidth(lineWidth);
    invmassMuonTracks2->Rebin(rebin);
    invmassMuonTracks2->GetXaxis()->SetRangeUser(massMin, massMax);
    invmassMuonTracks2->Draw("E");
    //invmassMuonTracks2ME->Draw("same");
    FitJPsi(invmassMuonTracks2, label, lineStyle, lineWidth);
    c.SaveAs("Plots/invmass_AO2D_MuonKine_GlobalMuonCuts.pdf");

    TH1* invmassMuonTracks3 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts");
    TH1* invmassMuonTracks3ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts");
    NormalizeInvmassME(invmassMuonTracks3, invmassMuonTracks3ME);
    invmassMuonTracks3->SetLineColor(kGreen - 2);
    invmassMuonTracks3->SetLineStyle(lineStyle);
    invmassMuonTracks3->SetLineWidth(lineWidth);
    invmassMuonTracks3->Rebin(rebin);
    invmassMuonTracks3->GetXaxis()->SetRangeUser(massMin, massMax);
    invmassMuonTracks3->Draw("E");
    //invmassMuonTracks3ME->Draw("same");
    FitJPsi(invmassMuonTracks3, label, lineStyle, lineWidth);
    c.SaveAs("Plots/invmass_AO2D_MuonKine_GlobalMatchesCuts.pdf");

    TH1* invmassMuonTracks4 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts");
    TH1* invmassMuonTracks4ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts");
    NormalizeInvmassME(invmassMuonTracks4, invmassMuonTracks4ME);
    invmassMuonTracks4->SetLineColor(kMagenta);
    invmassMuonTracks4->SetLineStyle(lineStyle);
    invmassMuonTracks4->SetLineWidth(lineWidth);
    invmassMuonTracks4->Rebin(rebin);
    invmassMuonTracks4->GetXaxis()->SetRangeUser(massMin, massMax);
    invmassMuonTracks4->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmassMuonTracks4, label, lineStyle, lineWidth);
    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_GlobalMatchesCuts.pdf");

    TH1* invmassMuonTracks5 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts");
    TH1* invmassMuonTracks5ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts");
    NormalizeInvmassME(invmassMuonTracks5, invmassMuonTracks5ME);
    invmassMuonTracks5->SetLineColor(kBlack);
    invmassMuonTracks5->SetLineStyle(lineStyle);
    invmassMuonTracks5->SetLineWidth(lineWidth);
    invmassMuonTracks5->Rebin(rebin);
    invmassMuonTracks5->GetXaxis()->SetRangeUser(massMin, massMax);
    invmassMuonTracks5->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmassMuonTracks5, label, lineStyle, lineWidth);
    c.SaveAs("Plots/invmass_AO2D_ScaledMftKine_GlobalMatchesCuts.pdf");

    // TODO: add top-top/top-bottom/bottom-bottom and left-right/left-left/right-right combinations here!!
    invmassMuonTracks5->Draw("hist");
    invmassMuonTracks4->Draw("same");
    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_vs_ScaledMftKine_and_GlobalMatchesCuts.pdf");

    c.Clear();
    c.Divide(2,2);
    // top-top
    c.cd(1);
    TH1* invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT");
    TH1* invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT");
    if (invmass) {
    NormalizeInvmassME(invmass, invmassME);
    invmass->SetLineColor(kBlack);
    invmass->SetLineStyle(lineStyle);
    invmass->SetLineWidth(lineWidth);
    invmass->Rebin(rebin);
    invmass->GetXaxis()->SetRangeUser(massMin, massMax);
    invmass->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmass, label, lineStyle, lineWidth);
    }
    // top-bottom
    c.cd(2);
    invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB");
    invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB");
    if (invmass) {
    NormalizeInvmassME(invmass, invmassME);
    invmass->SetLineColor(kBlack);
    invmass->SetLineStyle(lineStyle);
    invmass->SetLineWidth(lineWidth);
    invmass->Rebin(rebin);
    invmass->GetXaxis()->SetRangeUser(massMin, massMax);
    invmass->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmass, label, lineStyle, lineWidth);
    }
    // bottom-top
    c.cd(3);
    invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BT");
    invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BT");
    if (invmass) {
    NormalizeInvmassME(invmass, invmassME);
    invmass->SetLineColor(kBlack);
    invmass->SetLineStyle(lineStyle);
    invmass->SetLineWidth(lineWidth);
    invmass->Rebin(rebin);
    invmass->GetXaxis()->SetRangeUser(massMin, massMax);
    invmass->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmass, label, lineStyle, lineWidth);
    }
    // bottom-bottom
    c.cd(4);
    invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB");
    invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB");
    if (invmass) {
    NormalizeInvmassME(invmass, invmassME);
    invmass->SetLineColor(kBlack);
    invmass->SetLineStyle(lineStyle);
    invmass->SetLineWidth(lineWidth);
    invmass->Rebin(rebin);
    invmass->GetXaxis()->SetRangeUser(massMin, massMax);
    invmass->Draw("E");
    //invmassMuonTracks4ME->Draw("same");
    FitJPsi(invmass, label, lineStyle, lineWidth);
    }

    c.SaveAs("Plots/invmass_AO2D_ScaledMftKine_GlobalMatchesCuts_TB.pdf");

    c.Clear();

    TH1* histogram;
    TH2* histogram2;

    histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_vs_GlobalMuonKine");
    histogram2->GetXaxis()->SetRangeUser(massMin, massMax);
    histogram2->GetYaxis()->SetRangeUser(massMin, massMax);
    histogram2->SetLineStyle(lineStyle);
    histogram2->SetLineWidth(lineWidth);
    histogram2->Draw("col");

    c.SaveAs("Plots/invmass_AO2D_MuonKine_vs_GlobalMuonKine.pdf");

    histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_vs_GlobalMuonKine");
    histogram2->GetXaxis()->SetRangeUser(massMin, massMax);
    histogram2->GetYaxis()->SetRangeUser(massMin, massMax);
    histogram2->SetLineStyle(lineStyle);
    histogram2->SetLineWidth(lineWidth);
    histogram2->Draw("col");

    c.SaveAs("Plots/invmass_AO2D_ScaledMftKine_vs_GlobalMuonKine.pdf");

    c.SetLogy(kTRUE);
    histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMassFull_MuonKine_MuonCuts");
    histogram->Draw();
    histogram->SetLineStyle(lineStyle);
    histogram->SetLineWidth(lineWidth);
    c.SaveAs("Plots/invmass_AO2D_MuonKine_MuonCuts_log.pdf");
    histogram->GetXaxis()->SetRangeUser(1.0, 15.0);
    c.SaveAs("Plots/invmass_AO2D_MuonKine_MuonCuts_log_zoomed.pdf");

    c.SetLogy(kTRUE);
    histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts");
    histogram->Draw();
    histogram->SetLineStyle(lineStyle);
    histogram->SetLineWidth(lineWidth);
    c.SaveAs("Plots/invmass_AO2D_ScaledMftKine_GlobalMatchesCuts_log.pdf");
    histogram->GetXaxis()->SetRangeUser(1.0, 15.0);
    c.SaveAs("Plots/invmass_AO2D_ScaledMftKine_GlobalMatchesCuts_log_zoomed.pdf");

    TH1* histogramLSL = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_leading_subleading");
    histogramLSL->Rebin(4);
    histogramLSL->SetLineStyle(lineStyle);
    histogramLSL->SetLineWidth(lineWidth);
    histogramLSL->Draw();

    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_GlobalMatchesCuts_leading_subleading.pdf");

    TH1* histogramSLL = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_leading");
    histogramSLL->Rebin(4);
    histogramSLL->SetLineStyle(lineStyle);
    histogramSLL->SetLineWidth(lineWidth);
    histogramSLL->Draw();

    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_GlobalMatchesCuts_subleading_leading.pdf");

    histogramLSL->Add(histogramSLL);
    histogramLSL->SetTitle("M_{#mu^{+}#mu^{-}} - leading + sub-leading matches");
    histogramLSL->SetLineStyle(lineStyle);
    histogramLSL->SetLineWidth(lineWidth);
    histogramLSL->Draw();

    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_GlobalMatchesCuts_subleading_leading_no_rebin.pdf");

    histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading");
    histogram->SetTitle("M_{#mu^{+}#mu^{-}} - sub-leading matches");
    histogram->Rebin(4);
    histogram->SetLineStyle(lineStyle);
    histogram->SetLineWidth(lineWidth);
    histogram->Draw();

    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading.pdf");

    histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_subleading_vs_leading");
    histogram2->GetXaxis()->SetRangeUser(1, 5);
    histogram2->GetYaxis()->SetRangeUser(1, 5);
    histogram2->SetLineStyle(lineStyle);
    histogram2->SetLineWidth(lineWidth);
    histogram2->Draw("col");

    c.SaveAs("Plots/invmass_AO2D_GlobalMuonKine_subleading_vs_leading.pdf)");
  }
}

void plot_invmass_AO2D_Chi() 
{
  // TODO: make it read and save on the same canvas
  // Maybe easiest to define canvases above function (with different names) and read them in the function?
  std::vector<TString> vfNames;
  std::vector<TString> vLabels;
  std::vector<Style_t> vLineStyles;
  std::vector<Width_t> vLineWidths;

  vfNames.push_back("AnalysisResults_LHC24aq_pass1_small_muon_qa_test.root");
  vLabels.push_back("LHC24aq_pass1_small");
  vLineStyles.push_back(1);
  vLineWidths.push_back(1);

  vfNames.push_back("AnalysisResults_LHC24an_pass1_skimmed_small_muon_qa_test.root");
  vLabels.push_back("LHC24an_pass1_skimmed_small");
  vLineStyles.push_back(3);
  vLineWidths.push_back(2);

  // skeleton_plot_invmass("AnalysisResults_LHC24aq_pass1_small_muon_qa_test.root", "LHC24aq_pass1_small", 1, 1);
  // skeleton_plot_invmass("AnalysisResults_LHC24an_pass1_skimmed_small_muon_qa_test.root", "LHC24an_pass1_skimmed_small", 3, 2);
  skeleton_plot_invmass(vfNames, vLabels, vLineStyles, vLineWidths);
  c_MuonKine_MuonCuts->Draw();
  c_MuonKine_MuonCuts->SaveAs("Plots/invmass_AO2D_MuonKine_MuonCuts.pdf");

  for (size_t i = 0; i < all_texts.size(); ++i) 
  {
    TPaveText* t = all_texts[i];
    if (!t) {
      std::cout << "Text " << i << " is a null pointer!\n";
      continue;
    }
    std::cout << "Text " << i << ": name = " << t->GetName()
              << ", entries = " << t->GetListOfLines()->GetSize()
              << ", x1 = " << t->GetX1NDC()
              << ", x2 = " << t->GetX2NDC()
              << ", y1 = " << t->GetY1NDC()
              << ", y2 = " << t->GetY2NDC()
              << std::endl;
}

}