TFile* fAnalysisResults;

double scaleME = 1;

TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}


void plotMFTDCAProjection(TH2* histogram2, float xMin, float xMax, float yMin, float yMax, TCanvas& c, bool printFits = false)
{
  histogram2->GetXaxis()->SetRange(0, 0);

  TH1* histogram = histogram2->ProfileX();
  histogram->SetTitle(TString::Format("%s (profile)", histogram2->GetTitle()));
  histogram->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogram->GetXaxis()->SetRangeUser(xMin, xMax);
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetStats(kTRUE);
  histogram->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", histogram2->GetName()));
  histogramMean->SetTitle(TString::Format("%s (mean)", histogram2->GetTitle()));
  histogramMean->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", histogram2->GetName()));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = (histogram2->GetEntries() / histogram2->GetXaxis()->GetNbins()) / 2;
  std::cout << std::format("[\"{}\"] Nbins {} {}\n", histogram2->GetName(), histogram2->GetXaxis()->GetNbins(), histogramMean->GetXaxis()->GetNbins());
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", histogram2->GetName(), bin), bin, bin);
    std::cout << std::format("[\"{}\"] bin {} entries {} min {}\n", histogram2->GetName(), bin, proj->GetEntries(), entriesMin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > 1) {
      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      TF1 fgaus2("fgaus2", "gaus(0)", xPeak - 0.1, xPeak + 0.1);
      //TF1 fgaus("fgaus", "gaus(0)+gaus(3)");
      fgaus2.SetParameter(1, xPeak);
      fgaus2.SetParameter(2, 0.01);
      fgaus2.SetParLimits(2, 0, 0.1);
      //fgaus.SetParameter(5, 0.1);
      //fgaus.SetParLimits(5, 0.1, 100);
      proj->Fit("fgaus2", "NBRQ");

      xPeak = fgaus2.GetParameter(1);
      TF1 fgaus("fgaus", "gaus(0)+gaus(3)", xPeak - 1, xPeak + 1);
      //TF1 fgaus("fgaus", "gaus(0)+gaus(3)");
      fgaus.SetParameter(0, fgaus2.GetParameter(0));
      fgaus.SetParameter(1, xPeak);
      fgaus.SetParameter(2, fgaus2.GetParameter(2));
      fgaus.SetParLimits(2, fgaus2.GetParameter(2) / 2, fgaus2.GetParameter(2) * 2);
      fgaus.SetParameter(3, fgaus2.GetParameter(0) / 10);
      fgaus.SetParLimits(3, 0, fgaus2.GetParameter(0) * 100);
      fgaus.SetParameter(4, xPeak);
      fgaus.SetParLimits(4, xPeak - 0.1, xPeak + 0.1);
      fgaus.SetParameter(5, fgaus2.GetParameter(2) * 10);
      fgaus.SetParLimits(5, 0.1, 100);
      proj->Fit("fgaus", "BR");

      if (printFits) {
        proj->GetXaxis()->SetRangeUser(xPeak - 2, xPeak + 2);
        proj->Draw("E");
        c.SaveAs("alignment_AO2D.pdf");
      }

      mean = fgaus.GetParameter(1);
      meanErr = fgaus.GetParError(1);
      sigma = fgaus.GetParameter(2);
      sigmaErr = fgaus.GetParError(2);
    }

    //std::cout << std::format("histogramMean->SetBinContent({}, {:0.2f})\n", bin, mean);
    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->GetXaxis()->SetRangeUser(xMin, xMax);
  histogramMean->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
  histogramSigma->GetXaxis()->SetRangeUser(xMin, xMax);
  histogramSigma->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
}


void plotProjection(TH2* histogram2, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  TH1* histogram = histogram2->ProfileX();
  histogram->SetTitle(TString::Format("%s (profile)", histogram2->GetTitle()));
  histogram->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetStats(kTRUE);
  histogram->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", histogram2->GetName()));
  histogramMean->SetTitle(TString::Format("%s (mean)", histogram2->GetTitle()));
  histogramMean->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", histogram2->GetName()));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = (histogram2->GetEntries() / histogram2->GetXaxis()->GetNbins()) / 2;
  std::cout << std::format("[\"{}\"] Nbins {} {}\n", histogram2->GetName(), histogram2->GetXaxis()->GetNbins(), histogramMean->GetXaxis()->GetNbins());
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", histogram2->GetName(), bin), bin, bin);
    std::cout << std::format("[\"{}\"] bin {} entries {} min {}\n", histogram2->GetName(), bin, proj->GetEntries(), entriesMin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    proj->GetXaxis()->SetRangeUser(yMin, yMax);
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > 1) {
      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      TF1 fgaus("fgaus", "gaus");
      fgaus.SetParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "BQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("alignment_AO2D.pdf");
      }

      mean = fgaus.GetParameter(1);
      meanErr = fgaus.GetParError(1);
      sigma = fgaus.GetParameter(2);
      sigmaErr = fgaus.GetParError(2);
    }

    //std::cout << std::format("histogramMean->SetBinContent({}, {:0.2f})\n", bin, mean);
    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
  histogramSigma->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
}


void plotProjection(TH2* histogram2, TH2* histogram2ME, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  TH1* histogram = histogram2->ProfileX();
  histogram->SetTitle(TString::Format("%s (profile)", histogram2->GetTitle()));
  histogram->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetStats(kTRUE);
  histogram->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", histogram2->GetName()));
  histogramMean->SetTitle(TString::Format("%s (mean)", histogram2->GetTitle()));
  histogramMean->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", histogram2->GetName()));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->GetYaxis()->SetTitle(histogram2->GetYaxis()->GetTitle());
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = (histogram2->GetEntries() / histogram2->GetXaxis()->GetNbins()) / 2;
  std::cout << std::format("[\"{}\"] Nbins {} {}\n", histogram2->GetName(), histogram2->GetXaxis()->GetNbins(), histogramMean->GetXaxis()->GetNbins());
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", histogram2->GetName(), bin), bin, bin);
    TH1* projME = histogram2ME->ProjectionY(TString::Format("%s-%d-ME", histogram2->GetName(), bin), bin, bin);
    std::cout << std::format("[\"{}\"] bin {} entries {} min {}\n", histogram2->GetName(), bin, proj->GetEntries(), entriesMin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > 1) {
      //proj->Rebin(projRebin);

      if (printFits) {
        proj->SetLineColor(kRed);
        proj->Draw("E");
        projME->Draw("E same");
        c.SaveAs("alignment_AO2D.pdf");
      }

      proj->Add(projME, -1);

      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      TF1 fgaus("fgaus", "gaus(0)+pol1(3)");
      fgaus.SetParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "BQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("alignment_AO2D.pdf");
      }

      mean = fgaus.GetParameter(1);
      meanErr = fgaus.GetParError(1);
      sigma = fgaus.GetParameter(2);
      sigmaErr = fgaus.GetParError(2);
    }

    //std::cout << std::format("histogramMean->SetBinContent({}, {:0.2f})\n", bin, mean);
    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
  histogramSigma->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
}

void plotXYDistribution(const char* histName, float yMin, float yMax, int xRebin, int yRebin, TCanvas& c, bool printFits = false)
{
  std::string histNameFull = std::string("muon-extrap/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, histNameFull);
  std::cout << std::format("[TOTO] histNameFull={}  hist={}\n", histNameFull, (void*)histogram2);
  if (!histogram2) return;
  // Mixed events histogram
  std::string histNameFullME = std::string("muon-extrap/ME/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, histNameFullME);
  std::cout << std::format("[TOTO] histNameFullME={}  histME={}\n", histNameFullME, (void*)histogram2ME);
  if (!histogram2ME) return;
  histogram2->RebinX(xRebin);
  histogram2->RebinY(yRebin);
  histogram2ME->RebinX(xRebin);
  histogram2ME->RebinY(yRebin);
  histogram2ME->Scale(scaleME);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  histogram2ME->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  plotProjection(histogram2, histogram2ME, yMin, yMax, yRebin, c, printFits);
}

void plotPhiDistribution(const char* histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::string histNameFull = std::string("muon-extrap/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, histNameFull);
  std::cout << std::format("[TOTO] histNameFull={}  hist={}\n", histNameFull, (void*)histogram2);
  if (!histogram2) return;

  // keep only two bins for each phi quadrant
  int xRebin = histogram2->GetXaxis()->GetNbins() / 4;
  histogram2->RebinX(xRebin);
  histogram2->RebinY(projRebin);
  histogram2->GetXaxis()->SetBinLabel(1, "bottom-left");
  histogram2->GetXaxis()->SetBinLabel(2, "bottom-right");
  histogram2->GetXaxis()->SetBinLabel(3, "top-right");
  histogram2->GetXaxis()->SetBinLabel(4, "top-left");
  histogram2->GetXaxis()->SetTitle("");
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");

  // Mixed events histogram
  std::string histNameFullME = std::string("muon-extrap/ME/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, histNameFullME);
  std::cout << std::format("[TOTO] histNameFullME={}  histME={}\n", histNameFullME, (void*)histogram2ME);
  if (histogram2ME) {
    histogram2ME->RebinX(xRebin);
    histogram2ME->RebinY(projRebin);
    histogram2ME->Scale(scaleME);
    histogram2ME->Draw("col");
    c.SaveAs("alignment_AO2D.pdf");
    plotProjection(histogram2, histogram2ME, yMin, yMax, projRebin, c, printFits);
  } else {
    plotProjection(histogram2, yMin, yMax, projRebin, c, printFits);
  }
}

void plotPhiDistributionWithSign(const char* histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  plotPhiDistribution(histName, yMin, yMax, projRebin, c, printFits);
  plotPhiDistribution((std::string(histName) + "Plus").c_str(), yMin, yMax, projRebin, c, printFits);
  plotPhiDistribution((std::string(histName) + "Minus").c_str(), yMin, yMax, projRebin, c, printFits);
}

void plotPhiDistributionForResiduals(const char* histName, int deId, TCanvas& c, bool printFits = false)
{
  TH2* histogram2 = GetTH2(fAnalysisResults, histName);
  histogram2->RebinX(10);
  histogram2->RebinY(4);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  for (int xbin = 1; xbin <= histogram2->GetXaxis()->GetNbins(); xbin++) {
    double binCenter = histogram2->GetXaxis()->GetBinCenter(xbin);
    int quadrant = deId % 100;
    double phiMin;
    double phiMax;
    if (quadrant == 0) {
      // top-right, keep bins corresponding to phi from 0 to 90
      phiMin = 0;
      phiMax = 90;
    }
    if (quadrant == 1) {
      // top-left, keep bins corresponding to phi from 90 to 180
      phiMin = 90;
      phiMax = 180;
    }
    if (quadrant == 2) {
      // bottom-left, keep bins corresponding to phi from -180 to -90
      phiMin = -180;
      phiMax = -90;
    }
    if (quadrant == 3) {
      // bottom-right, keep bins corresponding to phi from -90 to 0
      phiMin = -90;
      phiMax = 0;
    }

    phiMin += 10;
    phiMax -= 10;

    if (binCenter > phiMin && binCenter < phiMax) {
      continue;
    }

    for (int ybin = 1; ybin <= histogram2->GetYaxis()->GetNbins(); ybin++) {
      histogram2->SetBinContent(xbin, ybin, 0);
      histogram2->SetBinError(xbin, ybin, 0);
    }
  }
  plotProjection(histogram2, -5.0, 5.0, 4, c, printFits);
}

TF1 plotProjectionForResiduals(TH2* histogram2, TH2* histogram2EM, float avgMult, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  TH1* histogram = histogram2->ProfileX();
  histogram->SetTitle(TString::Format("%s (profile)", histogram2->GetTitle()));
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetStats(kTRUE);
  histogram->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", histogram2->GetName()));
  histogramMean->SetTitle(TString::Format("%s (mean)", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", histogram2->GetName()));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  TH1* histogramFull = histogram2->ProjectionY();
  TH1* histogramFullEM = histogram2EM->ProjectionY();

  int valuePeak = histogramFull->GetMaximum();
  int binPeak = histogramFull->GetMaximumBin();
  double xPeak = histogramFull->GetXaxis()->GetBinCenter(binPeak);

  TF1 fgausFull("fgausFull", "gaus(0) + gaus(3) + pol0(6)");
  //fgaus.SetParLimits(0, 0, 100000000000.0);
  fgausFull.SetParameter(1, xPeak);
  fgausFull.SetParameter(2, 1);
  fgausFull.SetParLimits(2, 0, 100);
  fgausFull.SetParameter(4, xPeak);
  fgausFull.SetParameter(5, 10);
  fgausFull.SetParLimits(5, 2, 100);
  histogramFull->Fit("fgausFull", "QB");

  histogramFull->Draw("E");
  histogramFullEM->Scale(1.0 / avgMult);
  histogramFullEM->SetLineColor(kRed);
  //histogramFullEM->Draw("H same");
  c.SaveAs("alignment_AO2D.pdf");


  int entriesMin = histogram2->GetEntries() / 20;
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", histogram2->GetName(), bin), bin, bin);
    TH1* projEM = histogram2EM->ProjectionY(TString::Format("%s-%d", histogram2EM->GetName(), bin), bin, bin);
    projEM->Scale(1.0 / avgMult);
    std::cout << std::format("[\"{}\"] bin {} entries {} min {}\n", histogram2->GetName(), bin, proj->GetEntries(), entriesMin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > entriesMin) {

      if (printFits) {
        proj->Draw("E");
        projEM->SetLineColor(kRed);
        projEM->Draw("L same");
        c.SaveAs("alignment_AO2D.pdf");
      }

      proj->Add(projEM, -1);

      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);
      TF1 fgaus("fgaus", "gaus", xPeak - 1, xPeak + 1);
      //fgaus.SetParLimits(0, 0, 100000000000.0);
      fgaus.FixParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "RQBN");
      fgaus.ReleaseParameter(1);
      proj->Fit("fgaus", "RQBN");

      TF1 fgaus2("fgaus2", "gaus(0)+pol0(3)");
      fgaus2.SetParameter(0, fgaus.GetParameter(0));
      fgaus2.SetParameter(1, xPeak);
      fgaus2.SetParameter(2, fgaus.GetParameter(2));
      proj->Fit("fgaus2", "BQN");
      proj->Fit("fgaus2", "BQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("alignment_AO2D.pdf");
      }

      mean = fgaus2.GetParameter(1);
      meanErr = fgaus2.GetParError(1);
      sigma = fgaus2.GetParameter(2);
      sigmaErr = fgaus2.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
  histogramSigma->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");

  return fgausFull;
}

TF1 plotXYDistributionForResiduals(std::string histName, int deId, TCanvas& c, bool printFits = false)
{
  TH2* histogram2 = GetTH2(fAnalysisResults, histName.c_str());
  if (!histogram2) return {};
  histogram2->RebinX(10);
  histogram2->RebinY(4);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");

  TH2* histogram2EM = GetTH2(fAnalysisResults, (histName + "EM").c_str());
  histogram2EM->RebinX(10);
  histogram2EM->RebinY(4);
  histogram2EM->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");

  TH1* histogramMult = GetTH1(fAnalysisResults, std::format("muon-extrap/multForDE{}EM", deId).c_str());
  histogramMult->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  return plotProjectionForResiduals(histogram2, histogram2EM, histogramMult->GetMean(), -5.0, 5.0, 4, c, printFits);
}


void plot_alignment_AO2D()
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24l7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC22p/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23zk/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-4/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-7/AnalysisResultsFull.root");

  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  
  TH1* histogram;
  TH2* histogram2;

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtMFTVsPhiAtAbs");
  double integral1 = histogram2->Integral(1, histogram2->GetXaxis()->GetNbins(), 1, histogram2->GetYaxis()->FindBin(-5));
  double integral2 = histogram2->Integral(1, histogram2->GetXaxis()->GetNbins(), histogram2->GetYaxis()->FindBin(5), histogram2->GetYaxis()->GetNbins());
  TH2* histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtMFTVsPhiAtAbs");
  double integralME1 = histogram2ME->Integral(1, histogram2ME->GetXaxis()->GetNbins(), 1, histogram2ME->GetYaxis()->FindBin(-5));
  double integralME2 = histogram2ME->Integral(1, histogram2ME->GetXaxis()->GetNbins(), histogram2ME->GetYaxis()->FindBin(5), histogram2ME->GetYaxis()->GetNbins());

  scaleME = (integral1 + integral2) / (integralME1 + integralME2);

  c.SaveAs("alignment_AO2D.pdf(");

  c.Divide(2, 2);

  c.cd(1);
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/xVsyMFTAtMFT");
  histogram2->GetXaxis()->SetRangeUser(-20, 20);
  histogram2->GetYaxis()->SetRangeUser(-20, 20);
  histogram2->Draw("col");
  c.cd(3);
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/xVsyMCHAtMFT");
  histogram2->GetXaxis()->SetRangeUser(-20, 20);
  histogram2->GetYaxis()->SetRangeUser(-20, 20);
  histogram2->Draw("col");
  c.cd(2);
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/xVsyMFTAtAbs");
  histogram2->Draw("col");
  c.cd(4);
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/xVsyMCHAtAbs");
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");

  {
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtMFTVsPhiAtAbs");
  histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtMFTVsPhiAtAbs");
  c.cd(1);
  TH1* proj = histogram2->ProjectionY("dxAtMFTVsPhiAtAbs_py");
  TH1* projME = histogram2ME->ProjectionY("dxAtMFTVsPhiAtAbsME_py");
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");
  c.cd(3);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kBlack);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 100);
  proj2->Fit("fgaus", "BQ");
  proj2->Draw("E");
  }

  {
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtAbsVsPhiAtAbs");
  histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtAbsVsPhiAtAbs");
  c.cd(2);
  TH1* proj = histogram2->ProjectionY("dxAtAbsVsPhiAtAbs_py");
  TH1* projME = histogram2ME->ProjectionY("dxAtAbsVsPhiAtAbsME_py");
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");
  c.cd(4);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kBlack);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 100);
  proj2->Fit("fgaus", "BQ");
  proj2->Draw("E");
  }
  c.SaveAs("alignment_AO2D.pdf");

  {
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dyAtMFTVsPhiAtAbs");
  histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dyAtMFTVsPhiAtAbs");
  c.cd(1);
  TH1* proj = histogram2->ProjectionY("dyAtMFTVsPhiAtAbs_py");
  TH1* projME = histogram2ME->ProjectionY("dyAtMFTVsPhiAtAbsME_py");
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");
  c.cd(3);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kBlack);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 100);
  proj2->Fit("fgaus", "BQ");
  proj2->Draw("E");
  }

  {
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dyAtAbsVsPhiAtAbs");
  histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dyAtAbsVsPhiAtAbs");
  c.cd(2);
  TH1* proj = histogram2->ProjectionY("dyAtAbsVsPhiAtAbs_py");
  TH1* projME = histogram2ME->ProjectionY("dyAtAbsVsPhiAtAbsME_py");
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");
  c.cd(4);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kBlack);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 100);
  proj2->Fit("fgaus", "BQ");
  proj2->Draw("E");
  }
  c.SaveAs("alignment_AO2D.pdf");

  /*
  {
    c.Clear();
    c.Divide(2, 2);

    c.cd(1);
    histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtAbsVsDxAtMFT");
    histogram2->Draw("col");

    c.cd(2);
    histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtAbsVsDxAtMFT");
    histogram2ME->Scale(scaleME);
    histogram2ME->Draw("col");

    c.cd(3);
    TH1* proj = histogram2->ProjectionX("dxAtAbsVsDxAtMFT_px");
    TH1* projME = histogram2ME->ProjectionX("dxAtAbsVsDxAtMFT_pxME");
    //projME->Scale(scaleME);
    proj->SetLineColor(kRed);
    proj->Draw("E");
    projME->Draw("E same");

    c.cd(4);
    proj = histogram2->ProjectionY("dxAtAbsVsDxAtMFT_py");
    projME = histogram2ME->ProjectionY("dxAtAbsVsDxAtMFT_pyME");
    //projME->Scale(scaleME);
    proj->SetLineColor(kRed);
    proj->Draw("E");
    projME->Draw("E same");

    c.SaveAs("alignment_AO2D.pdf");

    c.Clear();
    c.Divide(2, 2);

    c.cd(1);
    histogram2->Add(histogram2ME, -1);
    histogram2->Draw("colz");

    c.cd(2);
    gPad->SetLogz(kTRUE);
    TH2* histogram2b = (TH2*)histogram2->Clone();
    histogram2b->RebinX(5);
    histogram2b->RebinY(5);
    histogram2b->Draw("colz");

    {
      c.cd(3);
      proj = histogram2->ProjectionX("dxAtAbsVsDxAtMFT_px");
      proj->SetLineColor(kBlack);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
      fgaus.SetParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "BQ");
      proj->Draw("E");
    }

    {
      c.cd(4);
      proj = histogram2->ProjectionY("dxAtAbsVsDxAtMFT_py");
      proj->SetLineColor(kBlack);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      TF1 fgaus("fgaus", "gausn(0)+pol1(3)");
      fgaus.SetParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "BQ");
      proj->Draw("E");
    }

    c.SaveAs("alignment_AO2D.pdf");
  }
  */
  c.Clear();

  /*
  //histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtMFTVsyAtAbs");
  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dxAtAbsVsyAtAbs");
  histogram2->RebinX(20);
  histogram2->RebinY(4);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  TH1* hproj = histogram2->ProjectionY("muon-extrap/dxAtMFTVsPhiAtAbs_py");
  hproj->Draw();
  std::cout << std::format("hproj->GetEntries() = {}\n", hproj->GetEntries());
  c.SaveAs("alignment_AO2D.pdf");

  //TH2* histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtMFTVsyAtAbs");
  TH2* histogram2ME = GetTH2(fAnalysisResults, "muon-extrap/ME/dxAtAbsVsyAtAbs");
  histogram2ME->RebinX(20);
  histogram2ME->RebinY(4);
  histogram2ME->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  TH1* hprojME = histogram2ME->ProjectionY("muon-extrap/ME/dxAtMFTVsPhiAtAbs_py");
  hprojME->Draw();
  std::cout << std::format("hproj->GetEntries() = {}\n", hproj->GetEntries());
  c.SaveAs("alignment_AO2D.pdf");

  std::cout << std::format("hproj->GetEntries() = {}\n", hproj->GetEntries());
  double integral = hproj->Integral(1, hproj->GetXaxis()->FindBin(-4));
  //hproj->Scale(1.0 / integral);
  hproj->SetLineColor(kRed);
  hproj->Draw();
  double integralME = hprojME->Integral(1, hprojME->GetXaxis()->FindBin(-4));
  hprojME->Scale(integral / integralME);
  hprojME->Draw("same");
  c.SaveAs("alignment_AO2D.pdf");

  hproj->Add(hprojME, -1);
  TF1 fitFunc("fitFunc", "gaus(0)+pol0(3)");
  fitFunc.SetParameter(1, 0);
  fitFunc.SetParameter(2, 1);
  hproj->Fit("fitFunc", "B");
  hproj->Draw();
  c.SaveAs("alignment_AO2D.pdf");

  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    hproj = histogram2->ProjectionY("muon-extrap/dxAtMFTVsPhiAtAbs_py", bin, bin);
    //hproj->Scale(1.0 / integral);
    hprojME = histogram2ME->ProjectionY("muon-extrap/ME/dxAtMFTVsPhiAtAbs_py", bin, bin);
    hprojME->Scale(integral / integralME);
    hproj->SetLineColor(kRed);
    hproj->Draw();
    hprojME->Draw("same");
    c.SaveAs("alignment_AO2D.pdf");

    hproj->Add(hprojME, -1);
    hproj->SetLineColor(kBlack);
    TF1 fitFunc("fitFunc", "gaus(0)+pol0(3)");
    fitFunc.SetParameter(1, 0);
    fitFunc.SetParameter(2, 1);
    hproj->Fit("fitFunc", "B");
    hproj->Draw();
    c.SaveAs("alignment_AO2D.pdf");
  }

  c.Clear();
  c.SaveAs("alignment_AO2D.pdf)");

  return;
  */

  // =====================================================================================
  // DCA
  
  //plotPhiDistribution("MFT/dcaxVsPhi", -2.0, 2.0, 1, c, true);
  //plotPhiDistribution("MFT/dcayVsPhi", -2.0, 2.0, 1, c, true);

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/MFT/dcaxVsPhi");
  if (histogram2) {
    histogram2->RebinX(histogram2->GetXaxis()->GetNbins() / 4);
    histogram2->GetYaxis()->SetRangeUser(-2,2);
    //histogram2->RebinY(projRebin);
    histogram2->Draw("col");
    c.SaveAs("alignment_AO2D.pdf");
    plotMFTDCAProjection(histogram2, -180, 180, -0.5, 0.5, c, true);
  }

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/MFT/dcayVsPhi");
  if (histogram2) {
    histogram2->RebinX(histogram2->GetXaxis()->GetNbins() / 4);
    histogram2->GetYaxis()->SetRangeUser(-2,2);
    //histogram2->RebinY(projRebin);
    histogram2->Draw("col");
    c.SaveAs("alignment_AO2D.pdf");
    plotMFTDCAProjection(histogram2, -180, 180, -0.5, 0.5, c, true);
  }

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/MFT/dcaxVsy");
  if (histogram2) {
    histogram2->RebinX(10);
    histogram2->GetXaxis()->SetRangeUser(-15,15);
    //histogram2->GetYaxis()->SetRangeUser(-2,2);
    //histogram2->RebinY(projRebin);
    histogram2->Draw("col");
    c.SaveAs("alignment_AO2D.pdf");
    plotMFTDCAProjection(histogram2, -15, 15, -0.1, 0.1, c, false);
  }

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/MFT/dcayVsx");
  if (histogram2) {
    histogram2->RebinX(10);
    histogram2->GetXaxis()->SetRangeUser(-15,15);
    //histogram2->GetYaxis()->SetRangeUser(-2,2);
    //histogram2->RebinY(projRebin);
    histogram2->Draw("col");
    c.SaveAs("alignment_AO2D.pdf");
    plotMFTDCAProjection(histogram2, -15, 15, -0.1, 0.1, c, false);
  }


  plotPhiDistribution("dcaxVsPhiAtAbs", -10.0, 10.0, 4, c, true);
  plotPhiDistribution("dcayVsPhiAtAbs", -10.0, 10.0, 4, c, true);

  // =====================================================================================
  // Delta X-Y at MFT matching plane

  plotXYDistribution("dxAtMFTVsxAtAbs", -2.0, 2.0, 20, 4, c, false);
  plotXYDistribution("dxAtMFTVsyAtAbs", -2.0, 2.0, 20, 4, c, false);
  plotPhiDistribution("dxAtMFTVsPhiAtAbs", -2.0, 2.0, 4, c, true);

  plotXYDistribution("dyAtMFTVsxAtAbs", -2.0, 2.0, 10, 4, c);
  plotXYDistribution("dyAtMFTVsyAtAbs", -2.0, 2.0, 10, 4, c);
  plotPhiDistribution("dyAtMFTVsPhiAtAbs", -2.0, 2.0, 4, c, true);

  // =====================================================================================
  // Delta theta at MFT matching plane

  plotPhiDistribution("dThetaxAtMFTVsPhiAtAbs", -0.5, 0.5, 4, c, true);
  plotPhiDistribution("dThetayAtMFTVsPhiAtAbs", -0.5, 0.5, 4, c, true);

  // =====================================================================================
  // Delta X-Y at end of absorber

  plotXYDistribution("dxAtAbsVsxAtAbs", -2.0, 2.0, 10, 4, c, false);
  plotXYDistribution("dxAtAbsVsyAtAbs", -2.0, 2.0, 10, 4, c, false);
  plotPhiDistribution("dxAtAbsVsPhiAtAbs", -2.0, 2.0, 4, c, true);

  plotXYDistribution("dyAtAbsVsxAtAbs", -2.0, 2.0, 10, 4, c);
  plotXYDistribution("dyAtAbsVsyAtAbs", -2.0, 2.0, 10, 4, c);
  plotPhiDistribution("dyAtAbsVsPhiAtAbs", -2.0, 2.0, 4, c, true);

  // =====================================================================================
  // Delta theta at end of absorber

  plotXYDistribution("dThetaxAtAbsVsxAtAbs", -0.5, 0.5, 20, 4, c, true);
  plotXYDistribution("dThetaxAtAbsVsyAtAbs", -0.5, 0.5, 20, 4, c, false);
  plotPhiDistribution("dThetaxAtAbsVsPhiAtAbs", -0.5, 0.5, 4, c, true);
  plotPhiDistribution("dThetayAtAbsVsPhiAtAbs", -0.5, 0.5, 4, c, true);

  // =====================================================================================
  // Delta phi at MFT matching plane

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dPhiAtMFTVsPhiAtAbs");
  if (histogram2) {
  histogram2->RebinX(20);
  histogram2->RebinY(4);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  plotProjection(histogram2, -5, 5, 1, c, false);
  }

  // =====================================================================================
  // Delta phi at end of absorber

  histogram2 = GetTH2(fAnalysisResults, "muon-extrap/dPhiAtAbsVsPhiAtAbs");
  if (histogram2) {
  histogram2->RebinX(20);
  histogram2->RebinY(4);
  histogram2->Draw("col");
  c.SaveAs("alignment_AO2D.pdf");
  plotProjection(histogram2, -5, 5, 1, c, false);
  }

  if (true) {

  // =====================================================================================
  // MFT track / MCH cluster residuals
  std::array<int, 16> quadrants {
    100, 101, 102, 103,
    200, 201, 202, 203,
    300, 301, 302, 303,
    400, 401, 402, 403
  };
  std::array<int, 1> quadrants2 {
    100
  };

  std::unordered_map<int, TH2*> dxResidualPlots;
  std::unordered_map<int, TH2*> dyResidualPlots;

  std::array<TGraphErrors*, 16> quadrantShifts;
  for (auto& gr : quadrantShifts) {
    gr = new TGraphErrors();
  }

  std::array<std::tuple<TMultiGraph*, TGraphErrors*, TGraphErrors*>, 4> quadrantsXY;
  for (auto& graphs : quadrantsXY) {
    auto* mgr = new TMultiGraph();
    auto* grx = new TGraphErrors();
    auto* gry = new TGraphErrors();
    gry->SetLineColor(kRed);
    mgr->Add(grx);
    mgr->Add(gry);

    std::get<0>(graphs) = mgr;
    std::get<1>(graphs) = grx;
    std::get<2>(graphs) = gry;
  }


  for (int i = 0; i < quadrants.size(); i++) {
    int deId = quadrants[i];
    //plotPhiDistributionForResiduals(std::format("muon-extrap/dxVsPhiAtClusDE{}", deId).c_str(), deId, c, false);
    //plotPhiDistributionForResiduals(std::format("muon-extrap/dyVsPhiAtClusDE{}", deId).c_str(), deId, c, false);

    auto fitx = plotXYDistributionForResiduals(std::format("muon-extrap/dxVsxAtClusDE{}", deId).c_str(), deId, c, false);
    plotXYDistributionForResiduals(std::format("muon-extrap/dxVsyAtClusDE{}", deId).c_str(), deId, c, false);
    plotXYDistributionForResiduals(std::format("muon-extrap/dyVsxAtClusDE{}", deId).c_str(), deId, c, false);
    auto fity = plotXYDistributionForResiduals(std::format("muon-extrap/dyVsyAtClusDE{}", deId).c_str(), deId, c, false);

    double xMean = -fitx.GetParameter(1);
    double xError = fitx.GetParError(1);
    double yMean = -fity.GetParameter(1);
    double yError = fity.GetParError(1);
    quadrantShifts[i]->AddPoint(xMean, yMean);
    quadrantShifts[i]->SetPointError(0, xError, yError);

    int quadrantId = deId % 100;
    std::get<1>(quadrantsXY[quadrantId])->AddPoint(deId, xMean);
    std::get<2>(quadrantsXY[quadrantId])->AddPoint(deId, yMean);

    std::cout << std::format("DE{} {} {}\n", deId, xMean, yMean);
  }

  c.Clear();
  c.Divide(2, 2);
  double range = 3;
  for (int chamberIndex = 0; chamberIndex < 4; chamberIndex++) {
    for (int qIndex = 0; qIndex < 4; qIndex++) {
      c.cd(qIndex + 1);
      int graphIndex = chamberIndex * 4 + qIndex;
      quadrantShifts[graphIndex]->Draw("AEP");
      quadrantShifts[graphIndex]->GetHistogram()->GetXaxis()->SetLimits(-range, range);
      quadrantShifts[graphIndex]->SetMinimum(-range);
      quadrantShifts[graphIndex]->SetMaximum(range);

      TLine* l1 = new TLine(-range, 0, range, 0);
      l1->SetLineStyle(kDashed);
      l1->Draw();

      TLine* l2 = new TLine(0, -range, 0, range);
      l2->SetLineStyle(kDashed);
      l2->Draw();
    }
    c.SaveAs("alignment_AO2D.pdf");
  }

  for (int qIndex = 0; qIndex < 4; qIndex++) {
    c.cd(qIndex + 1);
    auto* mgr = std::get<0>(quadrantsXY[qIndex]);
    mgr->Draw("AL*");
    mgr->SetMinimum(-range);
    mgr->SetMaximum(range);
  }
  c.SaveAs("alignment_AO2D.pdf");

  }

  // =====================================================================================
  c.Clear();
  c.SaveAs("alignment_AO2D.pdf)");
}
