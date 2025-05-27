#include <MFTTracking/Constants.h>
// NOTE: not taken into account by git diff:
// 1. alignment (Andrea) -> Alignment (Chi)
// 2. mixed-event, same-event seperate on another hierarchy level c.f. Andrea c.f. Chi

TFile* fAnalysisResults;

TH1* histDxVsY;

double scaleME = 1;

constexpr int nPoints = 4;

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

int GetLR(int chamber, int deIndex)
{
  if (deIndex < 0)
    return -1;

  if (chamber >= 1 && chamber <= 4) {
    if (deIndex >= 1 && deIndex <= 2)
      return 0;
    else if (deIndex <= 3)
      return 1;
  } else if (chamber <= 6) {
    if (deIndex >= 5 && deIndex <= 13)
      return 0;
    else if (deIndex <= 17)
      return 1;
  } else {
    if (deIndex >= 7 && deIndex <= 19)
      return 0;
    else if (deIndex <= 25)
      return 1;
  }

  return -1;
}

Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

double DoubleSidedCBwithLinBgd(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] + x[0] * par[8]);
}

void PlotDCAProjection(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::string fullHistName = std::string("qa-muon/Alignment/same-event") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram2)
    return;

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName.c_str()));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  int peakMin = 50;
  // skip first bin because it is biased?
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);

    proj->Rebin(projRebin);

    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;

    int valuePeak = proj->GetMaximum();
    int binPeak = proj->GetMaximumBin();
    double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);

    if (valuePeak > peakMin) {

      if (printFits) {
        std::cout << std::format("[\"{}\"] bin {} max={}  binPeak={}  xPeak={}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);

      TF1 fcb("fcb", DoubleSidedCB, -0.2, 0.2, 7);
      fcb.SetNpx(1000);
      fcb.SetLineColor(kBlack);
      //fcb.SetParameter(0, valuePeak);
      double par[7];
      par[0]=35000;
      par[1]=xPeak;
      par[2]=0.01;
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      fcb.SetParameters(&par[0]);

      //fcb.FixParameter(1, xPeak);
      fcb.SetParLimits(2, 0.01, 0.1);
      //proj->Fit("fcb", "BRN");
      //fcb.ReleaseParameter(1);
      if (printFits)
        proj->Fit("fcb", "BR");
      else
        proj->Fit("fcb", "BRQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("residuals_AO2D.pdf");
      }

      mean = fcb.GetParameter(1);
      meanErr = fcb.GetParError(1);
      sigma = fcb.GetParameter(2);
      sigmaErr = fcb.GetParError(2);

      //mean = fgaus.GetParameter(1);
      //meanErr = fgaus.GetParError(1);
      //sigma = fgaus.GetParameter(2);
      //sigmaErr = fgaus.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);
  }

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  histogramMean->Fit("linFit", "Q");
  histogramMean->Draw("E");

  std::cout << std::format("Slope: {:0.4f} mm / 10 m", linFit.GetParameter(1) * 1000 * 10) << std::endl;
  //c.SaveAs("residuals_AO2D.pdf");
  //histogramSigma->Draw("E");
  //c.SaveAs("residuals_AO2D.pdf");

  histDxVsY = histogramMean;
}

std::pair<double, double> PlotDCAMFT(std::string histName)
{
  std::string fullHistName = std::string("qa-muon/Alignment/same-event/") + histName;
  TH2* histogram = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return {};

  histogram->GetYaxis()->SetRangeUser(-0.1, 0.1);
  histogram->Draw("col");

  return {};
}

std::pair<double, double> PlotDCAMCH(std::string histName)
{
  std::string fullHistName = std::string("qa-muon/Alignment/same-event/") + histName;
  TH1* histogram = GetTH1(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return {};
  histogram->Rebin(4);

  histogram->Draw("E");

  int valuePeak = histogram->GetMaximum();
  int binPeak = histogram->GetMaximumBin();
  double xPeak = histogram->GetXaxis()->GetBinCenter(binPeak);

  TF1 fcb("fcb", DoubleSidedCB, -6, 6, 7);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //fcb.SetParameter(0, valuePeak);
  double par[7];
  par[0]=35000;
  par[1]=xPeak;
  par[2]=0.5;
  par[3]=1;
  par[4]=1;
  par[5]=1;
  par[6]=1;
  fcb.SetParameters(&par[0]);

  fcb.FixParameter(1, xPeak);
  fcb.SetParLimits(2, 0.5, 10);
  //histogram->Fit("fcb", "BRN");
  fcb.ReleaseParameter(1);
  //histogram->Fit("fcb", "BRN");
  //histogram->Fit("fcb", "BR");

  TF1 fgaus2("fgaus2", "gaus(0)");
  fgaus2.SetNpx(1000);
  fgaus2.SetLineColor(kBlack);
  fgaus2.SetParameter(0, valuePeak / 10);
  //fgaus2.SetParLimits(0, 0, valuePeak * 10);
  fgaus2.SetParameter(1, xPeak);
  fgaus2.FixParameter(1, xPeak);
  fgaus2.SetParameter(2, 1);
  fgaus2.SetParLimits(2, 0, 10);
  histogram->Fit("fgaus2", "BQN");
  fgaus2.ReleaseParameter(1);
  fgaus2.SetParLimits(1, -10, 10);
  histogram->Fit("fgaus2", "LBQN");

  TF1 fgaus("fgaus", "gaus", fgaus2.GetParameter(1) - 2.0 * fgaus2.GetParameter(2), fgaus2.GetParameter(1) + 2.0 * fgaus2.GetParameter(2));
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  //fgaus.SetParLimits(0, 0, 100000000000.0);
  fgaus.SetParameter(0, fgaus2.GetParameter(0));
  fgaus.SetParameter(1, fgaus2.GetParameter(1));
  fgaus.SetParameter(2, fgaus2.GetParameter(2));
  fgaus.SetParLimits(2, 0, 100);
  histogram->Fit("fgaus", "RQB");

  histogram->Draw("E");

  //return {fcb.GetParameter(1), fcb.GetParError(1)};
  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
}

void PlotDXYProjection(const char* fullHistName, const char* fullHistNameME, TH2* histogram2, TH2* histogram2ME, double scaleME, float yMin, float yMax, int projRebin, TCanvas& c, bool subtractBackground, bool printFits = false)
{
  c.Clear();

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", fullHistName));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  // skip first bin because it is biased?
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName, bin), bin, bin);
    TH1* projME = histogram2ME->ProjectionY(TString::Format("%s-%d", fullHistNameME, bin), bin, bin);

    proj->Rebin(projRebin);
    projME->Rebin(projRebin);

    projME->Scale(scaleME);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > entriesMin) {

      if (printFits && subtractBackground) {
        proj->SetLineColor(kRed);
        proj->Draw("E");
        projME->Draw("L same");
        c.SaveAs("residuals_AO2D.pdf");
      }

      if (subtractBackground)
        proj->Add(projME, -1);

      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      if (printFits) {
        std::cout << std::format("[\"{}\"] bin {} max={}  binPeak={}  xPeak={}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);

      TF1 fgaus2("fgaus2", "gaus(0)+pol0(3)");
      fgaus2.SetNpx(1000);
      fgaus2.SetLineColor(kBlack);
      fgaus2.SetParameter(0, valuePeak / 10);
      //fgaus2.SetParLimits(0, 0, valuePeak * 10);
      fgaus2.SetParameter(1, xPeak);
      fgaus2.FixParameter(1, xPeak);
      fgaus2.SetParameter(2, 1);
      fgaus2.SetParLimits(2, 0, 10);
      fgaus2.SetParameter(3, 0);
      fgaus2.SetParameter(4, 0);
      proj->Fit("fgaus2", "BQN");
      fgaus2.ReleaseParameter(1);
      fgaus2.SetParLimits(1, -10, 10);
      proj->Fit("fgaus2", "BQN");

      TF1 fgaus("fgaus", "gaus", fgaus2.GetParameter(1) - 2.0 * fgaus2.GetParameter(2), fgaus2.GetParameter(1) + 2.0 * fgaus2.GetParameter(2));
      fgaus.SetNpx(1000);
      fgaus.SetLineColor(kBlack);
      //fgaus.SetParLimits(0, 0, 100000000000.0);
      fgaus.SetParameter(0, fgaus2.GetParameter(0));
      fgaus.SetParameter(1, fgaus2.GetParameter(1));
      fgaus.SetParameter(2, fgaus2.GetParameter(2));
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "RQBN");

      TF1 fcb("fcb", DoubleSidedCBwithLinBgd, -30, 30, 9);
      fcb.SetNpx(1000);
      fcb.SetLineColor(kBlack);
      //fcb.SetParameter(0, valuePeak);
      double par[9];
      par[0]=fgaus2.GetParameter(0);
      par[1]=fgaus2.GetParameter(1);
      par[2]=fgaus2.GetParameter(2);
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      par[7]=0;
      par[8]=0;
      fcb.SetParameters(&par[0]);

      fcb.SetParLimits(2, 0.5, 10);
      fcb.FixParameter(3, 1);
      fcb.FixParameter(4, 1);
      fcb.FixParameter(5, 1);
      fcb.FixParameter(6, 1);
      //histogram->Fit("fcb", "BRN");
      //fcb.ReleaseParameter(1);
      proj->Fit("fcb", "BRQN");
      fcb.ReleaseParameter(3);
      fcb.ReleaseParameter(4);
      fcb.ReleaseParameter(5);
      fcb.ReleaseParameter(6);
      if (printFits)
        proj->Fit("fcb", "BR");
      else
        proj->Fit("fcb", "BRQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("residuals_AO2D.pdf");
      }

      mean = fcb.GetParameter(1);
      meanErr = fcb.GetParError(1);
      sigma = fcb.GetParameter(2);
      sigmaErr = fcb.GetParError(2);

      //mean = fgaus.GetParameter(1);
      //meanErr = fgaus.GetParError(1);
      //sigma = fgaus.GetParameter(2);
      //sigmaErr = fgaus.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  //histogramMean->Draw("E");
  //c.SaveAs("residuals_AO2D.pdf");
  //histogramSigma->Draw("E");
  //c.SaveAs("residuals_AO2D.pdf");

  histDxVsY = histogramMean;
}

std::pair<double, double> PlotDXY(std::string histName, TCanvas& c, bool printFits = false)
{
  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/Alignment/same-event/residuals/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2) return {};

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/Alignment/mixed-event/residuals/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj = histogram2->ProjectionY((fullHistName + "_py").c_str());
  TH1* projME = histogram2ME->ProjectionY((fullHistNameME + "_py").c_str());

  double integral1 = proj->Integral(1, proj->GetXaxis()->FindBin(-15));
  double integral2 = proj->Integral(proj->GetXaxis()->FindBin(15), proj->GetXaxis()->GetNbins());
  double integralME1 = projME->Integral(1, projME->GetXaxis()->FindBin(-15));
  double integralME2 = projME->Integral(projME->GetXaxis()->FindBin(15), projME->GetXaxis()->GetNbins());
  double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

  c.cd(3);
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");

  c.cd(4);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kRed);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol0(3)");
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(0, valuePeak / 10);
  fgaus.SetParLimits(0, 0, valuePeak * 10);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 10);
  fgaus.SetParameter(3, 0);
  fgaus.SetParameter(4, 0);
  proj2->Fit("fgaus", "BQ");


  TF1 fcb("fcb","[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]");
  fcb.SetParameters(100, 0.6, -2.13903e+06, 1, xPeak);
  fcb.SetParLimits(3, 0, 100);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //proj2->Fit("fcb", "BNQ");
  //proj2->Fit("fcb", "BQ");

  proj2->GetXaxis()->SetRangeUser(-15.0, 15.0);
  proj2->Draw("E");

  c.SaveAs("residuals_AO2D.pdf");

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 1, c, true, printFits);

  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
}

std::array<std::pair<double, double>, 2> PlotDXYvsDE(std::string histName, int chamber, TCanvas& c, bool printFits = false)
{
  std::array<std::pair<double, double>, 2> result;

  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/Alignment/same-event/residuals/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2) return {};

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/Alignment/mixed-event/residuals/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj[2];
  TH1* projME[2];

  proj[0] = histogram2->ProjectionY((fullHistName + "_pyL").c_str());
  projME[0] = histogram2ME->ProjectionY((fullHistNameME + "_pyL").c_str());
  proj[1] = histogram2->ProjectionY((fullHistName + "_pyR").c_str());
  projME[1] = histogram2ME->ProjectionY((fullHistNameME + "_pyR").c_str());

  for (int lr = 0; lr < 2; lr++) {
    proj[lr]->Reset();
    proj[lr]->SetLineColor(kRed);
    projME[lr]->Reset();
  }

  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* p = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);
    TH1* pME = histogram2ME->ProjectionY(TString::Format("%s-%d", fullHistNameME.c_str(), bin), bin, bin);

    int lr = GetLR(chamber, bin - 1);
    //std::cout << std::format("CH{} DE{} LR={}", chamber, bin - 1, lr) << std::endl;
    if (lr < 0)
      continue;

    proj[lr]->Add(p);
    projME[lr]->Add(pME);
  }

  for (int lr = 0; lr < 2; lr++) {
    proj[lr]->Rebin(4);
    projME[lr]->Rebin(4);

    double integral1 = proj[lr]->Integral(1, proj[lr]->GetXaxis()->FindBin(-15));
    double integral2 = proj[lr]->Integral(proj[lr]->GetXaxis()->FindBin(15), proj[lr]->GetXaxis()->GetNbins());
    double integralME1 = projME[lr]->Integral(1, projME[lr]->GetXaxis()->FindBin(-15));
    double integralME2 = projME[lr]->Integral(projME[lr]->GetXaxis()->FindBin(15), projME[lr]->GetXaxis()->GetNbins());
    double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

    projME[lr]->Scale(scaleME);

    c.cd(3 + lr);
    TH1* proj2 = (TH1*)proj[lr]->Clone((fullHistName + "_py_" + std::to_string(lr)).c_str());
    proj2->Add(projME[lr], -1);
    proj2->SetLineColor(kRed);
    int valuePeak = proj2->GetMaximum();
    int binPeak = proj2->GetMaximumBin();
    double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
    TF1 fgaus("fgaus", "gausn(0)+pol0(3)");
    fgaus.SetNpx(1000);
    fgaus.SetLineColor(kBlack);
    fgaus.SetParameter(0, valuePeak / 10);
    fgaus.SetParLimits(0, 0, valuePeak * 10);
    fgaus.SetParameter(1, xPeak);
    fgaus.SetParameter(2, 1);
    fgaus.SetParLimits(2, 0, 10);
    fgaus.SetParameter(3, 0);
    fgaus.SetParameter(4, 0);
    proj2->Fit("fgaus", "BQN");

    TF1 fcb("fcb", DoubleSidedCBwithLinBgd, proj2->GetXaxis()->GetXmin(), proj2->GetXaxis()->GetXmax(), 9);
    fcb.SetNpx(1000);
    fcb.SetLineColor(kBlack);
    //fcb.SetParameter(0, valuePeak);
    double par[9];
    par[0]=fgaus.GetParameter(0);
    par[1]=fgaus.GetParameter(1);
    par[2]=fgaus.GetParameter(2);
    par[3]=1;
    par[4]=1;
    par[5]=1;
    par[6]=1;
    par[7]=0;
    par[8]=0;
    fcb.SetParameters(&par[0]);

    fcb.FixParameter(1, fgaus.GetParameter(1));
    fcb.SetParLimits(2, fgaus.GetParameter(2) * 0.5, fgaus.GetParameter(2) * 2);
    proj2->Fit("fcb", "BRNQ");
    fcb.ReleaseParameter(1);
    proj2->Fit("fcb", "BRQ");

    proj2->Draw("E");
    //proj[lr]->Draw("E");

    //proj[lr]->Draw("E");
    //projME[lr]->Draw("E same");

    //result[lr].first = fgaus.GetParameter(1);
    //result[lr].second = fgaus.GetParError(1);
    result[lr].first = fcb.GetParameter(1);
    result[lr].second = fcb.GetParError(1);
  }

  c.SaveAs("residuals_AO2D.pdf");

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 8, c, false, printFits);

  return result;
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
}

void PlotZTrend(int n, double* xv, std::array<std::array<std::pair<double, double>, 10>, 4>& values, std::array<std::pair<double, double>, 4>& dca, const char* title, double ymin, double ymax, TCanvas& c)
{
  double exv[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double yv[10];
  double eyv[10];

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  int colors[4] = {kBlue, kRed, kOrange, kCyan};
  int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};

  c.Clear();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  TGraphErrors* gr = new TGraphErrors();
  gr->AddPoint(0, 1000);
  gr->AddPoint(1500, 1000);
  mg->Add(gr,"lp");

  for (int j = 0; j < quadrants.size(); j++) {
    TGraphErrors* gr = new TGraphErrors();
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(markers[2]);
    gr->SetMarkerSize(2);
    gr->AddPoint(0, dca[j].first);
    mg->Add(gr,"p");
  }

  TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend->SetNColumns(4);
  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < 10; i++) {
      yv[i] = values[j][i].first;
      eyv[i] = values[j][i].second;
    }
    TGraphErrors* gr = new TGraphErrors(n, xv, yv, exv, eyv);
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(markers[1]);
    gr->SetMarkerSize(2);
    mg->Add(gr,"lp");
    auto* entry = legend->AddEntry(gr, (quadrants[j] + " ").c_str(), "P");
    entry->SetTextColor(colors[j]);
  }
  mg->Draw("a");
  mg->SetMinimum(ymin);
  mg->SetMaximum(ymax);
  legend->Draw();
  c.SaveAs("residuals_AO2D.pdf");
}

void PlotZTrendPNLR(int n, double* xv, std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>& values, std::array<std::array<std::pair<double, double>, 2>, 2>& dca, const char* title, double ymin, double ymax, TCanvas& c, bool print = false)
{
  double exv[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double yv[10];
  double eyv[10];

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};
  std::array<std::string, 2> leftRight{"MCH_left", "MCH_right"};
  std::array<std::string, 2> posNeg{"pos", "neg"};

  int colors[4] = {kBlue, kCyan, kRed, kOrange};
  int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};

  c.Clear();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  TGraphErrors* gr = new TGraphErrors();
  gr->AddPoint(0, 1000);
  gr->AddPoint(1500, 1000);
  mg->Add(gr,"lp");

  // loop on left/right
  for (int j = 0; j < 2; j++) {
    // loop on positive/negative
    for (int i = 0; i < 2; i++) {
      int colorIndex = j * 2 + i;
      TGraphErrors* gr = new TGraphErrors();
      gr->SetLineColor(colors[colorIndex]);
      gr->SetLineStyle(kDashed);
      gr->SetMarkerColor(colors[colorIndex]);
      gr->SetMarkerStyle(markers[i]);
      gr->SetMarkerSize(2);
      gr->AddPoint(0, dca[j][i].first);
      gr->AddPoint(xv[0], values[i][j][0].first);
      mg->Add(gr,"lp");
    }
  }

  TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend->SetNColumns(4);
  // loop on left/right
  for (int j = 0; j < 2; j++) {
    // loop on positive/negative
    for (int i = 0; i < 2; i++) {
      // loop on chamber
      for (int k = 0; k < 10; k++) {
        yv[k] = values[i][j][k].first;
        eyv[k] = values[i][j][k].second;
        if (print)
          std::cout << std::format("CH{} {} / {} = {:0.4f} +/- {:0.4f}", k + 1, leftRight[j], posNeg[i], yv[k], eyv[k]) << std::endl;
      }
      int colorIndex = j * 2 + i;
      TGraphErrors* gr = new TGraphErrors(n, xv, yv, exv, eyv);
      gr->SetLineColor(colors[colorIndex]);
      gr->SetMarkerColor(colors[colorIndex]);
      gr->SetMarkerStyle(markers[i]);
      gr->SetMarkerSize(2);
      mg->Add(gr,"lp");
      auto* entry = legend->AddEntry(gr, (leftRight[j] + " " + posNeg[i]).c_str(), "P");
      entry->SetTextColor(colors[colorIndex]);
    }
  }
  mg->Draw("a");
  mg->SetMinimum(ymin);
  mg->SetMaximum(ymax);
  legend->Draw();
  c.SaveAs("residuals_AO2D.pdf");
}

void plot_residuals_AO2D_Chi()
{
  // fAnalysisResults = new TFile("AnalysisResults_Chi.root");
  fAnalysisResults = new TFile("AnalysisResultsFull.root");

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};
  //std::array<std::string, 1> quadrants = {"Q0"};


  std::array<std::array<std::array<std::pair<double, double>, 2>, 2>, 2> DCAx;
  std::array<std::array<std::array<std::pair<double, double>, 2>, 2>, 2> DCAy;

  std::array<std::array<std::pair<double, double>, 10>, 4> meanDx;
  std::array<std::array<std::pair<double, double>, 10>, 4> meanDy;

  std::array<std::array<TH1*, 4>, 10> dxVsXhistograms;
  std::array<std::array<TH1*, 4>, 10> dxVsYhistograms;
  std::array<std::array<TH1*, 4>, 10> dyVsXhistograms;
  std::array<std::array<TH1*, 4>, 10> dyVsYhistograms;

  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> meanDx_LR_TB_PN;
  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> meanDy_LR_TB_PN;

  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> mchMeanDx_LR_TB_PN;
  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> mchMeanDy_LR_TB_PN;

  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dxVsDEhistograms;
  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dyVsDEhistograms;

  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dxVsPhiHistograms;
  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dyVsPhiHistograms;

  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs("residuals_AO2D.pdf(");

  c.Clear();
  c.Divide(2, 2);

  for (int q = 0; q < quadrants.size(); q++) {
    if (q == 0) c.cd(2);
    if (q == 1) c.cd(1);
    if (q == 2) c.cd(3);
    if (q == 3) c.cd(4);
    // top/bottom
    int i = (q < 2) ? 0 : 1;
    // left/right
    int j = (q == 1 || q ==2) ? 0 : 1;
    PlotDCAMFT(std::string("DCA/MFT/") + quadrants[q] + "/DCA_x_vs_z");
  }
  c.SaveAs("residuals_AO2D.pdf");

  for (int q = 0; q < quadrants.size(); q++) {
    if (q == 0) c.cd(2);
    if (q == 1) c.cd(1);
    if (q == 2) c.cd(3);
    if (q == 3) c.cd(4);
    // top/bottom
    int i = (q < 2) ? 0 : 1;
    // left/right
    int j = (q == 1 || q ==2) ? 0 : 1;
    PlotDCAProjection(std::string("DCA/MFT/") + quadrants[q] + "/DCA_x_vs_z", -0.1, 0.1, 1, c, false);
  }
  c.SaveAs("residuals_AO2D.pdf");

  for (int q = 0; q < quadrants.size(); q++) {
    if (q == 0) c.cd(2);
    if (q == 1) c.cd(1);
    if (q == 2) c.cd(3);
    if (q == 3) c.cd(4);
    // top/bottom
    int i = (q < 2) ? 0 : 1;
    // left/right
    int j = (q == 1 || q ==2) ? 0 : 1;
    PlotDCAMFT(std::string("DCA/MFT/") + quadrants[q] + "/DCA_y_vs_z");
  }
  c.SaveAs("residuals_AO2D.pdf");

  for (int q = 0; q < quadrants.size(); q++) {
    if (q == 0) c.cd(2);
    if (q == 1) c.cd(1);
    if (q == 2) c.cd(3);
    if (q == 3) c.cd(4);
    // top/bottom
    int i = (q < 2) ? 0 : 1;
    // left/right
    int j = (q == 1 || q ==2) ? 0 : 1;
    PlotDCAProjection(std::string("DCA/MFT/") + quadrants[q] + "/DCA_y_vs_z", -0.1, 0.1, 1, c, false);
  }
  c.SaveAs("residuals_AO2D.pdf");

  //return;

  for (int k = 0; k < 2; k++) {
    for (int q = 0; q < quadrants.size(); q++) {
      if (q == 0) c.cd(2);
      if (q == 1) c.cd(1);
      if (q == 2) c.cd(3);
      if (q == 3) c.cd(4);
      // top/bottom
      int i = (q < 2) ? 0 : 1;
      // left/right
      int j = (q == 1 || q ==2) ? 0 : 1;
      DCAx[i][j][k] = PlotDCAMCH(std::string("DCA/MCH/") + quadrants[q] + ((k == 0) ? "/DCA_x_pos" : "/DCA_x_neg"));
    }
    c.SaveAs("residuals_AO2D.pdf");
  }

  for (int k = 0; k < 2; k++) {
    for (int q = 0; q < quadrants.size(); q++) {
      if (q == 0) c.cd(2);
      if (q == 1) c.cd(1);
      if (q == 2) c.cd(3);
      if (q == 3) c.cd(4);
      // top/bottom
      int i = (q < 2) ? 0 : 1;
      // left/right
      int j = (q == 1 || q ==2) ? 0 : 1;
      DCAy[i][j][k] = PlotDCAMCH(std::string("DCA/MCH/") + quadrants[q] + ((k == 0) ? "/DCA_y_pos" : "/DCA_y_neg"));
    }
    c.SaveAs("residuals_AO2D.pdf");
  }

  std::array<std::string, 2> topBottom{"MFT_top", "MFT_bottom"};
  std::array<std::string, 2> posNeg{"positive", "negative"};

  for (int chamber = 0; chamber < 10; chamber++) {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (chamber == 0 && i == 0 && j == 0) print = true;
        PlotDXY(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dx_vs_phi", c, print);
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        //histDxVsY->SetTitle("TOTO");
        dxVsPhiHistograms[chamber][i][j] = histDxVsY;
      }
    }
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (i == 9 && j == 0) print = true;
        PlotDXY(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dy_vs_phi", c, print);
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        dyVsPhiHistograms[chamber][i][j] = histDxVsY;
      }
    }
  }

  for (int k = 0; k < 10; k++) {
    c.Clear();
    c.Divide(2, 2);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == 0 && j == 0) c.cd(1);
        if (i == 0 && j == 1) c.cd(2);
        if (i == 1 && j == 0) c.cd(3);
        if (i == 1 && j == 1) c.cd(4);

        dxVsPhiHistograms[k][i][j]->SetTitle(std::format("#Deltax vs. #phi - {}-{} CH{}", topBottom[i], posNeg[j], (k+1)).c_str());
        dxVsPhiHistograms[k][i][j]->GetYaxis()->SetRangeUser(-2.0, 2.0);
        dxVsPhiHistograms[k][i][j]->Draw("E");
      }
    }
    c.SaveAs("residuals_AO2D.pdf");

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == 0 && j == 0) c.cd(1);
        if (i == 0 && j == 1) c.cd(2);
        if (i == 1 && j == 0) c.cd(3);
        if (i == 1 && j == 1) c.cd(4);

        dyVsPhiHistograms[k][i][j]->SetTitle(std::format("#Deltay vs. #phi - {}-{} CH{}", topBottom[i], posNeg[j], (k+1)).c_str());
        dyVsPhiHistograms[k][i][j]->GetYaxis()->SetRangeUser(-2.0, 2.0);
        dyVsPhiHistograms[k][i][j]->Draw("E");
      }
    }
    c.SaveAs("residuals_AO2D.pdf");
  }

  for (int chamber = 0; chamber < 10; chamber++) {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (chamber == 9 && i == 0 && j == 0) print = true;
        auto result = PlotDXYvsDE(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dx_vs_de", chamber + 1, c, print);
        meanDx_LR_TB_PN[i][j][0][chamber] = result[0];
        meanDx_LR_TB_PN[i][j][1][chamber] = result[1];
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        //histDxVsY->SetTitle("TOTO");
        dxVsDEhistograms[chamber][i][j] = histDxVsY;
      }
    }
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (chamber == 9 && i == 0 && j == 0) print = true;
        if (chamber == 9 && i == 0 && j == 1) print = true;
        auto result = PlotDXYvsDE(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dy_vs_de", chamber + 1, c, print);
        meanDy_LR_TB_PN[i][j][0][chamber] = result[0];
        meanDy_LR_TB_PN[i][j][1][chamber] = result[1];
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        dyVsDEhistograms[chamber][i][j] = histDxVsY;
      }
    }
  }

  for (int k = 0; k < 10; k++) {
    c.Clear();
    c.Divide(2, 2);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == 0 && j == 0) c.cd(1);
        if (i == 0 && j == 1) c.cd(2);
        if (i == 1 && j == 0) c.cd(3);
        if (i == 1 && j == 1) c.cd(4);

        dxVsDEhistograms[k][i][j]->SetTitle(std::format("#Deltax vs. DE - {}-{} CH{}", topBottom[i], posNeg[j], (k+1)).c_str());
        dxVsDEhistograms[k][i][j]->Draw("E");
      }
    }
    c.SaveAs("residuals_AO2D.pdf");

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == 0 && j == 0) c.cd(1);
        if (i == 0 && j == 1) c.cd(2);
        if (i == 1 && j == 0) c.cd(3);
        if (i == 1 && j == 1) c.cd(4);

        dyVsDEhistograms[k][i][j]->SetTitle(std::format("#Deltay vs. DE - {}-{} CH{}", topBottom[i], posNeg[j], (k+1)).c_str());
        dyVsDEhistograms[k][i][j]->Draw("E");
      }
    }
    c.SaveAs("residuals_AO2D.pdf");
  }

/*
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < quadrants.size(); j++) {
      bool print = false;
      if (i == 9 && j == 0) print = true;
      PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dx_vs_x", c, false);
      dxVsXhistograms[i][j] = histDxVsY;
      meanDx[j][i] = PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dx_vs_y", c, print);
      dxVsYhistograms[i][j] = histDxVsY;
      meanDy[j][i] = PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dy_vs_x", c, false);
      dyVsXhistograms[i][j] = histDxVsY;
      PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dy_vs_y", c, print);
      dyVsYhistograms[i][j] = histDxVsY;
    }
  }

  for (int i = 0; i < 10; i++) {
    c.Clear();
    c.Divide(2, 2);
    for (int j = 0; j < quadrants.size(); j++) {
      if (j == 0) c.cd(2);
      if (j == 1) c.cd(1);
      if (j == 2) c.cd(3);
      if (j == 3) c.cd(4);

      dxVsYhistograms[i][j]->SetTitle(std::format("#Deltax vs. y - CH{} Q{}", i + 1, j).c_str());
      dxVsYhistograms[i][j]->Draw("E");
    }
    c.SaveAs("residuals_AO2D.pdf");

    for (int j = 0; j < quadrants.size(); j++) {
      if (j == 0) c.cd(2);
      if (j == 1) c.cd(1);
      if (j == 2) c.cd(3);
      if (j == 3) c.cd(4);

      dyVsYhistograms[i][j]->SetTitle(std::format("#Deltay vs. y - CH{} Q{}", i + 1, j).c_str());
      dyVsYhistograms[i][j]->Draw("E");
    }
    c.SaveAs("residuals_AO2D.pdf");
  }
*/

  double xv[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double defaultChamberZ[10] = {526.16, 545.24, 676.4, 695.4, 967.5,
                                998.5, 1276.5, 1307.5, 1406.6, 1437.6};

  //PlotZTrend(10, defaultChamberZ, meanDx, DCAx, "#Delta(x) vs. chamber z;chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  //PlotZTrend(10, defaultChamberZ, meanDx, DCAx, "#Delta(x) vs. chamber z;chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  //PlotZTrend(10, defaultChamberZ, meanDy, DCAy, "#Delta(y) vs. chamber z;chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  //PlotZTrend(10, defaultChamberZ, meanDy, DCAy, "#Delta(y) vs. chamber z;chamber z (cm); #Delta(y) (cm)", 0.0, 3.0, c);

  PlotZTrendPNLR(10, defaultChamberZ, meanDx_LR_TB_PN[0], DCAx[0], "#Delta(x) vs. chamber z - MFT top;chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrendPNLR(10, defaultChamberZ, meanDx_LR_TB_PN[1], DCAx[1], "#Delta(x) vs. chamber z - MFT bottom;chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrendPNLR(10, defaultChamberZ, meanDy_LR_TB_PN[0], DCAy[0], "#Delta(y) vs. chamber z - MFT top;chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c, true);
  PlotZTrendPNLR(10, defaultChamberZ, meanDy_LR_TB_PN[1], DCAy[1], "#Delta(y) vs. chamber z - MFT bottom;chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);

  /*
  // MCH residuals

  topBottom = {"MCH_top", "MCH_bottom"};
  for (int chamber = 0; chamber < 10; chamber++) {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (i == 9 && j == 0) print = true;
        auto result = PlotDXYvsDE(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dx_vs_de", chamber + 1, c, print);
        mchMeanDx_LR_TB_PN[i][j][0][chamber] = result[0];
        mchMeanDx_LR_TB_PN[i][j][1][chamber] = result[1];
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        //histDxVsY->SetTitle("TOTO");
        dxVsDEhistograms[chamber][i][j] = histDxVsY;
      }
    }
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        bool print = false;
        //if (i == 9 && j == 0) print = true;
        auto result = PlotDXYvsDE(topBottom[i] + "/" + posNeg[j] + "/CH" + std::to_string(chamber + 1) + "/dy_vs_de", chamber + 1, c, print);
        mchMeanDy_LR_TB_PN[i][j][0][chamber] = result[0];
        mchMeanDy_LR_TB_PN[i][j][1][chamber] = result[1];
        //histDxVsY->SetTitle(std::format("{} {}-{} CH{}", histDxVsY->GetTitle(), topBottom[i], posNeg[j], (k+1)).c_str());
        dyVsDEhistograms[chamber][i][j] = histDxVsY;
      }
    }
  }

  PlotZTrendPNLR(10, defaultChamberZ, mchMeanDx_LR_TB_PN[0], DCAx[0], "#Delta(x) vs. chamber z - MCH top;chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrendPNLR(10, defaultChamberZ, mchMeanDx_LR_TB_PN[1], DCAx[1], "#Delta(x) vs. chamber z - MCH bottom;chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrendPNLR(10, defaultChamberZ, mchMeanDy_LR_TB_PN[0], DCAy[0], "#Delta(y) vs. chamber z - MCH top;chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrendPNLR(10, defaultChamberZ, mchMeanDy_LR_TB_PN[1], DCAy[1], "#Delta(y) vs. chamber z - MCH bottom;chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  */

  c.Clear();
  c.SaveAs("residuals_AO2D.pdf)");

  int top = 0;
  int bottom = 1;
  int left = 0;
  int right = 1;
  int pos = 0;
  int neg = 1;
  std::cout << "\nAverage displacement at CH10:" << std::endl;
  std::cout << "* MFT top:" << std::endl;
  std::cout << std::format("    MCH left: Dx={:0.4f} Dy={:0.4f}",
      (meanDx_LR_TB_PN[top][pos][left][9].first + meanDx_LR_TB_PN[top][neg][left][9].first) / 2.0,
      (meanDy_LR_TB_PN[top][pos][left][9].first + meanDy_LR_TB_PN[top][neg][left][9].first) / 2.0) << std::endl;
  std::cout << std::format("    MCH right: Dx={:0.4f} Dy={:0.4f}",
      (meanDx_LR_TB_PN[top][pos][right][9].first + meanDx_LR_TB_PN[top][neg][right][9].first) / 2.0,
      (meanDy_LR_TB_PN[top][pos][right][9].first + meanDy_LR_TB_PN[top][neg][right][9].first) / 2.0) << std::endl;
  std::cout << "* MFT bottom:" << std::endl;
  std::cout << std::format("    MCH left: Dx={:0.4f} Dy={:0.4f}",
      (meanDx_LR_TB_PN[bottom][pos][left][9].first + meanDx_LR_TB_PN[bottom][neg][left][9].first) / 2.0,
      (meanDy_LR_TB_PN[bottom][pos][left][9].first + meanDy_LR_TB_PN[bottom][neg][left][9].first) / 2.0) << std::endl;
  std::cout << std::format("    MCH right: Dx={:0.4f} Dy={:0.4f}",
      (meanDx_LR_TB_PN[bottom][pos][right][9].first + meanDx_LR_TB_PN[bottom][neg][right][9].first) / 2.0,
      (meanDy_LR_TB_PN[bottom][pos][right][9].first + meanDy_LR_TB_PN[bottom][neg][right][9].first) / 2.0) << std::endl;

}
