// compareRootFiles.C

TCanvas* GetCanvas_Primitives(TString fName, TString cName) {
    TFile *fInput = new TFile(fName.Data());
    TCanvas *c = (TCanvas*)fInput->Get(cName.Data());
    TCanvas *c_new = new TCanvas(Form("copy_%s_%s", cName.Data(), fName.Data()), fName, 1200, 800);
    TList* prims = c->GetListOfPrimitives();
    for (TObject* obj : *prims) {
        TObject* copy = obj->Clone();
        copy->Draw();  // redraw on new canvas
    }
    return c_new;
}

TCanvas* GetCanvas(TString fName, TString cName) {
    TFile *fInput = new TFile(fName.Data());
    TCanvas *c = (TCanvas*)fInput->Get(cName.Data());
    return c;
}

Int_t compareRootFiles() {
    std::vector<TString> vfNames;

    vfNames.push_back("output_residuals_LHC24an_pass1_skimmed_small.root");
    vfNames.push_back("output_residuals_LHC24aq_pass1_small.root");

    TCanvas dummy;
    dummy.SaveAs("comparison_summary.pdf(");

    for (Int_t i = 0; i < vfNames.size(); i++) {
        TString fName = vfNames[i];        
    
        // MFT DCA
        TCanvas *c_MFT_DCA_x = GetCanvas_Primitives(fName, "c;1");
        TCanvas *c_MFT_DCA_y = GetCanvas_Primitives(fName, "c;2");

        c_MFT_DCA_x->SaveAs("comparison_summary.pdf");
        c_MFT_DCA_y->SaveAs("comparison_summary.pdf");
    }

    TCanvas *c_test = GetCanvas(vfNames[0], "c;3");
    c_test->Draw();

    /*
    for (Int_t i = 0; i < 4; i++) { // dx_MFT_top, dx_MFT_bottom, dy_MFT_top, dy_MFT_bottom
        TCanvas c("c", "c", 1200, 800);
        c.Divide(2,2);
        for (Int_t j = 0; j < vfNames.size(); j++) {
            TString fName = vfNames[j];

            // PlotZTrendPNLR
            TCanvas *c_ZTrend = GetCanvas_Primitives(fName, Form("c;%i", i + 3));

            c.cd(1+j);
            TList* prims = c_ZTrend->GetListOfPrimitives();
            for (TObject* obj : *prims) {
                TObject* copy = obj->Clone();
                copy->Draw();
            }

            c_ZTrend->SaveAs("comparison_summary.pdf");
            c_ZTrend->Draw();
        }
        c.SaveAs("comparison_summary.pdf");
        c.Clear();
    }
    */

    dummy.SaveAs("comparison_summary.pdf)");

    return 0;
}