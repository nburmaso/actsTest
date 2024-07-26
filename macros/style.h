void SetPad(TVirtualPad* pad){
  pad->SetRightMargin(0.015);
  pad->SetLeftMargin(0.12);
  pad->SetBottomMargin(0.12);
  pad->SetTopMargin(0.06);
}

void SetHisto(TH1* h, TString title){
  h->SetTitleOffset(1.1);
  h->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->SetTitle(title.Data());
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
}