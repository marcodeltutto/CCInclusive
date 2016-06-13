#include "TFile.h"

void draw_TEM_MEC(int selection = 1, int choice = 0, double normalised = false) {
  
  // choice = 0    Track Range
  // choice = 1    CosTheta
  // choice = 2    Phi
  
  
  const double BASELINE_POT = 2.300468e+20;
  const double MEC_POT      = 2.075599e+20;
  const double TEM_POT      = 2.317981e+20;
  
  const double NOMINAL_POT = 6.6e20;
  
  
  const Color_t kTotalMCColor = kRed;
  const Color_t kTotalMCErrorBandColor = kRed-10;
  const Color_t kNueSignalColor = kViolet-6;
  const Color_t kNCBackgroundColor = kBlue;
  const Color_t kBeamNueBackgroundColor = kMagenta;
  const Color_t kNumuBackgroundColor = kGreen+2;
  const Color_t kNormalHierarchyColor = kBlue;
  const Color_t kInvertedHierarchyColor = kRed;
  const Style_t k90PercentConfidenceStyle = kSolid;
  const Style_t k68PercentConfidenceStyle = 7; ///< Dashed
  
  // for Selection I
  TFile* file = new TFile("./histograms_TEM_MEC_trkrange_costheta_phi.root");
  
  // for Selection II
  TFile* fileOriginal = new TFile("./MCOriginal.root");
  TFile* fileTEM = new TFile("./MCTEM.root");
  TFile* fileMEC = new TFile("./MCMEC.root");
  
  TH1F *histoPmu, *histoPmu_TEM, *histoPmu_MEC;
  
  if (selection == 1) {
  if (choice == 0) {
    file->GetObject("Track RangeMC Prodgenie BNB Nu Cosmic",     histoPmu);
    file->GetObject("Track RangeMC Prodgenie BNB Nu Cosmic TEM", histoPmu_TEM);
    file->GetObject("Track RangeMC Prodgenie BNB Nu Cosmic MEC", histoPmu_MEC);
  } else if (choice == 1) {
    file->GetObject("cos#theta-AngleMC Prodgenie BNB Nu Cosmic",     histoPmu);
    file->GetObject("cos#theta-AngleMC Prodgenie BNB Nu Cosmic TEM", histoPmu_TEM);
    file->GetObject("cos#theta-AngleMC Prodgenie BNB Nu Cosmic MEC", histoPmu_MEC);
  } else if (choice == 2) {
    file->GetObject("#phi-AngleMC Prodgenie BNB Nu Cosmic",     histoPmu);
    file->GetObject("#phi-AngleMC Prodgenie BNB Nu Cosmic TEM", histoPmu_TEM);
    file->GetObject("#phi-AngleMC Prodgenie BNB Nu Cosmic MEC", histoPmu_MEC);
  } else
    cout << "Not a valid choice." << endl;
  } else if (selection == 2) {
    if (choice == 0) {
      fileOriginal->GetObject("SelectedFinalTrackLengthAll",     histoPmu);
      fileTEM->GetObject("SelectedFinalTrackLengthAll", histoPmu_TEM);
      fileMEC->GetObject("SelectedFinalTrackLengthAll", histoPmu_MEC);
    } else if (choice == 1) {
      fileOriginal->GetObject("SelectedFinalTrackCorrectedCosZAll",     histoPmu);
      fileTEM->GetObject("SelectedFinalTrackCorrectedCosZAll", histoPmu_TEM);
      fileMEC->GetObject("SelectedFinalTrackCorrectedCosZAll", histoPmu_MEC);
    } else if (choice == 2) {
      fileOriginal->GetObject("SelectedFinalTrackCorrectedPhiAll",     histoPmu);
      fileTEM->GetObject("SelectedFinalTrackCorrectedPhiAll", histoPmu_TEM);
      fileMEC->GetObject("SelectedFinalTrackCorrectedPhiAll", histoPmu_MEC);
    } else
      cout << "Not a valid choice." << endl;
  } else
    cout << "Not a valid choice." << endl;

    
    
  
  histoPmu->Sumw2();
  histoPmu_TEM->Sumw2();
  histoPmu_MEC->Sumw2();
  
  histoPmu->Scale(NOMINAL_POT/BASELINE_POT);
  histoPmu_TEM->Scale(NOMINAL_POT/TEM_POT);
  histoPmu_MEC->Scale(NOMINAL_POT/MEC_POT);
  
  
  
  
  
   
   double histoPmu_Int = histoPmu->Integral();
   double histoPmu_TEM_Int = histoPmu_TEM->Integral();
   double histoPmu_MEC_Int = histoPmu_MEC->Integral();
   
   
   if (normalised) {
   histoPmu->Scale(1./histoPmu_Int);
   histoPmu_TEM->Scale(1./histoPmu_TEM_Int);
   histoPmu_MEC->Scale(1./histoPmu_MEC_Int);
   }
   
   // Calculate integrals
   cout << "Integral histoPmu:    " << histoPmu->Integral() << endl;
   cout << "Integral histoPmu_TEM: " << histoPmu_TEM->Integral() << endl;
   cout << "Integral histoPmu_MEC: " << histoPmu_MEC->Integral() << endl;
   
   cout << "Difference w.r.t. histoPmu (%):" << endl;
   cout << "histoPmu_TEM: " << (histoPmu_TEM->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
   cout << "histoPmu_MEC: " << (histoPmu_MEC->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
  
  cout << endl << endl << "Difference w.r.t. histoPmu (%) -- bin by bin -- TEM" << endl;
  for (int i = 0; i < histoPmu->GetNbinsX(); i++) {
    cout << "Bin: " << i << "   " << (histoPmu_TEM->GetBinContent(i)-histoPmu->GetBinContent(i))/(histoPmu->GetBinContent(i))*100. << endl;
  }
  
  cout << endl << endl << "Difference w.r.t. histoPmu (%) -- bin by bin -- MEC" << endl;
  for (int i = 0; i < histoPmu->GetNbinsX(); i++) {
    cout << "Bin: " << i << "   " << (histoPmu_MEC->GetBinContent(i)-histoPmu->GetBinContent(i))/(histoPmu->GetBinContent(i))*100. << endl;
  }
  
  
  
  
  // Define the Canvas
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  //c->SetFillStyle(4000);      // Transparent
  //c->SetFrameFillStyle(4000); // Transparent
  
  
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
  pad1->SetRightMargin(0.05);
  //pad1->SetFillStyle(4000);       // Transparent
  //pad1->SetFrameFillStyle(4000);  // Transparent
  pad1->SetBottomMargin(0);         // Upper and lower plot are joined
  pad1->SetGridx();                 // Vertical grid
  pad1->Draw();                     // Draw the upper pad: pad1
  pad1->cd();                       // pad1 becomes the current pad
  histoPmu_TEM->SetMinimum(0.0001); // Otherwise 0 label overlaps
  histoPmu_MEC->SetMinimum(0.0001); // Otherwise 0 label overlaps
  histoPmu->SetMinimum(0.0001);     // Otherwise 0 label overlaps
  histoPmu_TEM->SetStats(0);        // No statistics on upper plot
  histoPmu_MEC->SetStats(0);        // No statistics on upper plot
  histoPmu->SetStats(0);            // No statistics on upper plot
  if (choice == 0) histoPmu_TEM->GetXaxis()->SetRangeUser(0., 700.);
  if (choice == 0 && selection == 2) histoPmu_TEM->GetYaxis()->SetRangeUser(0.0001, 16000.);
  if (choice == 1) histoPmu_TEM->GetXaxis()->SetRangeUser(-1., 1.);
  //if (choice == 1 && selection == 2) histoPmu_TEM->GetXaxis()->SetRangeUser(-1., 1.);
  if (choice == 1) histoPmu_TEM->GetYaxis()->SetRangeUser(0.0001, 14000.);
  if (choice == 1 && selection == 2) histoPmu_TEM->GetYaxis()->SetRangeUser(0.0001, 17000.);
  if (choice == 2 && selection == 1) histoPmu_TEM->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  if (choice == 2 && selection == 2) histoPmu_TEM->GetXaxis()->SetRangeUser(-3., 3.);
  if (choice == 2 && selection == 1) histoPmu_TEM->GetYaxis()->SetRangeUser(700., 1800.);
  if (choice == 2 && selection == 2) histoPmu_TEM->GetYaxis()->SetRangeUser(1000.1, 4500.);

  
  histoPmu_TEM->Draw("E2");         // Draw h1
  
  TH1F * test2 = (TH1F*)histoPmu_TEM->Clone("test2");
  test2->Draw("same histo");
  test2->SetLineColor(kBlue);

  histoPmu_MEC->Draw("same E2");    // Draw h2 on top of h1
  
  TH1F * test3 = (TH1F*)histoPmu_MEC->Clone("test3");
  test3->Draw("same histo");
  test3->SetLineColor(kGreen+2);

  histoPmu->Draw("same E2");        // Draw h3 on top of h1
  
  TH1F * test = (TH1F*)histoPmu->Clone("test");
  test->Draw("same histo");
  test->SetLineColor(kTotalMCColor);
  

  
  histoPmu_TEM->GetYaxis()->SetTitle("Selected Events");
  histoPmu_TEM->SetTitle("");
  
  if (choice == 1 && selection == 1) histoPmu_TEM->GetXaxis()->SetTitle("cos#theta");
  
  if (choice == 0 && selection == 2) histoPmu_TEM->GetXaxis()->SetTitle("Track Range [cm]");
  if (choice == 1 && selection == 2) histoPmu_TEM->GetXaxis()->SetTitle("cos#theta");
  if (choice == 2 && selection == 2) histoPmu_TEM->GetXaxis()->SetTitle("#phi angle [rad]");
  
  uBooNESimulation_2();
  
  if (normalised) {
    // TLatex
    double x = 0.87;
    double y = 0.52;
    double size = 28;
    int color = 1;
    int font = 43;
    int align = 32;
    TLatex *latex = new TLatex( x, y, "Area Normalised" );
    latex->SetNDC();
    latex->SetTextSize(size);
    latex->SetTextColor(color);
    latex->SetTextFont(font);
    latex->SetTextAlign(align);
    
    latex->Draw();
    
  }
  
  
  
  
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  //histoPmu->GetYaxis()->SetLabelSize(0.);
  //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();
  
  // Legend for the upper plot
  leg = new TLegend(0.6541353,0.5750538,0.9135338,0.8451613,NULL,"brNDC");//0.65,0.6,.85,0.87);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  //leg->SetFillStyle(0);  // Transparent
  //leg->SetHeader("");
  //leg->SetTextFont(42);
  leg->AddEntry(histoPmu, "Baseline");
  leg->AddEntry(histoPmu_TEM, "ESF+TEM");
  leg->AddEntry(histoPmu_MEC, "MEC");
  leg->Draw();
  
  
  
  /*
  TPad *padImage = new TPad("padImage", "padImage", 0.8, 0.9, 0.9, 1.0);
  padImage->Draw();                     // Draw the upper pad: pad1
  padImage->cd();
  
  TImage *img = TImage::Open("logo.png");
  
  img->SetConstRatio(0);
  img->SetImageQuality(TAttImage::kImgBest);
  
  img->Draw("xxx");
  img->SetEditable(kTRUE);
   */
  
  
  
  
  // lower plot will be in pad
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.005, 1, 0.25);
  //pad2->SetFrameFillStyle(4000); // Transparent
  //pad2->SetFillStyle(4000); // Transparent
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetRightMargin(0.05);
  pad2->SetGridx(); // vertical grid
  //pad2->SetGridy(); // orizontal grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  
  // Define the first ratio plot
  TH1F *ratio_TEM = (TH1F*)histoPmu_TEM->Clone("ratio_TEM");
  ratio_TEM->SetMinimum(0.4);  // Define Y ..
  ratio_TEM->SetMaximum(2.1); // .. range
  ratio_TEM->Sumw2();
  ratio_TEM->SetStats(0);      // No statistics on lower plot
  ratio_TEM->Divide(histoPmu);
  ratio_TEM->SetLineWidth(2);
  ratio_TEM->SetLineColor(kBlue);
  if (choice == 0) ratio_TEM->GetXaxis()->SetRangeUser(0., 700.);
  if (choice == 1) ratio_TEM->GetXaxis()->SetRangeUser(-1., 1.);
  //if (choice == 1 && selection == 2) ratio_TEM->GetXaxis()->SetRangeUser(-1., 1.);
  if (choice == 2 && selection == 1) ratio_TEM->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  if (choice == 2 && selection == 2) ratio_TEM->GetXaxis()->SetRangeUser(-3., 3.);

  ratio_TEM->Draw("E2");       // Draw the ratio plot
  
  TH1F * test4 = (TH1F*)ratio_TEM->Clone("test4");
  test4->Draw("same histo");
  test4->SetLineColor(kBlue);
  
  ratio_TEM->SetFillColor(38);
  
  
  // Define the second ratio plot
  TH1F *ratio_MEC = (TH1F*)histoPmu_MEC->Clone("ratio_MEC");
  //ratio_MEC->SetMinimum(0.8);  // Define Y ..
  //ratio_MEC->SetMaximum(1.2); // .. range
  ratio_MEC->Sumw2();
  ratio_MEC->SetStats(0);      // No statistics on lower plot
  ratio_MEC->Divide(histoPmu);
  ratio_MEC->SetLineWidth(2);
  ratio_MEC->SetLineColor(kGreen+2);
  ratio_MEC->Draw("E2 same");       // Draw the ratio plot
  
  TH1F * test5 = (TH1F*)ratio_MEC->Clone("test5");
  test5->Draw("same histo");
  test5->SetLineColor(kGreen+2);
  
  ratio_MEC->SetFillColor(29);
  
  
  
  
  
  
  //**********************
  //
  // Settings
  //
  //**********************
  
  // h1 settings
  histoPmu->SetLineColor(kRed);
  histoPmu->SetLineWidth(2);
  histoPmu->SetFillColor(kTotalMCErrorBandColor);
  if(choice == 2) histoPmu->SetFillColorAlpha(kTotalMCErrorBandColor, 0.35);//histoPmu->SetFillStyle(3001);
  
  // Y axis h1 plot settings
  histoPmu_MEC->GetYaxis()->CenterTitle();
  histoPmu_MEC->GetYaxis()->SetTitleSize(25);
  histoPmu_MEC->GetYaxis()->SetTitleFont(43);
  histoPmu_MEC->GetYaxis()->SetTitleOffset(1.55);
  
  // h2 settings
  histoPmu_TEM->SetLineColor(kBlue);
  histoPmu_TEM->SetLineWidth(2);
  histoPmu_TEM->SetFillColor(38);
  histoPmu_TEM->GetYaxis()->CenterTitle();
  histoPmu_TEM->GetYaxis()->SetTitleSize(25);
  histoPmu_TEM->GetYaxis()->SetTitleFont(43);
  histoPmu_TEM->GetYaxis()->SetTitleOffset(1.55);
  //histoPmu_TEM->SetFillStyle(3001);
  
  
  // h3 settings
  histoPmu_MEC->SetLineColor(kGreen+2);
  histoPmu_MEC->SetLineWidth(2);
  //histoPmu_MEC->SetFillColor(29);
  //histoPmu_MEC->SetFillStyle(3001);
  histoPmu_MEC->SetFillColorAlpha(29, 0.60);
  
  // Ratio plot (ratio_TEM) settings
  ratio_TEM->SetTitle(""); // Remove the ratio title
  
  // Y axis ratio plot settings
  ratio_TEM->GetYaxis()->SetTitle("Ratio");
  ratio_TEM->GetYaxis()->CenterTitle();
  ratio_TEM->GetYaxis()->SetNdivisions(505);
  ratio_TEM->GetYaxis()->SetTitleSize(25);
  ratio_TEM->GetYaxis()->SetTitleFont(43);
  ratio_TEM->GetYaxis()->SetTitleOffset(1.0);
  ratio_TEM->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio_TEM->GetYaxis()->SetLabelSize(15);
  
  // X axis ratio plot settings
  ratio_TEM->GetXaxis()->CenterTitle();
  ratio_TEM->GetXaxis()->SetTitleSize(25);
  ratio_TEM->GetXaxis()->SetTitleFont(43);
  ratio_TEM->GetXaxis()->SetTitleOffset(3.5);
  ratio_TEM->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio_TEM->GetXaxis()->SetLabelSize(20);
  
  //ratio_MEC->SetFillStyle(3001);
  ratio_MEC->SetFillColorAlpha(29, 0.60);
  
  // Draw linea at 1 in ratio plot
  TLine *line;
  if (choice == 0 && selection == 1) line = new TLine(0,1,726,1);
  if (choice == 0 && selection == 2) line = new TLine(0,1,700,1);
  if (choice == 1) line = new TLine(-1,1,1,1);
  if (choice == 2) line = new TLine(-TMath::Pi(),1,TMath::Pi(),1);
  if (choice == 2 && selection == 2) line = new TLine(-TMath::Pi(),1,TMath::Pi(),1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(9); // dashed
  line->Draw();
  
  
  
  
  //**********************
  //
  // Area Normalised plot
  //
  //**********************
  /*
   TH1F *histoPmu_copy = (TH1F*)histoPmu->Clone("histoPmu_copy");
   TH1F *histoPmu_TEM_copy = (TH1F*)histoPmu_TEM->Clone("histoPmu_TEM_copy");
   TH1F *histoPmu_MEC_copy = (TH1F*)histoPmu_MEC->Clone("histoPmu_MEC_copy");
   
   // Area Norm
   histoPmu_copy->Scale(1./histoPmu_copy->Integral());
   histoPmu_TEM_copy->Scale(1./histoPmu_TEM_copy->Integral());
   histoPmu_MEC_copy->Scale(1./histoPmu_MEC_copy->Integral());
   
   // Draw them in a new canvas
   TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
   histoPmu_TEM_copy->Draw();
   histoPmu_copy->Draw("same");
   histoPmu_MEC_copy->Draw("same");
   
   //Settings
   histoPmu_TEM_copy->GetXaxis()->CenterTitle();
   histoPmu_TEM_copy->GetYaxis()->CenterTitle();
   //histoPmu_TEM_copy->GetYaxis()->SetTitleOffset(1.0);
   
   leg->Draw();
   
   // TLatex
   double x = 0.87;
   double y = 0.52;
   double size = 28;
   int color = 1;
   int font = 43;
   int align = 32;
   
   TLatex *latex = new TLatex( x, y, "Area Normalised" );
   latex->SetNDC();
   latex->SetTextSize(size);
   latex->SetTextColor(color);
   latex->SetTextFont(font);
   latex->SetTextAlign(align);
   
   latex->Draw();
   */
  
  
}
