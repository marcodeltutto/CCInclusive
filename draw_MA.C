#include "TFile.h"

void draw_MA(int selection = 1, int choice = 0, double normalised = false) {
  
  // choice = 0    Track Range
  // choice = 1    CosTheta
  // choice = 2    Phi
  
  // selection = 1      Christoph selection
  // selection = 2      Xiao selection
  
  
  const double BASELINE_POT = 2.300468e+20;
  const double MA_POT      = 2.000320e+20;
  
  const double NOMINAL_POT = 6.6e20;
  
  
  // Define colours...but I am not using them...
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
  TFile* file = new TFile("./histograms_MA_trkrange_costheta_phi.root");

  // for Selection II
  TFile* fileOriginal = new TFile("./MCOriginal.root");
  TFile* fileM_A = new TFile("./MCM_A.root");

  TH1F *histoPmu, *histoPmu_MA;
  
  if (selection == 1) {
    if (choice == 0) {
      file->GetObject("Track RangeMC Prodgenie BNB Nu Cosmic",     histoPmu);
      file->GetObject("Track RangeMC Prodgenie BNB Nu Cosmic M_A", histoPmu_MA);
    } else if (choice == 1) {
      file->GetObject("cos#theta-AngleMC Prodgenie BNB Nu Cosmic",     histoPmu);
      file->GetObject("cos#theta-AngleMC Prodgenie BNB Nu Cosmic M_A", histoPmu_MA);
    } else if (choice == 2) {
      file->GetObject("#phi-AngleMC Prodgenie BNB Nu Cosmic",     histoPmu);
      file->GetObject("#phi-AngleMC Prodgenie BNB Nu Cosmic M_A", histoPmu_MA);
    } else
      cout << "Not a valid choice." << endl;
  } else if (selection == 2) {
    if (choice == 0) {
      fileOriginal->GetObject("SelectedFinalTrackLengthAll",     histoPmu);
      fileM_A->GetObject("SelectedFinalTrackLengthAll", histoPmu_MA);
    } else if (choice == 1) {
      fileOriginal->GetObject("SelectedFinalTrackCorrectedCosZAll",     histoPmu);
      fileM_A->GetObject("SelectedFinalTrackCorrectedCosZAll", histoPmu_MA);
    } else if (choice == 2) {
      fileOriginal->GetObject("SelectedFinalTrackCorrectedPhiAll",     histoPmu);
      fileM_A->GetObject("SelectedFinalTrackCorrectedPhiAll", histoPmu_MA);
    } else
      cout << "Not a valid choice." << endl;
  } else
    cout << "Not a valid choice." << endl;
  
  
  histoPmu->Sumw2();     // just to be sure
  histoPmu_MA->Sumw2();  // just to be sure
  
  // POt scaling
  histoPmu->Scale(NOMINAL_POT/BASELINE_POT);
  histoPmu_MA->Scale(NOMINAL_POT/MA_POT);
  
  
  double histoPmu_Int = histoPmu->Integral();
  double histoPmu_MA_Int = histoPmu_MA->Integral();
  
  
  if (normalised) {
    histoPmu->Scale(1./histoPmu_Int);
    histoPmu_MA->Scale(1./histoPmu_MA_Int);
  }
  
  
  // Calculate integrals and print differences
  cout << "Integral histoPmu:    " << histoPmu->Integral() << endl;
  cout << "Integral histoPmu_MA: " << histoPmu_MA->Integral() << endl;
  
  cout << "Difference w.r.t. histoPmu (%):" << endl;
  cout << "histoPmu_MA: " << (histoPmu_MA->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
  
  cout << endl << endl << "Difference w.r.t. histoPmu (%) -- bin by bin -- MA" << endl;
  for (int i = 0; i < histoPmu->GetNbinsX(); i++) {
    cout << "Bin: " << i << "   " << (histoPmu_MA->GetBinContent(i)-histoPmu->GetBinContent(i))/(histoPmu->GetBinContent(i))*100. << endl;
  }
  
  
  
  
  
  
  // Define the Canvas
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  //c->SetFillStyle(4000);      // Transparent
  //c->SetFrameFillStyle(4000); // Transparent
  
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
  //pad1->SetFillStyle(4000);       // Transparent
  //pad1->SetFrameFillStyle(4000);  // Transparent
  pad1->SetBottomMargin(0);         // Upper and lower plot are joined
  pad1->SetRightMargin(0.05);
  pad1->SetGridx();                 // Vertical grid
  pad1->Draw();                     // Draw the upper pad: pad1
  pad1->cd();                       // pad1 becomes the current pad
  histoPmu_MA->SetMinimum(0.0001);  // Otherwise 0 label overlaps
  histoPmu->SetMinimum(0.0001);     // Otherwise 0 label overlaps
  histoPmu_MA->SetStats(0);         // No statistics on upper plot
  histoPmu->SetStats(0);            // No statistics on upper plot
  
  // Define plots range based on the quantity plotted
  if (choice == 0) histoPmu_MA->GetXaxis()->SetRangeUser(0., 700.);
  if (choice == 1 && selection == 1) histoPmu_MA->GetXaxis()->SetRangeUser(-1., 1.);
  if (choice == 1 && selection == 2) histoPmu_MA->GetXaxis()->SetRangeUser(-1., 1.);
  if (choice == 2 && selection == 1) histoPmu_MA->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  if (choice == 2 && selection == 2) histoPmu_MA->GetXaxis()->SetRangeUser(-3., 3.);
  if (choice == 1 && selection == 1) histoPmu_MA->GetYaxis()->SetRangeUser(0.0001, 13000.);
  if (choice == 1 && selection == 2) histoPmu_MA->GetYaxis()->SetRangeUser(0.0001, 15000.);
  if (choice == 2 && selection == 1) histoPmu_MA->GetYaxis()->SetRangeUser(700., 1800.);
  if (choice == 2 && selection == 2) histoPmu_MA->GetYaxis()->SetRangeUser(1000.1, 4500.);

  
  
  // Draw!
  histoPmu_MA->Draw("E2");                                 // Draw error bars only
  TH1F * test2 = (TH1F*)histoPmu_MA->Clone("test2");
  test2->SetLineColor(kGreen+2);
  test2->Draw("same histo");                               // Draw the histo line now
  
  histoPmu->Draw("same E2");                               // Draw error bars only
  TH1F * test = (TH1F*)histoPmu->Clone("test");
  test->Draw("same histo");
  test->SetLineColor(kTotalMCColor);                       // Draw the histo line now

  
  // Change titles
  histoPmu_MA->GetYaxis()->SetTitle("Selected Events");
  histoPmu_MA->SetTitle("");
  if (choice == 1 && selection == 1) histoPmu_MA->GetXaxis()->SetTitle("cos#theta");
  if (choice == 2 && selection == 1) histoPmu_MA->GetXaxis()->SetTitle("#phi angle [rad]");
  if (choice == 0 && selection == 2) histoPmu_MA->GetXaxis()->SetTitle("Track Range [cm]");
  if (choice == 1 && selection == 2) histoPmu_MA->GetXaxis()->SetTitle("cos#theta");
  if (choice == 2 && selection == 2) histoPmu_MA->GetXaxis()->SetTitle("#phi angle [rad]");

  uBooNESimulation_2();  // Simulation tag, this is defined in /nashome/m/mdeltutt/rootlogon.C
  
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
  // axis instead, in order to avoid the first label (0) to be clipped.   (THIS SHOULDN'T BE NECESSARY NOW)
  //histoPmu->GetYaxis()->SetLabelSize(0.);
  //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();
  
  // Legend for the upper plot
  leg = new TLegend(0.5964912,0.6266667,0.8922306,0.8348387,NULL,"brNDC");//0.65,0.6,.85,0.87);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  //leg->SetFillStyle(0);  // Transparent
  //leg->SetHeader("");
  leg->AddEntry(histoPmu,    "Baseline");
  leg->AddEntry(histoPmu_MA, "M_{A}^{CCQE} = 1.35 GeV");
  leg->Draw();
  
  
  
  
  
  
  
  
  // LOWER plot will be in pad
  c->cd();                            // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.005, 1, 0.25);
  //pad2->SetFrameFillStyle(4000);    // Transparent
  //pad2->SetFillStyle(4000);         // Transparent
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);         // Leave some space for the X axis title
  pad2->SetRightMargin(0.05);
  pad2->SetGridx();                   // vertical grid
  pad2->Draw();
  pad2->cd();                         // pad2 becomes the current pad
  
  // Define the first ratio plot
  TH1F *ratio_MA = (TH1F*)histoPmu_MA->Clone("ratio_MA");
  ratio_MA->SetMinimum(0.4);   // Define Y ..
  ratio_MA->SetMaximum(2.1);   // .. range
  ratio_MA->Sumw2();
  ratio_MA->SetStats(0);       // No statistics on lower plot
  ratio_MA->Divide(histoPmu);
  ratio_MA->SetLineWidth(2);
  ratio_MA->SetLineColor(kGreen+2);
  if (choice == 0) ratio_MA->GetXaxis()->SetRangeUser(0., 700.);
  if (choice == 1) ratio_MA->GetXaxis()->SetRangeUser(-1., 1.);
  if (choice == 2 && selection == 1) ratio_MA->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  if (choice == 2 && selection == 2) ratio_MA->GetXaxis()->SetRangeUser(-3., 3.);

  // Draw!
  ratio_MA->Draw("E2");                              // Draw error bars only
  TH1F * test4 = (TH1F*)ratio_MA->Clone("test4");
  test4->SetLineColor(kGreen+2);
  test4->Draw("same histo");                         // Draw the histo line now
  
  
  
  
  
  
  
  
  //**********************
  //
  // Settings
  //
  //**********************
  
  // h1 settings
  histoPmu->SetLineColor(kRed);
  histoPmu->SetLineWidth(2);
  histoPmu->SetFillColor(kTotalMCErrorBandColor);
  
  
  // h2 settings
  histoPmu_MA->SetLineColor(kGreen+2);
  histoPmu_MA->SetLineWidth(2);
  histoPmu_MA->SetFillColor(29);
  histoPmu_MA->GetYaxis()->CenterTitle();
  histoPmu_MA->GetYaxis()->SetTitleSize(25);
  histoPmu_MA->GetYaxis()->SetTitleFont(43);
  histoPmu_MA->GetYaxis()->SetTitleOffset(1.55);
  
  
  // Ratio plot (ratio_MA) settings
  ratio_MA->SetTitle(""); // Remove the ratio title
  ratio_MA->SetFillColor(29);
  
  // Y axis ratio plot settings
  ratio_MA->GetYaxis()->SetTitle("Ratio");
  ratio_MA->GetYaxis()->CenterTitle();
  ratio_MA->GetYaxis()->SetNdivisions(505);
  ratio_MA->GetYaxis()->SetTitleSize(25);
  ratio_MA->GetYaxis()->SetTitleFont(43);
  ratio_MA->GetYaxis()->SetTitleOffset(1.0);
  ratio_MA->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio_MA->GetYaxis()->SetLabelSize(15);
  
  // X axis ratio plot settings
  ratio_MA->GetXaxis()->CenterTitle();
  ratio_MA->GetXaxis()->SetTitleSize(25);
  ratio_MA->GetXaxis()->SetTitleFont(43);
  ratio_MA->GetXaxis()->SetTitleOffset(3.5);
  ratio_MA->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio_MA->GetXaxis()->SetLabelSize(20);
  
  
  
  
  // Draw line at 1 in ratio plot
  TLine *line;
  if (choice == 0 && selection == 1) line = new TLine(0,1,726,1);
  if (choice == 0 && selection == 2) line = new TLine(0,1,700,1);
  if (choice == 1) line = new TLine(-1,1,1,1);
  if (choice == 2) line = new TLine(-TMath::Pi(),1,TMath::Pi(),1);
  if (choice == 2 && selection == 2) line = new TLine(-TMath::Pi(),1,TMath::Pi(),1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(9); // dashed
  line->Draw();
  
  
  

  
}
