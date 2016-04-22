#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>

double beammin = 3.55-0.36; //us. Beam window start
double beammax = 5.15-0.36; //us. Beam window end

float GetMaximum(const std::vector<TH1F*>& HistVector);
void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight);

void HistoProducer()
{
    TGaxis::SetMaxDigits(2);
    std::vector<TChain*> ChainVec;

    std::vector<std::string> FileNameVec;
    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;

    std::vector<TH1F*> FlashTime;
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionPhi;
    std::vector<TH1F*> SelXVtxPosition;
    std::vector<TH1F*> SelYVtxPosition;
    std::vector<TH1F*> SelZVtxPosition;

    TH1F* SelectionCCTheta;
    TH1F* SelectionNCTheta;
    TH1F* SelectionCCPhi;
    TH1F* SelectionNCPhi;
    TH1F* SelectionCCRange;
    TH1F* SelectionNCRange;

    TLegend* LegendData = new TLegend(0.7,0.8,0.9,0.9);
    LegendData->SetHeader("Data Type");

    DataLabel.push_back("On-Beam");
    DataLabel.push_back("Off-Beam");

    TLegend* LegendMC = new TLegend(0.65,0.8,0.9,0.9);
//     LegendMC->SetHeader("Data Type");

    MCLabel.push_back("On-Beam Minus Off-Beam Sample");
    MCLabel.push_back("MC BNB Only with Stat. Error");
    MCLabel.push_back("Expected NC Contamination");

    TLegend* FlashLabel = new TLegend(0.7,0.7,0.9,0.9);
//     FlashLabel->SetHeader("Generator Type");

    GenLabel.push_back("Data On-Beam BNB");
    GenLabel.push_back("Data Off-Beam BNBEXT");
    GenLabel.push_back("MC Prodgenie BNB Nu");
    GenLabel.push_back("MC Prodgenie BNB Nu Cosmic");
    GenLabel.push_back("MC Prodcosmic Corsika in-Time");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
    ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");

    ChainVec.push_back(new TChain("analysistree/anatree"));
    ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1.root");
    ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2.root");

    std::vector<double> PotBNB(ChainVec.size());
    std::vector<int> TrackCandidate(ChainVec.size());
    std::vector<int> VertexCandidate(ChainVec.size());

    std::vector<std::vector<float>> TrackTheta(ChainVec.size());

    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        ChainVec.at(file_no)->SetBranchAddress("potbnb", &PotBNB.at(file_no));
        ChainVec.at(file_no)->SetBranchAddress("TrackCand", &TrackCandidate.at(file_no));
        ChainVec.at(file_no)->SetBranchAddress("VertexCand", &VertexCandidate.at(file_no));

        ChainVec.at(file_no)->SetBranchAddress("trktheta_pandoraNu", &TrackTheta.at(file_no));

        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no)->GetEntries(); tree_index++)
        {
            ChainVec.at(file_no)->GetEntry(tree_index);
        }


//             treenc -> SetBranchAddress("flash_zcenter", flash_zcenter);
//             treenc -> SetBranchAddress("flash_ycenter", flash_ycenter);
//             treenc -> SetBranchAddress("mcevts_truth", &mcevts_truth);
//             treenc -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
//             treenc -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
//             treenc -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
//             treenc -> SetBranchAddress("ccnc_truth", ccnc_truth);
//             treenc -> SetBranchAddress("nuPDG_truth", nuPDG_truth);
//             treenc -> SetBranchAddress("pdg", PDG_truth);
//             treenc -> SetBranchAddress("mode_truth", mode_truth);
//             treenc -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
//             treenc -> SetBranchAddress("StartPointx", XMCTrackStart);
//             treenc -> SetBranchAddress("StartPointy", YMCTrackStart);
//             treenc -> SetBranchAddress("StartPointz", ZMCTrackStart);
//             treenc -> SetBranchAddress("EndPointx", XMCTrackEnd);
//             treenc -> SetBranchAddress("EndPointy", YMCTrackEnd);
//             treenc -> SetBranchAddress("EndPointz", ZMCTrackEnd);
//
//             // Product specific stuff
//             treenc -> SetBranchAddress(("ntracks_"+TrackingName).c_str(),&ntracks_reco);
//             treenc -> SetBranchAddress(("trkstartx_"+TrackingName).c_str(),trkstartx);
//             treenc -> SetBranchAddress(("trkstarty_"+TrackingName).c_str(),trkstarty);
//             treenc -> SetBranchAddress(("trkstartz_"+TrackingName).c_str(),trkstartz);
//             treenc -> SetBranchAddress(("trkendx_"+TrackingName).c_str(),trkendx);
//             treenc -> SetBranchAddress(("trkendy_"+TrackingName).c_str(),trkendy);
//             treenc -> SetBranchAddress(("trkendz_"+TrackingName).c_str(),trkendz);
//             treenc -> SetBranchAddress(("trktheta_"+TrackingName).c_str(),trktheta);
//             treenc -> SetBranchAddress(("trkphi_"+TrackingName).c_str(),trkphi);
//
//                 treenc -> SetBranchAddress(("nvtx_"+VertexingName).c_str(), &nvtx);
//                 treenc -> SetBranchAddress(("vtxx_"+VertexingName).c_str(), vtxx);
//                 treenc -> SetBranchAddress(("vtxy_"+VertexingName).c_str(), vtxy);
//                 treenc -> SetBranchAddress(("vtxz_"+VertexingName).c_str(), vtxz);
    }
    
    


//     FileNameVec.push_back("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_v05_08_00.root");
//     FileNameVec.push_back("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_v05_08_00.root");
//     FileNameVec.push_back("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_prodcosmics_corsika_inTime_v05_08_00.root");
//
//     TH1F *hSelectionTheta = new TH1F("#theta-Angle of Selected Track","#theta-Angle of Selected Track",90,0,3.142);
//     SelectionTheta
//     hSelectionTheta->SetStats(0);
//     hSelectionTheta->GetXaxis()->SetTitle("#theta angle [rad]");
//     hSelectionTheta->GetYaxis()->SetTitle("Number of Tracks [ ]");
//
//     TLine FlashTimeMinCut(beammin, 0, beammin, 0.6);
//     TLine FlashTimeMaxCut(beammax, 0, beammax, 0.6);
//     FlashTimeMinCut.SetLineColor(1);
//     FlashTimeMaxCut.SetLineColor(1);
//     TLine FlashTimeMinCut_1(beammin+0.36, 0, beammin+0.36, 0.6);
//     TLine FlashTimeMaxCut_1(beammax+0.36, 0, beammax+0.36, 0.6);
//     FlashTimeMinCut_1.SetLineColor(2);
//     FlashTimeMaxCut_1.SetLineColor(2);
//
//     TF1* Const = new TF1("const","0.1",0,8.2);
//     TF1* SinTheta = new TF1("const","sin(x)",0,3.142);
//
//     unsigned int FileNumber = 0;
//
//     for(const auto& FileName : FileNameVec)
//     {
//         TFile* File = new TFile( (FileName).c_str(),"READ" ) ;
//         File->cd();
//
//         FlashTime.push_back((TH1F*)File->Get("Flash Time"));
//         SelectionTrackRange.push_back((TH1F*)File->Get("Track Range of Selected Track"));
//         SelectionTheta.push_back((TH1F*)File->Get("#theta-Angle of Selected Track"));
//         SelectionPhi.push_back((TH1F*)File->Get("#phi-Angle of Selected Track"));
//
//         if(FileNumber == 3)
//         {
//             SelectionCCTheta = (TH1F*)File->Get("#theta-Angle of Selected CC Track");
//             SelectionNCTheta = (TH1F*)File->Get("#theta-Angle of Selected NC Track");
//             SelectionCCPhi = (TH1F*)File->Get("#phi-Angle of Selected CC Track");
//             SelectionNCPhi = (TH1F*)File->Get("#phi-Angle of Selected NC Track");
//             SelectionCCRange = (TH1F*)File->Get("Track Range of Selected CC Track");
//             SelectionNCRange = (TH1F*)File->Get("Track Range of Selected NC Track");
//         }
//
//     FileNumber++;
// }
//
//
// for(unsigned int hist_no = 0; hist_no < SelectionTrackRange.size(); hist_no++)
// {
//     SelectionTrackRange.at(hist_no)->Sumw2();
//     SelectionTheta.at(hist_no)->Sumw2();
//     if(hist_no != 1)SelectionPhi.at(hist_no)->Sumw2();
// }
//
// AddFirstTwoHistograms(FlashTime,1.);
// AddFirstTwoHistograms(SelectionTrackRange,1.);
// AddFirstTwoHistograms(SelectionTheta,1.);
// SelectionPhi.erase(SelectionPhi.begin()+1);
// //     AddFirstTwoHistograms(SelectionPhi,1.);
//
// for(unsigned int hist_no = 0; hist_no < FlashTime.size(); hist_no++)
// {
//     FlashTime.at(hist_no)->Rebin(5);
//     SelectionTrackRange.at(hist_no)->Rebin(5);
//     SelectionTheta.at(hist_no)->Rebin(5);
//     SelectionPhi.at(hist_no)->Rebin(10);
//
//     FlashLabel->AddEntry( FlashTime.at(hist_no), (GenLabel.at(hist_no)).c_str(),"l" );
//     if(hist_no < DataLabel.size())
//     {
//         LegendData->AddEntry( SelectionTrackRange.at(hist_no), (DataLabel.at(hist_no)).c_str(),"l" );
//     }
// }
//
// SelectionCCTheta->Rebin(5);
// SelectionNCTheta->Rebin(5);
// SelectionCCPhi->Rebin(10);
// SelectionNCPhi->Rebin(10);
// SelectionCCRange->Rebin(5);
// SelectionNCRange->Rebin(5);
//
// SelectionCCTheta->Scale(1/SelectionTheta.at(2)->Integral());
// SelectionNCTheta->Scale(1/SelectionTheta.at(2)->Integral());
// SelectionCCPhi->Scale(1/SelectionPhi.at(2)->Integral());
// SelectionNCPhi->Scale(1/SelectionPhi.at(2)->Integral());
// SelectionCCRange->Scale(1/SelectionTrackRange.at(2)->Integral());
// SelectionNCRange->Scale(1/SelectionTrackRange.at(2)->Integral());
//
// for(unsigned int index = 0; index < 2; index++)
// {
//     SelectionTrackRange.pop_back();
//     SelectionTheta.pop_back();
//     SelectionPhi.pop_back();
// }
//
// //     FlashTime.pop_back();
//
// FlashTime.front()->Scale(3./FlashTime.front()->Integral());
// FlashTime.at(1)->Scale(3./FlashTime.at(1)->Integral());
// FlashTime.at(2)->Scale(2./FlashTime.at(2)->Integral());
// FlashTime.at(3)->Scale(2./FlashTime.at(3)->Integral());
// FlashTime.at(4)->Scale(2./FlashTime.at(4)->Integral());
// SelectionTrackRange.front()->Scale(1./567157.);
// SelectionTrackRange.at(1)->Scale(1./130967);
// SelectionTrackRange.at(2)->Scale(1./10880.);
// SelectionTheta.front()->Scale(1./567157.);
// SelectionTheta.at(1)->Scale(1./130967.);
// SelectionTheta.at(2)->Scale(1./10880.);
// SelectionPhi.front()->Scale(1./300000.);
// SelectionPhi.at(1)->Scale(1./130967.);
// SelectionPhi.at(2)->Scale(1./10880.);
//
// for(unsigned int hist_no = 0; hist_no < FlashTime.size(); hist_no++)
// {
//     FlashTime.at(hist_no)->Add(Const,FlashTime.size()-hist_no-1);
// }
//
// TCanvas *Canvas1 = new TCanvas("Flash Time", "Flash Time", 1400, 1000);
// Canvas1->cd();
// FlashTime.front()->SetMaximum(1.1*GetMaximum(FlashTime));
// FlashTime.front()->SetMinimum(0);
// FlashTime.front()->GetYaxis()->SetTitle("#propto # Flashes");
// FlashTime.front()->Draw();
// FlashTime.at(1)->SetLineColor(2);
// FlashTime.at(1)->Draw("SAME");
// FlashTime.at(2)->SetLineColor(3);
// FlashTime.at(2)->Draw("SAME");
// FlashTime.at(3)->SetLineColor(6);
// FlashTime.at(3)->Draw("SAME");
// FlashTime.at(4)->SetLineColor(7);
// FlashTime.at(4)->Draw("SAME");
// FlashTimeMinCut.Draw("SAME");
// FlashTimeMaxCut.Draw("SAME");
// FlashTimeMinCut_1.Draw("SAME");
// FlashTimeMaxCut_1.Draw("SAME");
// FlashLabel->Draw();
// Canvas1->SaveAs("DataFlashTime.png");
//
// TCanvas *Canvas2 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
// Canvas2->cd();
// SelectionTrackRange.front()->SetMaximum(0.0035);//SetMaximum(1.1*GetMaximum(SelectionTrackRange));
// SelectionTrackRange.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
// SelectionTrackRange.front()->GetYaxis()->SetTitleOffset(1.3);
// SelectionTrackRange.front()->Draw();
// SelectionTrackRange.at(1)->SetLineColor(2);
// SelectionTrackRange.at(1)->Draw("SAME");
// LegendData->Draw();
// Canvas2->SaveAs("DataSelRange.png");
//
// TCanvas *Canvas3 = new TCanvas("Theta-Angle of Selected Track", "Theta-Angle of Selected Track", 1400, 1000);
// Canvas3->cd();
// SelectionTheta.front()->SetMaximum(/*1.1*GetMaximum(SelectionTheta)*/0.0016);
// SelectionTheta.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
// SelectionTheta.front()->GetYaxis()->SetTitleOffset(1.3);
// SelectionTheta.front()->Draw();
// SelectionTheta.at(1)->SetLineColor(2);
// SelectionTheta.at(1)->Draw("SAME");
// LegendData->Draw();
// Canvas3->SaveAs("DataSelTheta.png");
//
// TCanvas *Canvas4 = new TCanvas("Phi-Angle of Selected Track", "Phi-Angle of Selected Track", 1400, 1000);
// Canvas4->cd();
// SelectionPhi.front()->SetMaximum(/*1.1*GetMaximum(SelectionPhi)*/0.002);
// SelectionPhi.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
// SelectionPhi.front()->GetYaxis()->SetTitleOffset(1.3);
// SelectionPhi.front()->Draw();
// SelectionPhi.at(1)->SetLineColor(2);
// SelectionPhi.at(1)->Draw("SAME");
// LegendData->Draw();
// Canvas4->SaveAs("DataSelPhi.png");
//
// AddFirstTwoHistograms(SelectionTrackRange,-1.);
// AddFirstTwoHistograms(SelectionTheta,-1.);
// AddFirstTwoHistograms(SelectionPhi,-1.);
//
// LegendMC->AddEntry( SelectionTrackRange.at(0), (MCLabel.at(0)).c_str(),"lep" );
// LegendMC->AddEntry( SelectionTrackRange.at(1), (MCLabel.at(1)).c_str(),"f" );
// LegendMC->AddEntry( SelectionNCRange, (MCLabel.at(2)).c_str(),"l" );
//
// TCanvas *Canvas8 = new TCanvas("MC and On-Off Beam Track Length", "MC and On-Off Beam Track Length", 1400, 1000);
// Canvas8->cd();
// SelectionTrackRange.front()->SetMaximum(1.2*GetMaximum(SelectionTrackRange));
// SelectionTrackRange.at(1)->SetFillColorAlpha(46,0.5);
// SelectionTrackRange.at(1)->GetYaxis()->SetTitle("Normalized #frac{dn}{dx}");
// SelectionTrackRange.at(1)->GetYaxis()->SetTitleOffset(1.3);
// SelectionTrackRange.at(1)->DrawNormalized("E2");
// SelectionTrackRange.front()->SetLineColor(1);
// SelectionTrackRange.front()->SetLineWidth(2);
// SelectionTrackRange.front()->DrawNormalized("SAME");
// SelectionNCRange->SetLineColor(2);
// SelectionNCRange->SetLineWidth(2);
// SelectionNCRange->Draw("SAME");
// LegendMC->Draw();
// Canvas8->SaveAs("OffBeam-OnBeamRange.png");
//
// TCanvas *Canvas9 = new TCanvas("MC and On-Off Beam Theta 1", "MC and On-Off Beam Theta 1", 1400, 1000);
// Canvas9->cd();
// SelectionTheta.at(1)->SetMaximum(1.2*GetMaximum(SelectionTheta));
// SelectionTheta.at(1)->SetFillColorAlpha(46,0.5);
// SelectionTheta.at(1)->GetYaxis()->SetTitle("Normalized #frac{dn}{d#theta}");
// SelectionTheta.at(1)->GetYaxis()->SetTitleOffset(1.3);
// SelectionTheta.at(1)->DrawNormalized("E2");
// SelectionTheta.front()->SetLineColor(1);
// SelectionTheta.front()->SetLineWidth(2);
// SelectionTheta.front()->DrawNormalized("SAME");
// SelectionNCTheta->SetLineColor(2);
// SelectionNCTheta->SetLineWidth(2);
// SelectionNCTheta->Draw("SAME");
// LegendMC->Draw();
// Canvas9->SaveAs("OffBeam-OnBeamTheta_1.png");
//
// for(auto& ThetaHistogram : SelectionTheta)
// {
//     ThetaHistogram->Divide(SinTheta,1.);
// }
//
// TCanvas *Canvas10 = new TCanvas("MC and On-Off Beam Theta 2", "MC and On-Off Beam Theta 2", 1400, 1000);
// Canvas10->cd();
// SelectionTheta.at(1)->SetMaximum(1.2*GetMaximum(SelectionTheta));
// SelectionTheta.at(1)->SetFillColorAlpha(46,0.5);
// SelectionTheta.at(1)->GetYaxis()->SetTitle("Normalized #frac{dn}{d#Omega}");
// SelectionTheta.at(1)->GetYaxis()->SetTitleOffset(1.3);
// SelectionTheta.at(1)->DrawNormalized("E2");
// SelectionTheta.front()->SetLineColor(1);
// SelectionTheta.front()->SetLineWidth(2);
// SelectionTheta.front()->DrawNormalized("SAME");
// SelectionNCTheta->SetLineColor(2);
// SelectionNCTheta->SetLineWidth(2);
// SelectionNCTheta->Draw("SAME");
// LegendMC->Draw();
// Canvas10->SaveAs("OffBeam-OnBeamTheta_2.png");
//
// TCanvas *Canvas11 = new TCanvas("MC and On-Off Beam Phi", "MC and On-Off Beam Phi", 1400, 1000);
// Canvas11->cd();
// SelectionPhi.at(1)->SetMaximum(2*GetMaximum(SelectionPhi));
// SelectionPhi.at(1)->SetMinimum(0.0);
// SelectionPhi.at(1)->SetFillColorAlpha(46,0.5);
// SelectionPhi.at(1)->GetYaxis()->SetTitle("Normalized #frac{dn}{d#phi}");
// SelectionPhi.at(1)->GetYaxis()->SetTitleOffset(1.3);
// SelectionPhi.at(1)->DrawNormalized("E2");
// SelectionPhi.front()->SetLineColor(1);
// SelectionPhi.front()->SetLineWidth(2);
// SelectionPhi.front()->DrawNormalized("SAME");
// SelectionNCPhi->SetLineColor(2);
// SelectionNCPhi->SetLineWidth(2);
// SelectionNCPhi->Draw("SAME");
// LegendMC->Draw();
// Canvas11->SaveAs("OffBeam-OnBeamPhi.png");

//
//     TCanvas *Canvas6 = new TCanvas("Flash Track Distance", "Flash Track Distance", 1400, 1000);
//     Canvas6->cd();
//     Canvas6->SetLogy();
//     FlashTrackHistVec.front()->SetMaximum(1.1*GetMaximum(FlashTrackHistVec));
//     FlashTrackHistVec.front()->Draw();
//
//     TCanvas *Canvas7 = new TCanvas("Vertex Track Distance", "Vertex Track Distance", 1400, 1000);
//     Canvas7->cd();
//     Canvas7->SetLogy();
//     VertexTrackHistVec.front()->SetMaximum(1.1*GetMaximum(VertexTrackHistVec));
//     VertexTrackHistVec.front()->Draw();
//
//     TCanvas *Canvas8 = new TCanvas("Track Length", "Track Length", 1400, 1000);
//     Canvas8->cd();
//     Canvas8->SetLogy();
//     TrackLengthHistVec.front()->SetMaximum(1.1*GetMaximum(TrackLengthHistVec));
//     TrackLengthHistVec.front()->Draw();
//
//     TCanvas *Canvas9 = new TCanvas("Track Multiplicity", "Track Multiplicity", 1400, 1000);
//     Canvas9->cd();
//     Canvas9->SetLogy();
//     TrackMultipHistVec.front()->SetMaximum(1.1*GetMaximum(TrackMultipHistVec));
//     TrackMultipHistVec.front()->Draw();
//
//     for(unsigned int no_entries = 1; no_entries < FileNameVec.size(); no_entries++)
//     {
//         Canvas1->cd();
//         XVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//         XVtxHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas2->cd();
//         YVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//         YVtxHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas3->cd();
//         ZVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//         ZVtxHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas6->cd();
//         FlashTrackHistVec.at(no_entries)->SetLineColor(no_entries);
//         FlashTrackHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas7->cd();
//         VertexTrackHistVec.at(no_entries)->SetLineColor(no_entries);
//         VertexTrackHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas8->cd();
//         TrackLengthHistVec.at(no_entries)->SetLineColor(no_entries);
//         TrackLengthHistVec.at(no_entries)->Draw("SAME");
//
//         Canvas9->cd();
//         TrackMultipHistVec.at(no_entries)->SetLineColor(no_entries);
//         TrackMultipHistVec.at(no_entries)->Draw("SAME");
//     }
//
//     // Fill legend
//     for(unsigned int no_entries = 0; no_entries < FileNameVec.size(); no_entries++)
//     {
//         Legend->AddEntry( FlashTrackHistVec.at(no_entries), (FileNameVec.at(no_entries)).c_str(),"l" );
//     }
//
//     Canvas1->cd();
//     Legend->Draw();
//     Canvas1->SaveAs(("XVtxPosition_"+Version+".png").c_str());
//
//     Canvas2->cd();
//     Legend->Draw();
//     Canvas2->SaveAs(("YVtxPosition_"+Version+".png").c_str());
//
//     Canvas3->cd();
//     Legend->Draw();
//     Canvas3->SaveAs(("ZVtxPosition_"+Version+".png").c_str());
//
//     Canvas4->SaveAs(("FlashTime_"+Version+".png").c_str());
//     Canvas5->SaveAs(("PECount_"+Version+".png").c_str());
//
//     Canvas6->cd();
//     Legend->Draw();
//     Canvas6->SaveAs(("FlashTrackDist_"+Version+".png").c_str());
//
//     Canvas7->cd();
//     Legend->Draw();
//     Canvas7->SaveAs(("VtxTrackDist_"+Version+".png").c_str());
//
//     Canvas8->cd();
//     Legend->Draw();
//     Canvas8->SaveAs(("TrackLength_"+Version+".png").c_str());
//
//     Canvas9->cd();
//     Legend->Draw();
//     Canvas9->SaveAs(("TrackMultiplicity_"+Version+".png").c_str());
//
// //   TChain *XVtxHist = new TChain("X Vertex Position");

}

float GetMaximum(const std::vector<TH1F*>& HistVector)
{
    float Maximum = 0;

    for(const auto& Histogram : HistVector)
    {
        float TempMax = Histogram->GetBinContent(Histogram->GetMaximumBin());

        if(TempMax > Maximum)
        {
            Maximum = TempMax;
        }
    }

    return Maximum;
}

void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight)
{
    if(HistVector.size() > 1)
    {
        HistVector.at(0)->Add(HistVector.at(1), Weight);
        delete HistVector.at(1);
        HistVector.erase(HistVector.begin()+1);
    }
    else
    {
        std::cout << "Histograms not added!" << std::endl;
    }
}
