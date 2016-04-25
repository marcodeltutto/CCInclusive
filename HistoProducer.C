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
float CalcLength(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);

void HistoProducer()
{
    TGaxis::SetMaxDigits(4);
    std::vector<TChain*> ChainVec;

    std::vector<std::string> FileNameVec;
    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;

    std::vector<float> ScalingFactors;
    ScalingFactors.push_back(1/383519.);
    ScalingFactors.push_back(1/400675.);

    std::vector<TH1F*> FlashTime;
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionEnergy;
    std::vector<TH1F*> SelectionTheta;
    std::vector<TH1F*> SelectionPhi;
    
    std::vector<TH1F*> SelXTrackStartEnd;
    std::vector<TH1F*> SelYTrackStartEnd;
    std::vector<TH1F*> SelZTrackStartEnd;
    
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
//     GenLabel.push_back("MC Prodgenie BNB Nu");
//     GenLabel.push_back("MC Prodgenie BNB Nu Cosmic");
//     GenLabel.push_back("MC Prodcosmic Corsika in-Time");

    ChainVec.push_back(new TChain("anatree"));
//     ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
//     ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");
    ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
    ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");

    ChainVec.push_back(new TChain("anatree"));
//     ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1.root");
//     ChainVec.back() -> Add("rootfiles/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2.root");
    ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1.root");
    ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2.root");

    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",20,0,1000));
        SelectionTrackRange.back()->SetStats(0);
        SelectionTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
        SelectionTrackRange.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");

        SelectionTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta-Angle of Selected Track",20,0,3.142));
        SelectionTheta.back()->SetStats(0);
        SelectionTheta.back()->GetXaxis()->SetTitle("#theta angle [rad]");
        SelectionTheta.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",20,-3.142,3.142));
        SelectionPhi.back()->SetStats(0);
        SelectionPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        SelectionPhi.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");

        SelectionEnergy.push_back(new TH1F(("Energy"+Label).c_str(),"Energy of Selected Track",20,0,3000));
        SelectionEnergy.back()->SetStats(0);
        SelectionEnergy.back()->GetXaxis()->SetTitle("Muon Kinetic Energy [MeV]");
        SelectionEnergy.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");
        
        SelXTrackStartEnd.push_back(new TH1F(("XTrack"+Label).c_str(),"X Track Start & End Positions",20,0,256));
        SelXTrackStartEnd.back()->SetStats(0);
        SelXTrackStartEnd.back()->GetXaxis()->SetTitle("x [cm]");
        SelXTrackStartEnd.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");
        
        SelYTrackStartEnd.push_back(new TH1F(("YTrack"+Label).c_str(),"Y Track Start & End Positions",20,-233/2,233/2));
        SelYTrackStartEnd.back()->SetStats(0);
        SelYTrackStartEnd.back()->GetXaxis()->SetTitle("y [cm]");
        SelYTrackStartEnd.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");
        
        SelZTrackStartEnd.push_back(new TH1F(("ZTrack"+Label).c_str(),"Z Track Start & End Positions",20,0,1000));
        SelZTrackStartEnd.back()->SetStats(0);
        SelZTrackStartEnd.back()->GetXaxis()->SetTitle("z [cm]");
        SelZTrackStartEnd.back()->GetYaxis()->SetTitle("Number of Tracks [ ]");
        
        SelXVtxPosition.push_back(new TH1F(("XVertex"+Label).c_str(),"X Vertex Position",20,0,256));
        SelXVtxPosition.back()->SetStats(0);
        SelXVtxPosition.back()->GetXaxis()->SetTitle("x [cm]");
        SelXVtxPosition.back()->GetYaxis()->SetTitle("Number of Vertices [ ]");
        
        SelYVtxPosition.push_back(new TH1F(("YVertex"+Label).c_str(),"Y Vertex Position",20,-233/2,233/2));
        SelYVtxPosition.back()->SetStats(0);
        SelYVtxPosition.back()->GetXaxis()->SetTitle("y [cm]");
        SelYVtxPosition.back()->GetYaxis()->SetTitle("Number of Vertices [ ]");
        
        SelZVtxPosition.push_back(new TH1F(("ZVertex"+Label).c_str(),"Z Vertex Position",20,0,1000));
        SelZVtxPosition.back()->SetStats(0);
        SelZVtxPosition.back()->GetXaxis()->SetTitle("z [cm]");
        SelZVtxPosition.back()->GetYaxis()->SetTitle("Number of Vertices [ ]");
    }

    int TrkID;
    int VtxID;
    
    int MCTrkID;
    int CCNCFlag[10];
    int PDGTruth[5000];
    int MCTrkOrigin[5000][3];

    float TrackTheta[5000];
    float TrackPhi[5000];

    float XTrackStart[5000];
    float YTrackStart[5000];
    float ZTrackStart[5000];

    float XTrackEnd[5000];
    float YTrackEnd[5000];
    float ZTrackEnd[5000];
    
    float XVertexPosition[500];
    float YVertexPosition[500];
    float ZVertexPosition[500];

    float KineticEnergy[5000][3];

    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
        ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);
        
        ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
        ChainVec.at(file_no) -> SetBranchAddress("ccnc_truth", CCNCFlag);
        ChainVec.at(file_no) -> SetBranchAddress("pdg", PDGTruth);
        ChainVec.at(file_no) -> SetBranchAddress("trkorigin_pandoraNu", MCTrkOrigin);

        ChainVec.at(file_no) -> SetBranchAddress("trkke_pandoraNu", KineticEnergy);
        ChainVec.at(file_no) -> SetBranchAddress("trktheta_pandoraNu", TrackTheta);
        ChainVec.at(file_no) -> SetBranchAddress("trkphi_pandoraNu",TrackPhi);

        ChainVec.at(file_no) -> SetBranchAddress("trkstartx_pandoraNu",XTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress("trkstarty_pandoraNu",YTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress("trkstartz_pandoraNu",ZTrackStart);

        ChainVec.at(file_no) -> SetBranchAddress("trkendx_pandoraNu",XTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress("trkendy_pandoraNu",YTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress("trkendz_pandoraNu",ZTrackEnd);
        
        ChainVec.at(file_no) -> SetBranchAddress("vtxx_pandoraNu", XVertexPosition);
        ChainVec.at(file_no) -> SetBranchAddress("vtxy_pandoraNu", YVertexPosition);
        ChainVec.at(file_no) -> SetBranchAddress("vtxz_pandoraNu", ZVertexPosition);

        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no)->GetEntries(); tree_index++)
        {
            if(!(tree_index % 100)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no)->GetEntries() << std::endl;

            ChainVec.at(file_no)->GetEntry(tree_index);
            
            SelectionTrackRange.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
            SelectionTheta.at(file_no)->Fill(TrackTheta[TrkID]);
            SelectionPhi.at(file_no)->Fill(TrackPhi[TrkID]);
            SelectionEnergy.at(file_no)->Fill(KineticEnergy[TrkID][2]);
            
            SelXTrackStartEnd.at(file_no)->Fill(XTrackStart[TrkID]);
            SelXTrackStartEnd.at(file_no)->Fill(XTrackEnd[TrkID]);
            SelYTrackStartEnd.at(file_no)->Fill(YTrackStart[TrkID]);
            SelYTrackStartEnd.at(file_no)->Fill(YTrackEnd[TrkID]);
            SelZTrackStartEnd.at(file_no)->Fill(ZTrackStart[TrkID]);
            SelZTrackStartEnd.at(file_no)->Fill(ZTrackEnd[TrkID]);
            
            SelXVtxPosition.at(file_no)->Fill(XVertexPosition[VtxID]);
            SelYVtxPosition.at(file_no)->Fill(YVertexPosition[VtxID]);
            SelZVtxPosition.at(file_no)->Fill(ZVertexPosition[VtxID]);
        }
        ChainVec.at(file_no)->ResetBranchAddresses();
    }


    for(unsigned int file_no = 0; file_no < ScalingFactors.size(); file_no++)
    {
        SelectionTrackRange.at(file_no)->Sumw2();
        SelectionTheta.at(file_no)->Sumw2();
        SelectionPhi.at(file_no)->Sumw2();
        SelectionEnergy.at(file_no)->Sumw2();
        SelXTrackStartEnd.at(file_no)->Sumw2();
        SelYTrackStartEnd.at(file_no)->Sumw2();
        SelZTrackStartEnd.at(file_no)->Sumw2();
        SelXVtxPosition.at(file_no)->Sumw2();
        SelYVtxPosition.at(file_no)->Sumw2();
        SelZVtxPosition.at(file_no)->Sumw2();

        SelectionTrackRange.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelectionTheta.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelectionPhi.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelectionEnergy.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelXTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelYTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelZTrackStartEnd.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelXVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelYVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
        SelZVtxPosition.at(file_no)->Scale(ScalingFactors.at(file_no));
    }

    for(unsigned int hist_no = 0; hist_no < DataLabel.size(); hist_no++)
    {
        LegendData->AddEntry( SelectionTrackRange.at(hist_no), (DataLabel.at(hist_no)).c_str(),"l" );
    }

    TCanvas *Canvas1 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
    Canvas1->cd();
    SelectionTrackRange.front()->SetMaximum(1.1*GetMaximum(SelectionTrackRange));
    SelectionTrackRange.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
    SelectionTrackRange.front()->GetYaxis()->SetTitleOffset(1.3);
    SelectionTrackRange.front()->Draw();
    SelectionTrackRange.at(1)->SetLineColor(2);
    SelectionTrackRange.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas1->SaveAs("DataSelRange.png");

    TCanvas *Canvas2 = new TCanvas("Theta-Angle of Selected Track", "Theta-Angle of Selected Track", 1400, 1000);
    Canvas2->cd();
    SelectionTheta.front()->SetMaximum(1.1*GetMaximum(SelectionTheta));
    SelectionTheta.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
    SelectionTheta.front()->GetYaxis()->SetTitleOffset(1.3);
    SelectionTheta.front()->Draw();
    SelectionTheta.at(1)->SetLineColor(2);
    SelectionTheta.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas2->SaveAs("DataSelTheta.png");

    TCanvas *Canvas3 = new TCanvas("Phi-Angle of Selected Track", "Phi-Angle of Selected Track", 1400, 1000);
    Canvas3->cd();
    SelectionPhi.front()->SetMaximum(1.1*GetMaximum(SelectionPhi));
    SelectionPhi.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
    SelectionPhi.front()->GetYaxis()->SetTitleOffset(1.3);
    SelectionPhi.front()->Draw();
    SelectionPhi.at(1)->SetLineColor(2);
    SelectionPhi.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas3->SaveAs("DataSelPhi.png");

    TCanvas *Canvas4 = new TCanvas("Energy of Selected Track", "Energy of Selected Track", 1400, 1000);
    Canvas4->cd();
    SelectionEnergy.front()->SetMaximum(1.1*GetMaximum(SelectionEnergy));
    SelectionEnergy.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dE}");
    SelectionEnergy.front()->GetYaxis()->SetTitleOffset(1.3);
    SelectionEnergy.front()->Draw();
    SelectionEnergy.at(1)->SetLineColor(2);
    SelectionEnergy.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas4->SaveAs("DataSelEnergy.png");
    
    TCanvas *Canvas5 = new TCanvas("X Start & End Point Selected Track", "X Start & End Point Selected Track", 1400, 1000);
    Canvas5->cd();
    SelXTrackStartEnd.front()->SetMaximum(1.3*GetMaximum(SelXTrackStartEnd));
    SelXTrackStartEnd.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
    SelXTrackStartEnd.front()->GetYaxis()->SetTitleOffset(1.3);
    SelXTrackStartEnd.front()->Draw();
    SelXTrackStartEnd.at(1)->SetLineColor(2);
    SelXTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas5->SaveAs("DataSelXTrack.png");
    
    TCanvas *Canvas6 = new TCanvas("Y Start & End Point Selected Track", "Y Start & End Point Selected Track", 1400, 1000);
    Canvas6->cd();
    SelYTrackStartEnd.front()->SetMaximum(1.3*GetMaximum(SelYTrackStartEnd));
    SelYTrackStartEnd.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
    SelYTrackStartEnd.front()->GetYaxis()->SetTitleOffset(1.3);
    SelYTrackStartEnd.front()->Draw();
    SelYTrackStartEnd.at(1)->SetLineColor(2);
    SelYTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas6->SaveAs("DataSelYTrack.png");
    
    TCanvas *Canvas7 = new TCanvas("Z Start & End Point Selected Track", "Z Start & End Point Selected Track", 1400, 1000);
    Canvas7->cd();
    SelZTrackStartEnd.front()->SetMaximum(1.3*GetMaximum(SelZTrackStartEnd));
    SelZTrackStartEnd.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
    SelZTrackStartEnd.front()->GetYaxis()->SetTitleOffset(1.3);
    SelZTrackStartEnd.front()->Draw();
    SelZTrackStartEnd.at(1)->SetLineColor(2);
    SelZTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas7->SaveAs("DataSelZTrack.png");
    
    TCanvas *Canvas8 = new TCanvas("X Vertex Postion", "X Vertex Postion", 1400, 1000);
    Canvas8->cd();
    SelXVtxPosition.front()->SetMaximum(1.3*GetMaximum(SelXVtxPosition));
    SelXVtxPosition.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
    SelXVtxPosition.front()->GetYaxis()->SetTitleOffset(1.3);
    SelXVtxPosition.front()->Draw();
    SelXVtxPosition.at(1)->SetLineColor(2);
    SelXVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas8->SaveAs("DataSelXVertex.png");
    
    TCanvas *Canvas9 = new TCanvas("Y Vertex Postion", "Y Vertex Postion", 1400, 1000);
    Canvas9->cd();
    SelYVtxPosition.front()->SetMaximum(1.3*GetMaximum(SelYVtxPosition));
    SelYVtxPosition.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
    SelYVtxPosition.front()->GetYaxis()->SetTitleOffset(1.3);
    SelYVtxPosition.front()->Draw();
    SelYVtxPosition.at(1)->SetLineColor(2);
    SelYVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas9->SaveAs("DataSelYVertex.png");
    
    TCanvas *Canvas10 = new TCanvas("Z Vertex Postion", "Z Vertex Postion", 1400, 1000);
    Canvas10->cd();
    SelZVtxPosition.front()->SetMaximum(1.3*GetMaximum(SelZVtxPosition));
    SelZVtxPosition.front()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
    SelZVtxPosition.front()->GetYaxis()->SetTitleOffset(1.3);
    SelZVtxPosition.front()->Draw();
    SelZVtxPosition.at(1)->SetLineColor(2);
    SelZVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas10->SaveAs("DataSelZVertex.png");

    AddFirstTwoHistograms(SelectionTrackRange,-1.);
    AddFirstTwoHistograms(SelectionTheta,-1.);
    AddFirstTwoHistograms(SelectionPhi,-1.);
    AddFirstTwoHistograms(SelectionEnergy,-1);
    AddFirstTwoHistograms(SelXTrackStartEnd,-1);
    AddFirstTwoHistograms(SelYTrackStartEnd,-1);
    AddFirstTwoHistograms(SelZTrackStartEnd,-1);
    AddFirstTwoHistograms(SelXVtxPosition,-1);
    AddFirstTwoHistograms(SelYVtxPosition,-1);
    AddFirstTwoHistograms(SelZVtxPosition,-1);
    
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

float CalcLength(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2)
{
    return sqrt(pow(x_1-x_2, 2) + pow(y_1-y_2, 2) + pow(z_1-z_2, 2));
}
