#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>

double beammin = 3.55-0.36; //us. Beam window start
double beammax = 5.15-0.36; //us. Beam window end

float GetMaximum(const std::vector<TH1F*>& HistVector);
void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight);

void HistoSubtractor()
{
    std::vector<std::string> FileNameVec;
    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;

    std::vector<TH1F*> FlashTime;
    std::vector<TH1F*> SelectionTrackRange;
    std::vector<TH1F*> SelectionTheta;

    TLegend* LegendData = new TLegend(0.7,0.7,0.9,0.9);
    LegendData->SetHeader("Data Type");

    DataLabel.push_back("On-Beam");
    DataLabel.push_back("Off-Beam");

    TLegend* LegendMC = new TLegend(0.7,0.7,0.9,0.9);
    LegendMC->SetHeader("Data Type");

    MCLabel.push_back("On- Minus Off-Beam");
    MCLabel.push_back("MC BNB Only");

    FileNameVec.push_back("Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
    FileNameVec.push_back("Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");
    FileNameVec.push_back("Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00.root");
    FileNameVec.push_back("Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_v05_08_00.root");

    TLine FlashTimeMinCut(beammin, 0, beammin, 0.1);
    TLine FlashTimeMaxCut(beammax, 0, beammax, 0.1);
    FlashTimeMinCut.SetLineColor(1);
    FlashTimeMaxCut.SetLineColor(1);


    for(const auto& FileName : FileNameVec)
    {
        TFile* File = new TFile( (FileName).c_str(),"READ" ) ;
        File->cd();

        FlashTime.push_back((TH1F*)File->Get("Flash Time"));
        SelectionTrackRange.push_back((TH1F*)File->Get("Track Range of Selected Track"));
        SelectionTheta.push_back((TH1F*)File->Get("#theta-Angle of Selected Track"));
    }

    AddFirstTwoHistograms(FlashTime,1.);
    AddFirstTwoHistograms(SelectionTrackRange,1.);
    AddFirstTwoHistograms(SelectionTheta,1.);

    for(unsigned int hist_no = 0; hist_no < FlashTime.size(); hist_no++)
    {
        FlashTime.at(hist_no)->Rebin(3);
        SelectionTrackRange.at(hist_no)->Rebin(3);
        SelectionTheta.at(hist_no)->Rebin(3);
        if(hist_no < DataLabel.size())
        {
            LegendData->AddEntry( FlashTime.at(hist_no), (DataLabel.at(hist_no)).c_str(),"l" );
        }
    }

    FlashTime.front()->Scale(1./567157.);
    FlashTime.at(1)->Scale(1./130967.);
    FlashTime.at(2)->Scale(1./10880./7);
    SelectionTrackRange.front()->Scale(1./567157.);
    SelectionTrackRange.at(1)->Scale(1./130967);
    SelectionTrackRange.at(2)->Scale(1./10880./7);
    SelectionTheta.front()->Scale(1./567157.);
    SelectionTheta.at(1)->Scale(1./130967.);
    SelectionTheta.at(2)->Scale(1./10880./7);

    TCanvas *Canvas1 = new TCanvas("Flash Time", "Flash Time", 1400, 1000);
    Canvas1->cd();
    FlashTime.front()->SetMaximum(1.1*GetMaximum(FlashTime));
    FlashTime.front()->Draw();
    FlashTime.at(1)->SetLineColor(2);
    FlashTime.at(1)->Draw("SAME");
    FlashTimeMinCut.Draw("same");
    FlashTimeMaxCut.Draw("same");
    LegendData->Draw();
    Canvas1->SaveAs("DataFlashTime.png");

    TCanvas *Canvas2 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
    Canvas2->cd();
    SelectionTrackRange.front()->SetMaximum(1.1*GetMaximum(SelectionTrackRange));
    SelectionTrackRange.front()->Draw();
    SelectionTrackRange.at(1)->SetLineColor(2);
    SelectionTrackRange.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas2->SaveAs("DataSelRange.png");

    TCanvas *Canvas3 = new TCanvas("Theta-Angle of Selected Track", "Theta-Angle of Selected Track", 1400, 1000);
    Canvas3->cd();
    SelectionTheta.front()->SetMaximum(1.1*GetMaximum(SelectionTheta));
    SelectionTheta.front()->Draw();
    SelectionTheta.at(1)->SetLineColor(2);
    SelectionTheta.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas3->SaveAs("DataSelTheta.png");

    AddFirstTwoHistograms(SelectionTrackRange,-1.);
    AddFirstTwoHistograms(SelectionTheta,-1.);

    for(unsigned int hist_no = 0; hist_no < SelectionTrackRange.size(); hist_no++)
    {
        LegendMC->AddEntry( SelectionTrackRange.at(hist_no), (MCLabel.at(hist_no)).c_str(),"l" );
    }

    TCanvas *Canvas4 = new TCanvas("MC and On-Off Beam Track Length", "MC and On-Off Beam Track Length", 1400, 1000);
    Canvas4->cd();
    SelectionTrackRange.front()->SetMaximum(1.1*GetMaximum(SelectionTrackRange));
    SelectionTrackRange.front()->Draw();
    SelectionTrackRange.at(1)->SetLineColor(2);
    SelectionTrackRange.at(1)->Draw("SAME");
    LegendMC->Draw();
    Canvas4->SaveAs("OffBeam-OnBeamRange.png");

    TCanvas *Canvas5 = new TCanvas("MC and On-Off Beam Theta", "MC and On-Off Beam Theta", 1400, 1000);
    Canvas5->cd();
    SelectionTheta.front()->SetMaximum(1.1*GetMaximum(SelectionTheta));
    SelectionTheta.front()->Draw();
    SelectionTheta.at(1)->SetLineColor(2);
    SelectionTheta.at(1)->Draw("SAME");
    LegendMC->Draw();
    Canvas5->SaveAs("OffBeam-OnBeamTheta.png");

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
//         XVtxHistVec.at(no_entries)->Draw("same");
//
//         Canvas2->cd();
//         YVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//         YVtxHistVec.at(no_entries)->Draw("same");
//
//         Canvas3->cd();
//         ZVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//         ZVtxHistVec.at(no_entries)->Draw("same");
//
//         Canvas6->cd();
//         FlashTrackHistVec.at(no_entries)->SetLineColor(no_entries);
//         FlashTrackHistVec.at(no_entries)->Draw("same");
//
//         Canvas7->cd();
//         VertexTrackHistVec.at(no_entries)->SetLineColor(no_entries);
//         VertexTrackHistVec.at(no_entries)->Draw("same");
//
//         Canvas8->cd();
//         TrackLengthHistVec.at(no_entries)->SetLineColor(no_entries);
//         TrackLengthHistVec.at(no_entries)->Draw("same");
//
//         Canvas9->cd();
//         TrackMultipHistVec.at(no_entries)->SetLineColor(no_entries);
//         TrackMultipHistVec.at(no_entries)->Draw("same");
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
