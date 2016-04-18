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

float GetMaximum(std::vector<TH1F*> HistVector);

void Draw_numuCC()
{
    //Define cut variables
    double FVx = 256.35; //cm. Active volume in x
    double FVy = 233; //cm. Active volume in y
    double FVz = 1036.8; //cm. Active volume in z
    double borderx = 10.; //cm. FV cut in x
    double bordery = 20.; //cm. FV cut in y
    double borderz = 10.; //cm. FV cut in y
    double flashwidth = 80; //cm. Distance flash-track
    double distcut = 5; //cm. Distance track start/end to vertex
    double lengthcut = 75; //cm. Length of longest track
    double beammin = 3.55/*-0.36*/; //us. Beam window start
    double beammax = 5.15/*-0.36*/; //us. Beam window end
    double PEthresh = 50; //PE
    double MCTrackToMCVtxDist = 1; //cm. distance between mc track start and mc vertex
    double TrackToMCDist = 5; //cm. Distance track start/end to mcvertex

    std::vector<std::string> DataSampleNameVec;
    std::vector<std::string> TrackingNameVec;
    std::vector<std::string> VertexingNameVec;

    std::vector<TFile*> FileVec;

    std::vector<TH1F*> XVtxHistVec;
    std::vector<TH1F*> YVtxHistVec;
    std::vector<TH1F*> ZVtxHistVec;
    std::vector<TH1F*> FlashTrackHistVec;
    std::vector<TH1F*> VertexTrackHistVec;
    std::vector<TH1F*> TrackLengthHistVec;
    std::vector<TH1F*> TrackMultipHistVec;
    std::vector<TH1F*> FlashVertexDistVec;
    std::vector<TH1F*> XTrackStartEndVec;
    std::vector<TH1F*> YTrackStartEndVec;
    std::vector<TH1F*> ZTrackStartEndVec;

    TH1F* FlashTimeHist;
    TH1F* FlashPEHist;

    std::string Version = "v05_08_00";

    DataSampleNameVec.push_back("prodgenie_bnb_nu");

    VertexingNameVec.push_back("nuvtx");
    VertexingNameVec.push_back("pandoraCosmic");
    VertexingNameVec.push_back("pandoraNu");
    VertexingNameVec.push_back("pmtrack");

    TrackingNameVec.push_back("pandoraNuKHit");
    TrackingNameVec.push_back("pandoraNu");
    TrackingNameVec.push_back("pandoraCosmic");
    TrackingNameVec.push_back("pmtrack");
    TrackingNameVec.push_back("pandoraNuPMA");
    TrackingNameVec.push_back("trackkalmanhit");

    // Loop over data sample names
    for(const auto& DataName : DataSampleNameVec)
    {

        // Vertex and track independent variables
        FileVec.push_back( new TFile(("Hist_Track_"+TrackingNameVec.front()+"_Vertex_"+VertexingNameVec.front()+"_"+DataName+"_"+Version+".root").c_str(),"READ") );
        FileVec.back()->cd();

        FlashTimeHist = (TH1F*)FileVec.back()->Get("Flash Time");
        FlashPEHist = (TH1F*)FileVec.back()->Get("Flash Intensity in PE");

        TCanvas *Canvas1 = new TCanvas("Flash Time Distribution", "Flash Time Distribution", 1400, 1000);
        Canvas1->cd();
        TLine FlashTimeMinCut(beammin, 0, beammin, 1200);
        TLine FlashTimeMaxCut(beammax, 0, beammax, 1200);
        FlashTimeMinCut.SetLineColor(kRed);
        FlashTimeMaxCut.SetLineColor(kRed);
        FlashTimeHist->Draw();
        FlashTimeMinCut.Draw("same");
        FlashTimeMaxCut.Draw("same");

        TCanvas *Canvas2 = new TCanvas("PE Count of Flash", "PE Count of Flash", 1400, 1000);
        Canvas2->cd();
        Canvas2->SetLogy();
        TLine FlashPEMinCut(PEthresh, 0, PEthresh, 20000);
        FlashPEMinCut.SetLineColor(kRed);
        FlashPEHist->Draw();
        FlashPEMinCut.Draw("same");

        Canvas1->SaveAs(("FlashTime_Vertex_"+DataName+"_"+Version+".png").c_str());
        Canvas2->SaveAs(("PECount_Vertex_"+DataName+"_"+Version+".png").c_str());

        // Clean up
        delete Canvas1;
        delete Canvas2;
        FileVec.back()->Close();
        FileVec.pop_back();
        // End vertex and track independent variables
        //---------------------------------------------------------------------------------------------------------------------------------------------

        // Only vertex dependent variables
        for(const auto& VertexName : VertexingNameVec)
        {
            FileVec.push_back( new TFile(("Hist_Track_"+TrackingNameVec.front()+"_Vertex_"+VertexName+"_"+DataName+"_"+Version+".root").c_str(),"READ") );
            FileVec.back()->cd();

            XVtxHistVec.push_back((TH1F*)FileVec.back()->Get("X Vertex Distribution"));
            YVtxHistVec.push_back((TH1F*)FileVec.back()->Get("Y Vertex Distribution"));
            ZVtxHistVec.push_back((TH1F*)FileVec.back()->Get("Z Vertex Distribution"));
            FlashVertexDistVec.push_back((TH1F*)FileVec.back()->Get("Flash To Vertex Distance"));
        }

        TCanvas *Canvas3 = new TCanvas("X Vertex Position", "X Vertex Position", 1400, 1000);
        Canvas3->cd();
        XVtxHistVec.front()->SetMaximum(1.1*GetMaximum(XVtxHistVec));
        XVtxHistVec.front()->Draw();

        TCanvas *Canvas4 = new TCanvas("Y Vertex Position", "Y Vertex Position", 1400, 1000);
        Canvas4->cd();
        YVtxHistVec.front()->SetMaximum(1.1*GetMaximum(YVtxHistVec));
        YVtxHistVec.front()->Draw();

        TCanvas *Canvas5 = new TCanvas("Z Vertex Position", "Z Vertex Position", 1400, 1000);
        Canvas5->cd();
        ZVtxHistVec.front()->SetMaximum(1.1*GetMaximum(ZVtxHistVec));
        ZVtxHistVec.front()->Draw();

        TCanvas *Canvas6 = new TCanvas("Vertex to Flash Distance", "Vertex to Flash Distance", 1400, 1000);
        Canvas6->cd();
        FlashVertexDistVec.front()->SetMaximum(1.1*GetMaximum(FlashVertexDistVec));
        FlashVertexDistVec.front()->Draw();

        for(unsigned int no_entries = 1; no_entries < VertexingNameVec.size(); no_entries++)
        {
            Canvas3->cd();
            XVtxHistVec.at(no_entries)->SetLineColor(no_entries);
            XVtxHistVec.at(no_entries)->Draw("same");

            Canvas4->cd();
            YVtxHistVec.at(no_entries)->SetLineColor(no_entries);
            YVtxHistVec.at(no_entries)->Draw("same");

            Canvas5->cd();
            ZVtxHistVec.at(no_entries)->SetLineColor(no_entries);
            ZVtxHistVec.at(no_entries)->Draw("same");

            Canvas6->cd();
            FlashVertexDistVec.at(no_entries)->SetLineColor(no_entries);
            FlashVertexDistVec.at(no_entries)->Draw("same");
        }

        TLegend* Legend = new TLegend(0.7,0.7,0.9,0.9);
        Legend->SetHeader("Vertex Reco Method");

        for(unsigned int no_entries = 0; no_entries < VertexingNameVec.size(); no_entries++)
        {
            Legend->AddEntry( XVtxHistVec.at(no_entries), (VertexingNameVec.at(no_entries)).c_str(),"l" );
        }

        Canvas3->cd();
        TLine XVertexMinCut(borderx, 0, borderx, GetMaximum(XVtxHistVec));
        TLine XVertexMaxCut(FVx-borderx, 0, FVx-borderx, GetMaximum(XVtxHistVec));
        XVertexMinCut.SetLineColor(kRed);
        XVertexMaxCut.SetLineColor(kRed);
        XVertexMinCut.Draw("same");
        XVertexMaxCut.Draw("same");
        Legend->Draw();
        Canvas3->SaveAs(("XVtxPosition_Vertex_"+DataName+"_"+Version+".png").c_str());

        Canvas4->cd();
        TLine YVertexMinCut(-FVy/2+bordery, 0, -FVy/2+bordery, GetMaximum(YVtxHistVec));
        TLine YVertexMaxCut(FVy/2-bordery, 0, FVy/2-bordery, GetMaximum(YVtxHistVec));
        YVertexMinCut.SetLineColor(kRed);
        YVertexMaxCut.SetLineColor(kRed);
        YVertexMinCut.Draw("same");
        YVertexMaxCut.Draw("same");
        Legend->Draw();
        Canvas4->SaveAs(("YVtxPosition_Vertex_"+DataName+"_"+Version+".png").c_str());

        Canvas5->cd();
        Legend->Draw();
        TLine ZVertexMinCut(borderz, 0, borderz, GetMaximum(ZVtxHistVec));
        TLine ZVertexMaxCut(FVz-borderz, 0, FVz-borderz, GetMaximum(ZVtxHistVec));
        ZVertexMinCut.SetLineColor(kRed);
        ZVertexMaxCut.SetLineColor(kRed);
        ZVertexMinCut.Draw("same");
        ZVertexMaxCut.Draw("same");
        Canvas5->SaveAs(("ZVtxPosition_Vertex_"+DataName+"_"+Version+".png").c_str());

        Canvas6->cd();
        Legend->Draw();
        Canvas6->SaveAs(("FlashVtxDist_Vertex_"+DataName+"_"+Version+".png").c_str());

        delete Canvas3;
        delete Canvas4;
        delete Canvas5;
        delete Canvas6;
        delete Legend;

        for(const auto& VertexName : VertexingNameVec)
        {
            delete XVtxHistVec.back();
            delete YVtxHistVec.back();
            delete ZVtxHistVec.back();
            delete FlashVertexDistVec.back();
            FileVec.back()->Close();

            XVtxHistVec.pop_back();
            YVtxHistVec.pop_back();
            ZVtxHistVec.pop_back();
            FlashVertexDistVec.pop_back();
            FileVec.pop_back();
        }
        // End only vertex dependent variables
        //---------------------------------------------------------------------------------------------------------------------------------------------


        // Begin of vertex and track dependent variables
        // Loop over vertex reco method names
        for(const auto& VertexName : VertexingNameVec)
        {
            unsigned int VertexNumber = 0;

            // Loop over track reco method names
            for(const auto& TrackName : TrackingNameVec)
            {
                FileVec.push_back( new TFile(("Hist_Track_"+TrackName+"_Vertex_"+VertexName+"_"+DataName+"_"+Version+".root").c_str(),"READ") );
                FileVec.back()->cd();

                FlashTrackHistVec.push_back((TH1F*)FileVec.back()->Get("Distance from Flash to Track"));
                VertexTrackHistVec.push_back((TH1F*)FileVec.back()->Get("Distance from Vertex to Track"));
                TrackLengthHistVec.push_back((TH1F*)FileVec.back()->Get("Track Length of Longest Tracks"));
                TrackMultipHistVec.push_back((TH1F*)FileVec.back()->Get("Track Multiplicity"));
                XTrackStartEndVec.push_back((TH1F*)FileVec.back()->Get("Start & End X-Position of longest Track"));
                YTrackStartEndVec.push_back((TH1F*)FileVec.back()->Get("Start & End Y-Position of longest Track"));
                ZTrackStartEndVec.push_back((TH1F*)FileVec.back()->Get("Start & End Z-Position of longest Track"));

            } // loop over tracking method names


            TCanvas *Canvas7 = new TCanvas("Flash Track Distance", "Flash Track Distance", 1400, 1000);
            Canvas7->cd();
            Canvas7->SetLogy();
            FlashTrackHistVec.front()->SetMaximum(1.1*GetMaximum(FlashTrackHistVec));
            FlashTrackHistVec.front()->Draw();

            TCanvas *Canvas8 = new TCanvas("Vertex Track Distance", "Vertex Track Distance", 1400, 1000);
            Canvas8->cd();
            Canvas8->SetLogy();
            VertexTrackHistVec.front()->SetMaximum(1.1*GetMaximum(VertexTrackHistVec));
            VertexTrackHistVec.front()->Draw();

            TCanvas *Canvas9 = new TCanvas("Track Length", "Track Length", 1400, 1000);
            Canvas9->cd();
            Canvas9->SetLogy();
            TrackLengthHistVec.front()->SetMaximum(1.1*GetMaximum(TrackLengthHistVec));
            TrackLengthHistVec.front()->Draw();

            TCanvas *Canvas10 = new TCanvas("Track Multiplicity", "Track Multiplicity", 1400, 1000);
            Canvas10->cd();
            Canvas10->SetLogy();
            TrackMultipHistVec.front()->SetMaximum(1.1*GetMaximum(TrackMultipHistVec));
            TrackMultipHistVec.front()->Draw();

            TCanvas *Canvas11 = new TCanvas("X Track Start & End Position", "X Track Start & End Position", 1400, 1000);
            Canvas11->cd();
            XTrackStartEndVec.front()->SetMaximum(1.1*GetMaximum(XTrackStartEndVec));
            XTrackStartEndVec.front()->Draw();

            TCanvas *Canvas12 = new TCanvas("Y Track Start & End Position", "Y Track Start & End Position", 1400, 1000);
            Canvas12->cd();
            YTrackStartEndVec.front()->SetMaximum(1.1*GetMaximum(YTrackStartEndVec));
            YTrackStartEndVec.front()->Draw();

            TCanvas *Canvas13 = new TCanvas("Z Track Start & End Position", "Z Track Start & End Position", 1400, 1000);
            Canvas13->cd();
            ZTrackStartEndVec.front()->SetMaximum(1.1*GetMaximum(ZTrackStartEndVec));
            ZTrackStartEndVec.front()->Draw();

            for(unsigned int no_entries = 1; no_entries < TrackingNameVec.size(); no_entries++)
            {
                Canvas7->cd();
                FlashTrackHistVec.at(no_entries)->SetLineColor(no_entries);
                FlashTrackHistVec.at(no_entries)->Draw("same");

                Canvas8->cd();
                VertexTrackHistVec.at(no_entries)->SetLineColor(no_entries);
                VertexTrackHistVec.at(no_entries)->Draw("same");

                Canvas9->cd();
                TrackLengthHistVec.at(no_entries)->SetLineColor(no_entries);
                TrackLengthHistVec.at(no_entries)->Draw("same");

                Canvas10->cd();
                TrackMultipHistVec.at(no_entries)->SetLineColor(no_entries);
                TrackMultipHistVec.at(no_entries)->Draw("same");

                Canvas11->cd();
                XTrackStartEndVec.at(no_entries)->SetLineColor(no_entries);
                XTrackStartEndVec.at(no_entries)->Draw("same");

                Canvas12->cd();
                YTrackStartEndVec.at(no_entries)->SetLineColor(no_entries);
                YTrackStartEndVec.at(no_entries)->Draw("same");

                Canvas13->cd();
                ZTrackStartEndVec.at(no_entries)->SetLineColor(no_entries);
                ZTrackStartEndVec.at(no_entries)->Draw("same");
            }

            TLegend* Legend = new TLegend(0.7,0.7,0.9,0.9);
            Legend->SetHeader("Track Reco Method");

            for(unsigned int no_entries = 0; no_entries < TrackingNameVec.size(); no_entries++)
            {
                Legend->AddEntry( FlashTrackHistVec.at(no_entries), (TrackingNameVec.at(no_entries)).c_str(),"l" );
            }

            Canvas7->cd();
            TLine FlashTrackMinCut(flashwidth, 0, flashwidth, GetMaximum(FlashTrackHistVec));
            FlashTrackMinCut.SetLineColor(kRed);
            FlashTrackMinCut.Draw("same");
            Legend->Draw();
            Canvas7->SaveAs(("FlashTrackDist_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas8->cd();
            TLine VertexTrackMinCut(distcut, 0, distcut, GetMaximum(VertexTrackHistVec));
            VertexTrackMinCut.SetLineColor(kRed);
            VertexTrackMinCut.Draw("same");
            Legend->Draw();
            Canvas8->SaveAs(("VtxTrackDist_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas9->cd();
            Legend->Draw();
            TLine TrackLengthCut(lengthcut, 0, lengthcut, GetMaximum(TrackLengthHistVec));
            TrackLengthCut.SetLineColor(kRed);
            TrackLengthCut.Draw("same");
            Canvas9->SaveAs(("TrackLength_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas10->cd();
            Legend->Draw();
            Canvas10->SaveAs(("TrackMultiplicity_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas11->cd();
            TLine XTrackEndMinCut(borderx, 0, borderx, GetMaximum(XTrackStartEndVec));
            TLine XTrackEndMaxCut(FVx-borderx, 0, FVx-borderx, GetMaximum(XTrackStartEndVec));
            XTrackEndMinCut.SetLineColor(kRed);
            XTrackEndMaxCut.SetLineColor(kRed);
            XTrackEndMinCut.Draw("same");
            XTrackEndMaxCut.Draw("same");
            Legend->Draw();
            Canvas11->SaveAs(("XTrackEndPos_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas12->cd();
            TLine YTrackEndMinCut(-FVy/2+bordery, 0, -FVy/2+bordery, GetMaximum(YTrackStartEndVec));
            TLine YTrackEndMaxCut(FVy/2-bordery, 0, FVy/2-bordery, GetMaximum(YTrackStartEndVec));
            YTrackEndMinCut.SetLineColor(kRed);
            YTrackEndMaxCut.SetLineColor(kRed);
            YTrackEndMinCut.Draw("same");
            YTrackEndMaxCut.Draw("same");
            Legend->Draw();
            Canvas12->SaveAs(("YTrackEndPos_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            Canvas13->cd();
            TLine ZTrackEndMinCut(borderz, 0, borderz, GetMaximum(ZTrackStartEndVec));
            TLine ZTrackEndMaxCut(FVz-borderz, 0, FVz-borderz, GetMaximum(ZTrackStartEndVec));
            ZTrackEndMinCut.SetLineColor(kRed);
            ZTrackEndMaxCut.SetLineColor(kRed);
            ZTrackEndMinCut.Draw("same");
            ZTrackEndMaxCut.Draw("same");
            Legend->Draw();
            Canvas13->SaveAs(("ZTrackEndPos_Vertex_"+VertexName+"_"+DataName+"_"+Version+".png").c_str());

            delete Canvas7;
            delete Canvas8;
            delete Canvas9;
            delete Canvas10;
            delete Canvas11;
            delete Canvas12;
            delete Canvas13;
            delete Legend;

            // Garbage control
            for(const auto& TrackName : TrackingNameVec)
            {
                delete FlashTrackHistVec.back();
                delete VertexTrackHistVec.back();
                delete TrackLengthHistVec.back();
                delete TrackMultipHistVec.back();
                delete XTrackStartEndVec.back();
                delete YTrackStartEndVec.back();
                delete ZTrackStartEndVec.back();
                FileVec.back()->Close();

                FlashTrackHistVec.pop_back();
                VertexTrackHistVec.pop_back();
                TrackLengthHistVec.pop_back();
                TrackMultipHistVec.pop_back();
                XTrackStartEndVec.pop_back();
                YTrackStartEndVec.pop_back();
                ZTrackStartEndVec.pop_back();
                FileVec.pop_back();
            }
            // End of vertex and track dependent variables
            //-------------------------------------------------------------------------------------------------------------------------
        } // loop over vetex reco method names
    } // loop over data sample names

//   TChain *XVtxHist = new TChain("X Vertex Position");

}

float GetMaximum(std::vector< TH1F* > HistVector)
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
