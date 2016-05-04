#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>


float GetMaximum(const std::vector<TH1F*>& HistVector);
void AddFirstTwoHistograms(std::vector<TH1F*>& HistVector, float Weight);
void AddFirstTwoHistograms2D(std::vector<TH2F*>& HistVector, float Weight);
float CalcLength(const float& x_1, const float& y_1, const float& z_1, const float& x_2, const float& y_2, const float& z_2);
double FlashTrackDist(double flash, double start, double end);

void PEChecker()
{
    TGaxis::SetMaxDigits(4);
    std::vector<TChain*> ChainVec;

    std::vector<float> PECutValueVec;
    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;
    std::vector<std::string> BgrLabel;

    std::vector<float> ScalingFactors;
    ScalingFactors.push_back(1/383519.);
    ScalingFactors.push_back(1/400675.);
    ScalingFactors.push_back(1/550000.);



    std::vector<std::vector<TH1F*>> SelectionTrackRange;

//     std::string TrackProdName="pandoraNuKHit";
//     std::string TrackProdName = "pandoraCosmic";
    std::string TrackProdName="pandoraNu";
//     std::string TrackProdName="pmtrack";
//     std::string TrackProdName="pandoraNuPMA";
//     std::string TrackProdName="trackkalmanhit";

//     std::string  VertexProdName="nuvtx";
//     std::string VertexProdName="pandoraCosmic";
    std::string VertexProdName = "pandoraNu";
//     std::string VertexProdName = "pmtrack";

    TF1* SinTheta = new TF1("const","sin(x)",0,3.142);

    PECutValueVec = {50, 70, 100, 250};

    THStack* StackBgrTrackRange = new THStack("Bgr Track Range","Bgr Track Range");

    TLegend* LegendData = new TLegend(0.7,0.8,0.9,0.9);
    LegendData->SetHeader("PE Cut Value");
    
    for(const auto& CutValue : PECutValueVec)
    {
        DataLabel.push_back("Number of PE > "+std::to_string(CutValue));
    }

    TLegend* LegendMC = new TLegend(0.6,0.72,0.9,0.9);
//     LegendMC->SetHeader("Data Type");

    MCLabel.push_back("On-Beam Minus Off-Beam Sample");
    MCLabel.push_back("Selection on MC BNB+Cosmic with Stat. Error");

    TLegend* FlashLabel = new TLegend(0.7,0.7,0.9,0.9);
//     FlashLabel->SetHeader("Generator Type");



    GenLabel.push_back("Data On-Beam BNB");
    GenLabel.push_back("Data Off-Beam BNBEXT");
    GenLabel.push_back("MC Prodgenie BNB Nu Cosmic");
//     GenLabel.push_back("MC Prodgenie BNB Nu");
//     GenLabel.push_back("MC Prodcosmic Corsika in-Time");

    BgrLabel.push_back("Bgr #bar{#nu}_{#mu} Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr #nu_{e} Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr NC Events MC BNB+Cosmic");
    BgrLabel.push_back("Bgr Cosmic Events MC BNB+Cosmic");

    std::vector<unsigned int> ColorMap = {28,42,30,38};

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_1.root").c_str());
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_2.root").c_str());
//     ChainVec.back() -> Add(("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_1.root").c_str());
//     ChainVec.back() -> Add(("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_onbeam_bnb_v05_08_00_2.root").c_str());

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_1.root").c_str());
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_2.root").c_str());
//     ChainVec.back() -> Add(("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_1.root").c_str());
//     ChainVec.back() -> Add(("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_data_offbeam_bnbext_v05_08_00_2.root").c_str());

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add(("/lheppc46/data/uBData/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root").c_str());
//     ChainVec.back() -> Add(("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_"+ TrackProdName +"_Vertex_"+ VertexProdName +"_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root").c_str());

    for(const auto& CutValue : PECutValueVec)
    {
        SelectionTrackRange.push_back(std::vector<TH1F*>());
        for(const auto& Label : GenLabel)
        {
            SelectionTrackRange.back().push_back(new TH1F(("Track Range"+Label+std::to_string(CutValue)).c_str(),"Track Range of Selected Track",20,0,1000));
            SelectionTrackRange.back().back()->SetStats(0);
            SelectionTrackRange.back().back()->GetXaxis()->SetTitle("Track Range [cm]");
            SelectionTrackRange.back().back()->GetYaxis()->SetTitle("Proportional to Number of Events [ ]");
        }
    }

//     unsigned int BgrCount = 0;
//     for(const auto& Label : BgrLabel)
//     {
//         BgrTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",20,0,1000));
//         BgrTrackRange.back()->SetStats(0);
//         BgrTrackRange.back()->SetFillColor(ColorMap.at(BgrCount));
//         BgrTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
//         BgrTrackRange.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
// 
//         BgrCount++;
//     }

    int TrkID;
    int VtxID;

    float FlashPE[5000];
    int NumberOfFlashes;
    float FlashTime[5000];
    float ZFlashCenter[5000];

    int MCTrkID;
    int CCNCFlag[10];
    int PDGTruth[5000];
    short TrkOrigin[5000][3];

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

    double beammin;
    double beammax;

    for(unsigned int file_no = 0; file_no < ChainVec.size(); file_no++)
    {
        if(!file_no)
        {
            beammin = 3.55 - 0.36;
            beammax = 5.15 - 0.36;
        }
        else
        {
            beammin = 3.55;
            beammax = 5.15;
        }

        ChainVec.at(file_no) -> SetBranchAddress("TrackCand", &TrkID);
        ChainVec.at(file_no) -> SetBranchAddress("VertexCand", &VtxID);

        ChainVec.at(file_no) -> SetBranchAddress("flash_pe", FlashPE);
        ChainVec.at(file_no) -> SetBranchAddress("no_flashes", &NumberOfFlashes);
        ChainVec.at(file_no) -> SetBranchAddress("flash_time", FlashTime);
        ChainVec.at(file_no) -> SetBranchAddress("flash_zcenter", ZFlashCenter);

        ChainVec.at(file_no) -> SetBranchAddress("MCTrackCand", &MCTrkID);
        ChainVec.at(file_no) -> SetBranchAddress("ccnc_truth", CCNCFlag);
        ChainVec.at(file_no) -> SetBranchAddress("pdg", PDGTruth);
        ChainVec.at(file_no) -> SetBranchAddress(("trkorigin_"+TrackProdName).c_str(), TrkOrigin);

        ChainVec.at(file_no) -> SetBranchAddress(("trkke_"+TrackProdName).c_str(), KineticEnergy);
        ChainVec.at(file_no) -> SetBranchAddress(("trktheta_"+TrackProdName).c_str(), TrackTheta);
        ChainVec.at(file_no) -> SetBranchAddress(("trkphi_"+TrackProdName).c_str(),TrackPhi);

        ChainVec.at(file_no) -> SetBranchAddress(("trkstartx_"+TrackProdName).c_str(),XTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstarty_"+TrackProdName).c_str(),YTrackStart);
        ChainVec.at(file_no) -> SetBranchAddress(("trkstartz_"+TrackProdName).c_str(),ZTrackStart);

        ChainVec.at(file_no) -> SetBranchAddress(("trkendx_"+TrackProdName).c_str(),XTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendy_"+TrackProdName).c_str(),YTrackEnd);
        ChainVec.at(file_no) -> SetBranchAddress(("trkendz_"+TrackProdName).c_str(),ZTrackEnd);

        if(VertexProdName != "nuvtx")
        {
            ChainVec.at(file_no) -> SetBranchAddress(("vtxx_"+VertexProdName).c_str(), XVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress(("vtxy_"+VertexProdName).c_str(), YVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress(("vtxz_"+VertexProdName).c_str(), ZVertexPosition);
        }
        else
        {
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxx", XVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxy", YVertexPosition);
            ChainVec.at(file_no) -> SetBranchAddress("nuvtxz", ZVertexPosition);
        }

        unsigned int nubar = 0;
        unsigned int nue = 0;
        unsigned int NCnu = 0;
        unsigned int Cosmic = 0;

        unsigned int Signal = 0;

        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no)->GetEntries(); tree_index++)
        {
            if(!(tree_index % 100)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no)->GetEntries() << std::endl;

            ChainVec.at(file_no)->GetEntry(tree_index);

            float FlashMax = 0.0;
            float ZFlashCenterMax = 0.0;

            for(int flash_no = 0; flash_no < NumberOfFlashes; flash_no++)
            {
                if( FlashTime[flash_no] > beammin && FlashTime[flash_no] < beammax && FlashPE[flash_no] > FlashMax)
                {
                    FlashMax = FlashPE[flash_no];
                    ZFlashCenterMax = ZFlashCenter[flash_no];
                }
            }

            for(unsigned int cut_index = 0; cut_index < PECutValueVec.size(); cut_index++)
            {

                if(FlashMax > PECutValueVec.at(cut_index))
                {
                    Signal++;

                    SelectionTrackRange.at(cut_index).at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));

//                     if(file_no == 2 && MCTrkID > -1 && CCNCFlag[0] == 0 && TrkOrigin[TrkID][2] == 1)
//                     {
//                         if(PDGTruth[MCTrkID] == -13)
//                         {
//                             nubar++;
//                             BgrTrackRange.at(0)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
//                         }
//                         else if(abs(PDGTruth[MCTrkID]) == 11)
//                         {
//                             nue++;
//                             BgrTrackRange.at(1)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
//                         }
//                     }
//                     else if(file_no == 2 && CCNCFlag[0] == 1 && TrkOrigin[TrkID][2] == 1)
//                     {
//                         NCnu++;
//                         BgrTrackRange.at(2)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
//                     }
//                     else if(file_no == 2 && TrkOrigin[TrkID][2] != 1)
//                     {
//                         Cosmic++;
//                         BgrTrackRange.at(3)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
//                     }
                }
            }

        }
        ChainVec.at(file_no)->ResetBranchAddresses();
    }

    for(unsigned int cut_index = 0; cut_index < PECutValueVec.size(); cut_index++)
    {

//         for(unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++)
//         {
//             BgrTrackRange.at(bgrhist_no)->Scale(1/SelectionTrackRange.at(2)->Integral());
//             StackBgrTrackRange->Add(BgrTrackRange.at(bgrhist_no));
//         }

        for(unsigned int file_no = 0; file_no < ScalingFactors.size(); file_no++)
        {
            SelectionTrackRange.at(cut_index).at(file_no)->Sumw2();
            SelectionTrackRange.at(cut_index).at(file_no)->Scale(ScalingFactors.at(file_no));
        }
    }
    
    std::cout << DataLabel.size() << " " << SelectionTrackRange.size() << " " << PECutValueVec.size() << std::endl;
    
    for(auto& SelectionHistoSet : SelectionTrackRange)
    {
        AddFirstTwoHistograms(SelectionHistoSet,-1.);
    }
    
    for(unsigned int hist_no = 0; hist_no < DataLabel.size(); hist_no++)
    {
        LegendData->AddEntry( SelectionTrackRange.at(hist_no).at(1), (DataLabel.at(hist_no)).c_str(),"l" );
    }
    
    TCanvas *Canvas1 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
    Canvas1->cd();
//     SelectionTrackRange.at(0).at(1)->SetMaximum(1.1*GetMaximum(SelectionTrackRange.at(0)));
    SelectionTrackRange.at(0).at(1)->SetMinimum(0.0);
    SelectionTrackRange.at(0).at(1)->Draw();
    SelectionTrackRange.at(0).at(1)->SetLineColor(1);
    for(unsigned int cut_index = 1; cut_index < PECutValueVec.size(); cut_index++)
    {
        SelectionTrackRange.at(cut_index).at(1)->SetLineColor(cut_index+1);
        SelectionTrackRange.at(cut_index).at(1)->Draw("SAME");
    }
    LegendData->Draw();
    Canvas1->SaveAs("MCRangeByPE.png");
    
    
    for(unsigned int cut_index = PECutValueVec.size()-1; cut_index > 0; cut_index--)
    {
        SelectionTrackRange.at(cut_index).at(1)->Divide(SelectionTrackRange.at(0).at(1));
    }
    
    
    TCanvas *Canvas2 = new TCanvas("Track Range of Selected Track", "Track Range of Selected Track", 1400, 1000);
    Canvas2->cd();
//     SelectionTrackRange.at(1).at(1)->SetMaximum(1.1*GetMaximum(SelectionTrackRange.at(0)));
//     SelectionTrackRange.at(1).at(1)->SetMinimum(0.0);
    SelectionTrackRange.at(1).at(1)->GetYaxis()->SetTitle("Ratio of PE Cuts [ ]");
    SelectionTrackRange.at(1).at(1)->Draw();
    SelectionTrackRange.at(1).at(1)->SetLineColor(2);
    for(unsigned int cut_index = 2; cut_index < PECutValueVec.size(); cut_index++)
    {
        SelectionTrackRange.at(cut_index).at(1)->SetLineColor(cut_index+1);
        SelectionTrackRange.at(cut_index).at(1)->Draw("SAME");
    }
//     LegendData->Draw();
    Canvas2->SaveAs("MCRangeByPERatio.png");

//     LegendMC->AddEntry( SelectionTrackRange.at(0), (MCLabel.at(0)).c_str(),"lep" );
//     LegendMC->AddEntry( SelectionTrackRange.at(1), (MCLabel.at(1)).c_str(),"f" );
//     for(unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++)
//     {
//         LegendMC->AddEntry( BgrTrackRange.at(bgrhist_no), (BgrLabel.at(bgrhist_no)).c_str(),"f" );
//     }
// 
//     TCanvas *Canvas11 = new TCanvas("OnBeam Minus OffBeam Track Range", "OnBeam Minus OffBeam Track Range", 1400, 1000);
//     Canvas11->cd();
//     SelectionTrackRange.at(1)->SetMaximum(1.4*SelectionTrackRange.at(1)->GetBinContent(SelectionTrackRange.at(1)->GetMaximumBin()));
//     SelectionTrackRange.at(1)->SetMinimum(0.0);
//     SelectionTrackRange.at(1)->SetFillColorAlpha(46,0.5);
//     SelectionTrackRange.at(1)->DrawNormalized("E2");
//     StackBgrTrackRange->Draw("SAME");
//     SelectionTrackRange.at(0)->SetLineWidth(2);
//     SelectionTrackRange.at(0)->SetLineColor(1);
//     SelectionTrackRange.at(0)->DrawNormalized("SAME");
//     LegendMC->Draw();
//     Canvas11->SaveAs("On-OffBeamSelRange.png");
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

void AddFirstTwoHistograms2D(std::vector<TH2F*>& HistVector, float Weight)
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

double FlashTrackDist(double flash, double start, double end)
{
    if(end >= start)
    {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
    else
    {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
}

