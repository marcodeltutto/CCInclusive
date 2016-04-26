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
#include <THStack.h>
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

    std::vector<std::string> DataLabel;
    std::vector<std::string> MCLabel;
    std::vector<std::string> GenLabel;
    std::vector<std::string> BgrLabel;

    std::vector<float> ScalingFactors;
    ScalingFactors.push_back(1/383519.);
    ScalingFactors.push_back(1/400675.);
    ScalingFactors.push_back(1/550000.);

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

    std::vector<TH1F*> BgrTrackRange;
    std::vector<TH1F*> BgrEnergy;
    std::vector<TH1F*> BgrTheta;
    std::vector<TH1F*> BgrPhi;

    std::vector<TH1F*> BgrXTrackStartEnd;
    std::vector<TH1F*> BgrYTrackStartEnd;
    std::vector<TH1F*> BgrZTrackStartEnd;

    std::vector<TH1F*> BgrXVtxPosition;
    std::vector<TH1F*> BgrYVtxPosition;
    std::vector<TH1F*> BgrZVtxPosition;

    THStack* StackBgrTrackRange = new THStack("Bgr Track Range","Bgr Track Range");
    THStack* StackBgrEnergy = new THStack("Bgr Energy","Bgr Energy");
    THStack* StackBgrTheta = new THStack("Bgr Theta","Bgr Theta");
    THStack* StackBgrPhi = new THStack("Bgr Phi","Bgr Phi");
    THStack* StackBgrXTrackStartEnd = new THStack("Bgr X Start","Bgr X Start");
    THStack* StackBgrYTrackStartEnd = new THStack("Bgr Y Start","Bgr Y Start");
    THStack* StackBgrZTrackStartEnd = new THStack("Bgr Z Start","Bgr Z Start");
    THStack* StackBgrXVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrYVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");
    THStack* StackBgrZVtxPosition = new THStack("Bgr X Vertex","Bgr X Vertex");

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
    GenLabel.push_back("MC Prodgenie BNB Nu Cosmic");
//     GenLabel.push_back("MC Prodgenie BNB Nu");
//     GenLabel.push_back("MC Prodcosmic Corsika in-Time");

    BgrLabel.push_back("Bgr #bar{#nu}_{#mu} Events");
    BgrLabel.push_back("Bgr #nu_{e} Events");
    BgrLabel.push_back("Bgr NC Events");
    BgrLabel.push_back("Bgr Cosmic Events");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");
//     ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_1.root");
//     ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_onbeam_bnb_v05_08_00_2.root");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1.root");
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2.root");
//     ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_1.root");
//     ChainVec.back() -> Add("/media/christoph/200EFBDA63AA160B/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_data_offbeam_bnbext_v05_08_00_2.root");

    ChainVec.push_back(new TChain("anatree"));
    ChainVec.back() -> Add("/lheppc46/data/uBData/anatrees/Hist_Track_pandoraNu_Vertex_pandoraNu_prodgenie_bnb_nu_cosmic_uboone_v05_08_00.root");

    for(const auto& Label : GenLabel)
    {
        SelectionTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",20,0,1000));
        SelectionTrackRange.back()->SetStats(0);
        SelectionTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
        SelectionTrackRange.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");

        SelectionTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta-Angle of Selected Track",20,0,3.142));
        SelectionTheta.back()->SetStats(0);
        SelectionTheta.back()->GetXaxis()->SetTitle("#theta angle [rad]");
        SelectionTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
        SelectionTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        SelectionPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",20,-3.142,3.142));
        SelectionPhi.back()->SetStats(0);
        SelectionPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        SelectionPhi.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
        SelectionPhi.back()->GetYaxis()->SetTitleOffset(1.3);

        SelectionEnergy.push_back(new TH1F(("Energy"+Label).c_str(),"Energy of Selected Track",20,0,3000));
        SelectionEnergy.back()->SetStats(0);
        SelectionEnergy.back()->GetXaxis()->SetTitle("Muon Kinetic Energy [MeV]");
        SelectionEnergy.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dE}");
        SelectionEnergy.back()->GetYaxis()->SetTitleOffset(1.3);

        SelXTrackStartEnd.push_back(new TH1F(("XTrack"+Label).c_str(),"X Track Start & End Positions",20,0,256));
        SelXTrackStartEnd.back()->SetStats(0);
        SelXTrackStartEnd.back()->GetXaxis()->SetTitle("x [cm]");
        SelXTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        SelXTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelYTrackStartEnd.push_back(new TH1F(("YTrack"+Label).c_str(),"Y Track Start & End Positions",20,-233/2,233/2));
        SelYTrackStartEnd.back()->SetStats(0);
        SelYTrackStartEnd.back()->GetXaxis()->SetTitle("y [cm]");
        SelYTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        SelYTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelZTrackStartEnd.push_back(new TH1F(("ZTrack"+Label).c_str(),"Z Track Start & End Positions",20,0,1000));
        SelZTrackStartEnd.back()->SetStats(0);
        SelZTrackStartEnd.back()->GetXaxis()->SetTitle("z [cm]");
        SelZTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        SelZTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        SelXVtxPosition.push_back(new TH1F(("XVertex"+Label).c_str(),"X Vertex Position",20,0,256));
        SelXVtxPosition.back()->SetStats(0);
        SelXVtxPosition.back()->GetXaxis()->SetTitle("x [cm]");
        SelXVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        SelXVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        SelYVtxPosition.push_back(new TH1F(("YVertex"+Label).c_str(),"Y Vertex Position",20,-233/2,233/2));
        SelYVtxPosition.back()->SetStats(0);
        SelYVtxPosition.back()->GetXaxis()->SetTitle("y [cm]");
        SelYVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        SelYVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        SelZVtxPosition.push_back(new TH1F(("ZVertex"+Label).c_str(),"Z Vertex Position",20,0,1000));
        SelZVtxPosition.back()->SetStats(0);
        SelZVtxPosition.back()->GetXaxis()->SetTitle("z [cm]");
        SelZVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        SelZVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);
    }

    for(const auto& Label : BgrLabel)
    {
        BgrTrackRange.push_back(new TH1F(("Track Range"+Label).c_str(),"Track Range of Selected Track",20,0,1000));
        BgrTrackRange.back()->SetStats(0);
        BgrTrackRange.back()->GetXaxis()->SetTitle("Track Range [cm]");
        BgrTrackRange.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");

        BgrTheta.push_back(new TH1F(("#theta-Angle"+Label).c_str(),"#theta-Angle of Selected Track",20,-1,1));
        BgrTheta.back()->SetStats(0);
        BgrTheta.back()->GetXaxis()->SetTitle("#theta angle [rad]");
        BgrTheta.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#theta}");
        BgrTheta.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrPhi.push_back(new TH1F(("#phi-Angle"+Label).c_str(),"#phi-Angle of Selected Track",20,-3.142,3.142));
        BgrPhi.back()->SetStats(0);
        BgrPhi.back()->GetXaxis()->SetTitle("#phi angle [rad]");
        BgrPhi.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{d#phi}");
        BgrPhi.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrEnergy.push_back(new TH1F(("Energy"+Label).c_str(),"Energy of Selected Track",20,0,3000));
        BgrEnergy.back()->SetStats(0);
        BgrEnergy.back()->GetXaxis()->SetTitle("Muon Kinetic Energy [MeV]");
        BgrEnergy.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dE}");
        BgrEnergy.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrXTrackStartEnd.push_back(new TH1F(("XTrack"+Label).c_str(),"X Track Start & End Positions",20,0,256));
        BgrXTrackStartEnd.back()->SetStats(0);
        BgrXTrackStartEnd.back()->GetXaxis()->SetTitle("x [cm]");
        BgrXTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        BgrXTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrYTrackStartEnd.push_back(new TH1F(("YTrack"+Label).c_str(),"Y Track Start & End Positions",20,-233/2,233/2));
        BgrYTrackStartEnd.back()->SetStats(0);
        BgrYTrackStartEnd.back()->GetXaxis()->SetTitle("y [cm]");
        BgrYTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        BgrYTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrZTrackStartEnd.push_back(new TH1F(("ZTrack"+Label).c_str(),"Z Track Start & End Positions",20,0,1000));
        BgrZTrackStartEnd.back()->SetStats(0);
        BgrZTrackStartEnd.back()->GetXaxis()->SetTitle("z [cm]");
        BgrZTrackStartEnd.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        BgrZTrackStartEnd.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrXVtxPosition.push_back(new TH1F(("XVertex"+Label).c_str(),"X Vertex Position",20,0,256));
        BgrXVtxPosition.back()->SetStats(0);
        BgrXVtxPosition.back()->GetXaxis()->SetTitle("x [cm]");
        BgrXVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dx}");
        BgrXVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrYVtxPosition.push_back(new TH1F(("YVertex"+Label).c_str(),"Y Vertex Position",20,-233/2,233/2));
        BgrYVtxPosition.back()->SetStats(0);
        BgrYVtxPosition.back()->GetXaxis()->SetTitle("y [cm]");
        BgrYVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dy}");
        BgrYVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);

        BgrZVtxPosition.push_back(new TH1F(("ZVertex"+Label).c_str(),"Z Vertex Position",20,0,1000));
        BgrZVtxPosition.back()->SetStats(0);
        BgrZVtxPosition.back()->GetXaxis()->SetTitle("z [cm]");
        BgrZVtxPosition.back()->GetYaxis()->SetTitle("Weighted #frac{dn}{dz}");
        BgrZVtxPosition.back()->GetYaxis()->SetTitleOffset(1.3);
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

        unsigned int nubar = 0;
        unsigned int nue = 0;
        unsigned int NCnu = 0;
        unsigned int Cosmic = 0;
        
        for(unsigned int tree_index = 0; tree_index < ChainVec.at(file_no)->GetEntries(); tree_index++)
        {
            if(!(tree_index % 100)) std::cout << "Event\t" << tree_index << "\t of \t" << ChainVec.at(file_no)->GetEntries() << std::endl;

            ChainVec.at(file_no)->GetEntry(tree_index);

            SelectionTrackRange.at(file_no)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
            SelectionTheta.at(file_no)->Fill(cos(TrackTheta[TrkID]));
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
            
            if(file_no == 2 && MCTrkID > -1 && CCNCFlag[0] == 0 && MCTrkOrigin[TrkID][2] == 1)
            {
                if(PDGTruth[MCTrkID] == -13)
                {
                    nubar++;
                    BgrTrackRange.at(0)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    BgrTheta.at(0)->Fill(cos(TrackTheta[TrkID]));
                    BgrPhi.at(0)->Fill(TrackPhi[TrkID]);
                    BgrEnergy.at(0)->Fill(KineticEnergy[TrkID][2]);
                    BgrXTrackStartEnd.at(0)->Fill(XTrackStart[TrkID]);
                    BgrXTrackStartEnd.at(0)->Fill(XTrackEnd[TrkID]);
                    BgrYTrackStartEnd.at(0)->Fill(YTrackStart[TrkID]);
                    BgrYTrackStartEnd.at(0)->Fill(YTrackEnd[TrkID]);
                    BgrZTrackStartEnd.at(0)->Fill(ZTrackStart[TrkID]);
                    BgrZTrackStartEnd.at(0)->Fill(ZTrackEnd[TrkID]);
                    BgrXVtxPosition.at(0)->Fill(XVertexPosition[VtxID]);
                    BgrYVtxPosition.at(0)->Fill(YVertexPosition[VtxID]);
                    BgrZVtxPosition.at(0)->Fill(ZVertexPosition[VtxID]);
                }
                else if(abs(PDGTruth[MCTrkID]) == 11)
                {
                    nue++;
                    BgrTrackRange.at(1)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                    BgrTheta.at(1)->Fill(cos(TrackTheta[TrkID]));
                    BgrPhi.at(1)->Fill(TrackPhi[TrkID]);
                    BgrEnergy.at(1)->Fill(KineticEnergy[TrkID][2]);
                    BgrXTrackStartEnd.at(1)->Fill(XTrackStart[TrkID]);
                    BgrXTrackStartEnd.at(1)->Fill(XTrackEnd[TrkID]);
                    BgrYTrackStartEnd.at(1)->Fill(YTrackStart[TrkID]);
                    BgrYTrackStartEnd.at(1)->Fill(YTrackEnd[TrkID]);
                    BgrZTrackStartEnd.at(1)->Fill(ZTrackStart[TrkID]);
                    BgrZTrackStartEnd.at(1)->Fill(ZTrackEnd[TrkID]);
                    BgrXVtxPosition.at(1)->Fill(XVertexPosition[VtxID]);
                    BgrYVtxPosition.at(1)->Fill(YVertexPosition[VtxID]);
                    BgrZVtxPosition.at(1)->Fill(ZVertexPosition[VtxID]);
                }
            }
            else if(file_no == 2 && CCNCFlag[0] == 1 && MCTrkOrigin[TrkID][2] == 1)
            {
                NCnu++;
                BgrTrackRange.at(2)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                BgrTheta.at(2)->Fill(cos(TrackTheta[TrkID]));
                BgrPhi.at(2)->Fill(TrackPhi[TrkID]);
                BgrEnergy.at(2)->Fill(KineticEnergy[TrkID][2]);
                BgrXTrackStartEnd.at(2)->Fill(XTrackStart[TrkID]);
                BgrXTrackStartEnd.at(2)->Fill(XTrackEnd[TrkID]);
                BgrYTrackStartEnd.at(2)->Fill(YTrackStart[TrkID]);
                BgrYTrackStartEnd.at(2)->Fill(YTrackEnd[TrkID]);
                BgrZTrackStartEnd.at(2)->Fill(ZTrackStart[TrkID]);
                BgrZTrackStartEnd.at(2)->Fill(ZTrackEnd[TrkID]);
                BgrXVtxPosition.at(2)->Fill(XVertexPosition[VtxID]);
                BgrYVtxPosition.at(2)->Fill(YVertexPosition[VtxID]);
                BgrZVtxPosition.at(2)->Fill(ZVertexPosition[VtxID]);
            }
            else if(file_no == 2 && MCTrkOrigin[TrkID][2] == 2)
            {
                Cosmic++;
                BgrTrackRange.at(3)->Fill(CalcLength(XTrackStart[TrkID],YTrackStart[TrkID],ZTrackStart[TrkID],XTrackEnd[TrkID],YTrackEnd[TrkID],ZTrackEnd[TrkID]));
                BgrTheta.at(3)->Fill(cos(TrackTheta[TrkID]));
                BgrPhi.at(3)->Fill(TrackPhi[TrkID]);
                BgrEnergy.at(3)->Fill(KineticEnergy[TrkID][2]);
                BgrXTrackStartEnd.at(3)->Fill(XTrackStart[TrkID]);
                BgrXTrackStartEnd.at(3)->Fill(XTrackEnd[TrkID]);
                BgrYTrackStartEnd.at(3)->Fill(YTrackStart[TrkID]);
                BgrYTrackStartEnd.at(3)->Fill(YTrackEnd[TrkID]);
                BgrZTrackStartEnd.at(3)->Fill(ZTrackStart[TrkID]);
                BgrZTrackStartEnd.at(3)->Fill(ZTrackEnd[TrkID]);
                BgrXVtxPosition.at(3)->Fill(XVertexPosition[VtxID]);
                BgrYVtxPosition.at(3)->Fill(YVertexPosition[VtxID]);
                BgrZVtxPosition.at(3)->Fill(ZVertexPosition[VtxID]);
            }
        }
        std::cout << nubar << " " << nue << " " << NCnu << " " << Cosmic << std::endl;
        
        ChainVec.at(file_no)->ResetBranchAddresses();
    }

    for(unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++)
    {
        BgrTrackRange.at(bgrhist_no)->Scale(1/SelectionTrackRange.at(2)->Integral());
        BgrTheta.at(bgrhist_no)->Scale(1/SelectionTheta.at(2)->Integral());
        BgrPhi.at(bgrhist_no)->Scale(1/SelectionPhi.at(2)->Integral());
        BgrEnergy.at(bgrhist_no)->Scale(1/SelectionEnergy.at(2)->Integral());
        BgrXTrackStartEnd.at(bgrhist_no)->Scale(1/SelXTrackStartEnd.at(2)->Integral());
        BgrYTrackStartEnd.at(bgrhist_no)->Scale(1/SelYTrackStartEnd.at(2)->Integral());
        BgrZTrackStartEnd.at(bgrhist_no)->Scale(1/SelZTrackStartEnd.at(2)->Integral());
        BgrXVtxPosition.at(bgrhist_no)->Scale(1/SelXVtxPosition.at(2)->Integral());
        BgrYVtxPosition.at(bgrhist_no)->Scale(1/SelYVtxPosition.at(2)->Integral());
        BgrZVtxPosition.at(bgrhist_no)->Scale(1/SelZVtxPosition.at(2)->Integral());

        BgrTrackRange.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrTheta.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrPhi.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrEnergy.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrXTrackStartEnd.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrYTrackStartEnd.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrZTrackStartEnd.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrXVtxPosition.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrYVtxPosition.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        BgrZVtxPosition.at(bgrhist_no)->SetLineColor(bgrhist_no+1);
        
        StackBgrTrackRange->Add(BgrTrackRange.at(bgrhist_no));
        StackBgrTheta->Add(BgrTheta.at(bgrhist_no));
        StackBgrPhi->Add(BgrPhi.at(bgrhist_no));
        StackBgrEnergy->Add(BgrEnergy.at(bgrhist_no));
        StackBgrXTrackStartEnd->Add(BgrXTrackStartEnd.at(bgrhist_no));
        StackBgrYTrackStartEnd->Add(BgrYTrackStartEnd.at(bgrhist_no));
        StackBgrZTrackStartEnd->Add(BgrZTrackStartEnd.at(bgrhist_no));
        StackBgrXVtxPosition->Add(BgrXVtxPosition.at(bgrhist_no));
        StackBgrYVtxPosition->Add(BgrYVtxPosition.at(bgrhist_no));
        StackBgrZVtxPosition->Add(BgrZVtxPosition.at(bgrhist_no));
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
    SelectionTrackRange.at(0)->SetMaximum(1.1*GetMaximum(SelectionTrackRange));
    SelectionTrackRange.at(0)->Draw();
    SelectionTrackRange.at(1)->SetLineColor(2);
    SelectionTrackRange.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas1->SaveAs("DataSelRange.png");

    TCanvas *Canvas2 = new TCanvas("Theta-Angle of Selected Track", "Theta-Angle of Selected Track", 1400, 1000);
    Canvas2->cd();
    SelectionTheta.at(0)->SetMaximum(1.1*GetMaximum(SelectionTheta));
    SelectionTheta.at(0)->Draw();
    SelectionTheta.at(1)->SetLineColor(2);
    SelectionTheta.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas2->SaveAs("DataSelTheta.png");

    TCanvas *Canvas3 = new TCanvas("Phi-Angle of Selected Track", "Phi-Angle of Selected Track", 1400, 1000);
    Canvas3->cd();
    SelectionPhi.at(0)->SetMaximum(1.1*GetMaximum(SelectionPhi));
    SelectionPhi.at(0)->Draw();
    SelectionPhi.at(1)->SetLineColor(2);
    SelectionPhi.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas3->SaveAs("DataSelPhi.png");

    TCanvas *Canvas4 = new TCanvas("Energy of Selected Track", "Energy of Selected Track", 1400, 1000);
    Canvas4->cd();
    SelectionEnergy.at(0)->SetMaximum(1.1*GetMaximum(SelectionEnergy));
    SelectionEnergy.at(0)->Draw();
    SelectionEnergy.at(1)->SetLineColor(2);
    SelectionEnergy.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas4->SaveAs("DataSelEnergy.png");

    TCanvas *Canvas5 = new TCanvas("X Start & End Point Selected Track", "X Start & End Point Selected Track", 1400, 1000);
    Canvas5->cd();
    SelXTrackStartEnd.at(0)->SetMaximum(1.3*GetMaximum(SelXTrackStartEnd));
    SelXTrackStartEnd.at(0)->Draw();
    SelXTrackStartEnd.at(1)->SetLineColor(2);
    SelXTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas5->SaveAs("DataSelXTrack.png");

    TCanvas *Canvas6 = new TCanvas("Y Start & End Point Selected Track", "Y Start & End Point Selected Track", 1400, 1000);
    Canvas6->cd();
    SelYTrackStartEnd.at(0)->SetMaximum(1.3*GetMaximum(SelYTrackStartEnd));
    SelYTrackStartEnd.at(0)->Draw();
    SelYTrackStartEnd.at(1)->SetLineColor(2);
    SelYTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas6->SaveAs("DataSelYTrack.png");

    TCanvas *Canvas7 = new TCanvas("Z Start & End Point Selected Track", "Z Start & End Point Selected Track", 1400, 1000);
    Canvas7->cd();
    SelZTrackStartEnd.at(0)->SetMaximum(1.3*GetMaximum(SelZTrackStartEnd));
    SelZTrackStartEnd.at(0)->Draw();
    SelZTrackStartEnd.at(1)->SetLineColor(2);
    SelZTrackStartEnd.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas7->SaveAs("DataSelZTrack.png");

    TCanvas *Canvas8 = new TCanvas("X Vertex Postion", "X Vertex Postion", 1400, 1000);
    Canvas8->cd();
    SelXVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelXVtxPosition));
    SelXVtxPosition.at(0)->Draw();
    SelXVtxPosition.at(1)->SetLineColor(2);
    SelXVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas8->SaveAs("DataSelXVertex.png");

    TCanvas *Canvas9 = new TCanvas("Y Vertex Postion", "Y Vertex Postion", 1400, 1000);
    Canvas9->cd();
    SelYVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelYVtxPosition));
    SelYVtxPosition.at(0)->Draw();
    SelYVtxPosition.at(1)->SetLineColor(2);
    SelYVtxPosition.at(1)->Draw("SAME");
    LegendData->Draw();
    Canvas9->SaveAs("DataSelYVertex.png");

    TCanvas *Canvas10 = new TCanvas("Z Vertex Postion", "Z Vertex Postion", 1400, 1000);
    Canvas10->cd();
    SelZVtxPosition.at(0)->SetMaximum(1.3*GetMaximum(SelZVtxPosition));
    SelZVtxPosition.at(0)->Draw();
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

    LegendMC->AddEntry( SelectionTrackRange.at(0), (MCLabel.at(0)).c_str(),"lep" );
    LegendMC->AddEntry( SelectionTrackRange.at(1), (MCLabel.at(1)).c_str(),"f" );
    for(unsigned int bgrhist_no = 0; bgrhist_no < BgrLabel.size(); bgrhist_no++)
    {
        LegendMC->AddEntry( BgrTrackRange.at(bgrhist_no), (BgrLabel.at(bgrhist_no)).c_str(),"l" );
    }

    TCanvas *Canvas11 = new TCanvas("OnBeam Minus OffBeam Track Range", "OnBeam Minus OffBeam Track Range", 1400, 1000);
    Canvas11->cd();
    SelectionTrackRange.at(1)->SetMaximum(1.5*SelectionTrackRange.at(1)->GetBinContent(SelectionTrackRange.at(1)->GetMaximumBin()));
    SelectionTrackRange.at(1)->SetMinimum(0.0);
    SelectionTrackRange.at(1)->SetFillColorAlpha(46,0.5);
    SelectionTrackRange.at(1)->DrawNormalized("E2");
    SelectionTrackRange.at(0)->SetLineWidth(2);
    SelectionTrackRange.at(0)->SetLineColor(1);
    SelectionTrackRange.at(0)->DrawNormalized("SAME");
    StackBgrTrackRange->Draw("SAME");
    LegendMC->Draw();
    Canvas11->SaveAs("On-OffBeamSelRange.png");

    TCanvas *Canvas12 = new TCanvas("OnBeam Minus OffBeam Theta-Angle", "OnBeam Minus OffBeam Theta-Angle", 1400, 1000);
    Canvas12->cd();
    SelectionTheta.at(1)->SetMaximum(1.5*SelectionTheta.at(1)->GetBinContent(SelectionTheta.at(1)->GetMaximumBin()));
    SelectionTheta.at(1)->SetMinimum(0.0);
    SelectionTheta.at(1)->SetFillColorAlpha(46,0.5);
    SelectionTheta.at(1)->DrawNormalized("E2");
    SelectionTheta.at(0)->SetLineWidth(2);
    SelectionTheta.at(0)->SetLineColor(1);
    SelectionTheta.at(0)->DrawNormalized("SAME");
    StackBgrTheta->Draw("SAME");
    LegendMC->Draw();
    Canvas12->SaveAs("On-OffBeamSelTheta.png");

    TCanvas *Canvas13 = new TCanvas("OnBeam Minus OffBeam Phi-Angle", "OnBeam Minus OffBeam Phi-Angle", 1400, 1000);
    Canvas13->cd();
    SelectionPhi.at(1)->SetMaximum(2*SelectionPhi.at(1)->GetBinContent(SelectionPhi.at(1)->GetMaximumBin()));
    SelectionPhi.at(1)->SetMinimum(0.0);
    SelectionPhi.at(1)->SetFillColorAlpha(46,0.5);
    SelectionPhi.at(1)->DrawNormalized("E2");
    SelectionPhi.at(0)->SetLineWidth(2);
    SelectionPhi.at(0)->SetLineColor(1);
    SelectionPhi.at(0)->DrawNormalized("SAME");
    StackBgrPhi->Draw("SAME");
    LegendMC->Draw();
    Canvas13->SaveAs("On-OffBeamSelPhi.png");

    TCanvas *Canvas14 = new TCanvas("Energy", "Energy", 1400, 1000);
    Canvas14->cd();
    SelectionEnergy.at(1)->SetMaximum(1.5*SelectionEnergy.at(1)->GetBinContent(SelectionEnergy.at(1)->GetMaximumBin()));
    SelectionEnergy.at(1)->SetMinimum(0.0);
    SelectionEnergy.at(1)->SetFillColorAlpha(46,0.5);
    SelectionEnergy.at(1)->DrawNormalized("E2");
    SelectionEnergy.at(0)->SetLineWidth(2);
    SelectionEnergy.at(0)->SetLineColor(1);
    SelectionEnergy.at(0)->DrawNormalized("SAME");
    StackBgrEnergy->Draw("SAME");
    LegendMC->Draw();
    Canvas14->SaveAs("On-OffBeamSelEnergy.png");

    TCanvas *Canvas15 = new TCanvas("OnBeam Minus OffBeam X Start & End Point ", "OnBeam Minus OffBeam X Start & End Point ", 1400, 1000);
    Canvas15->cd();
    SelXTrackStartEnd.at(1)->SetMaximum(1.5*SelXTrackStartEnd.at(1)->GetBinContent(SelXTrackStartEnd.at(1)->GetMaximumBin()));
    SelXTrackStartEnd.at(1)->SetMinimum(0.0);
    SelXTrackStartEnd.at(1)->SetFillColorAlpha(46,0.5);
    SelXTrackStartEnd.at(1)->DrawNormalized("E2");
    SelXTrackStartEnd.at(0)->SetLineWidth(2);
    SelXTrackStartEnd.at(0)->SetLineColor(1);
    SelXTrackStartEnd.at(0)->DrawNormalized("SAME");
    StackBgrXTrackStartEnd->Draw("SAME");
    LegendMC->Draw();
    Canvas15->SaveAs("On-OffBeamSelXTrack.png");

    TCanvas *Canvas16 = new TCanvas("OnBeam Minus OffBeam Y Start & End Point ", "OnBeam Minus OffBeam Y Start & End Point ", 1400, 1000);
    Canvas16->cd();
    SelYTrackStartEnd.at(1)->SetMaximum(1.5*SelYTrackStartEnd.at(1)->GetBinContent(SelYTrackStartEnd.at(1)->GetMaximumBin()));
    SelYTrackStartEnd.at(1)->SetMinimum(0.0);
    SelYTrackStartEnd.at(1)->SetFillColorAlpha(46,0.5);
    SelYTrackStartEnd.at(1)->DrawNormalized("E2");
    SelYTrackStartEnd.at(0)->SetLineWidth(2);
    SelYTrackStartEnd.at(0)->SetLineColor(1);
    SelYTrackStartEnd.at(0)->DrawNormalized("SAME");
    StackBgrYTrackStartEnd->Draw("SAME");
    LegendMC->Draw();
    Canvas16->SaveAs("On-OffBeamSelYTrack.png");

    TCanvas *Canvas17 = new TCanvas("OnBeam Minus OffBeam Z Start & End Point ", "OnBeam Minus OffBeam Z Start & End Point ", 1400, 1000);
    Canvas17->cd();
    SelZTrackStartEnd.at(1)->SetMaximum(1.5*SelZTrackStartEnd.at(1)->GetBinContent(SelZTrackStartEnd.at(1)->GetMaximumBin()));
    SelZTrackStartEnd.at(1)->SetMinimum(0.0);
    SelZTrackStartEnd.at(1)->SetFillColorAlpha(46,0.5);
    SelZTrackStartEnd.at(1)->DrawNormalized("E2");
    SelZTrackStartEnd.at(0)->SetLineWidth(2);
    SelZTrackStartEnd.at(0)->SetLineColor(1);
    SelZTrackStartEnd.at(0)->DrawNormalized("SAME");
    StackBgrZTrackStartEnd->Draw("SAME");
    LegendMC->Draw();
    Canvas17->SaveAs("On-OffBeamSelZTrack.png");

    TCanvas *Canvas18 = new TCanvas("OnBeam Minus OffBeam X Vertex Postion", "OnBeam Minus OffBeam X Vertex Postion", 1400, 1000);
    Canvas18->cd();
    SelXVtxPosition.at(1)->SetMaximum(1.5*SelXVtxPosition.at(1)->GetBinContent(SelXVtxPosition.at(1)->GetMaximumBin()));
    SelXVtxPosition.at(1)->SetMinimum(0.0);
    SelXVtxPosition.at(1)->SetFillColorAlpha(46,0.5);
    SelXVtxPosition.at(1)->DrawNormalized("E2");
    SelXVtxPosition.at(0)->SetLineWidth(2);
    SelXVtxPosition.at(0)->SetLineColor(1);
    SelXVtxPosition.at(0)->DrawNormalized("SAME");
    StackBgrXVtxPosition->Draw("SAME");
    LegendMC->Draw();
    Canvas18->SaveAs("On-OffBeamSelXVertex.png");

    TCanvas *Canvas19 = new TCanvas("OnBeam Minus OffBeam Y Vertex Postion", "OnBeam Minus OffBeam Y Vertex Postion", 1400, 1000);
    Canvas19->cd();
    SelYVtxPosition.at(1)->SetMaximum(1.5*SelYVtxPosition.at(1)->GetBinContent(SelYVtxPosition.at(1)->GetMaximumBin()));
    SelYVtxPosition.at(1)->SetMinimum(0.0);
    SelYVtxPosition.at(1)->SetFillColorAlpha(46,0.5);
    SelYVtxPosition.at(1)->DrawNormalized("E2");
    SelYVtxPosition.at(0)->SetLineWidth(2);
    SelYVtxPosition.at(0)->SetLineColor(1);
    SelYVtxPosition.at(0)->DrawNormalized("SAME");
    StackBgrYVtxPosition->Draw("SAME");
    LegendMC->Draw();
    Canvas19->SaveAs("On-OffBeamSelYVertex.png");

    TCanvas *Canvas20 = new TCanvas("OnBeam Minus OffBeam Z Vertex Postion", "OnBeam Minus OffBeam Z Vertex Postion", 1400, 1000);
    Canvas20->cd();
    SelZVtxPosition.at(1)->SetMaximum(1.5*SelZVtxPosition.at(1)->GetBinContent(SelZVtxPosition.at(1)->GetMaximumBin()));
    SelZVtxPosition.at(1)->SetMinimum(0.0);
    SelZVtxPosition.at(1)->SetFillColorAlpha(46,0.5);
    SelZVtxPosition.at(1)->DrawNormalized("E2");
    SelZVtxPosition.at(0)->SetLineWidth(2);
    SelZVtxPosition.at(0)->SetLineColor(1);
    SelZVtxPosition.at(0)->DrawNormalized("SAME");
    StackBgrZVtxPosition->Draw("SAME");
    LegendMC->Draw();
    Canvas20->SaveAs("On-OffBeamSelZVertex.png");
    
    TCanvas *Canvas21 = new TCanvas("Bgr Z Vertex Postion", "Bgr Z Vertex Postion", 1400, 1000);
    StackBgrZVtxPosition->Draw();
    
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
