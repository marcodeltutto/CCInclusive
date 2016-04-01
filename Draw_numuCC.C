#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>

float GetMaximum(std::vector<TH1F*> HistVector);

void Draw_numuCC()
{
  
  std::vector<std::string> FileNameVec;
//   std::vector<TFile*> FileVec;
  
  std::vector<TH1F*> XVtxHistVec;
  std::vector<TH1F*> YVtxHistVec;
  std::vector<TH1F*> ZVtxHistVec;
  std::vector<TH1F*> FlashTrackHistVec;
  std::vector<TH1F*> VertexTrackHistVec;
  std::vector<TH1F*> TrackLengthHistVec;
  std::vector<TH1F*> TrackMultipHistVec;
  
  std::vector<float> XVtxMax;
  std::vector<float> XVtYMax;
  std::vector<float> XVtZMax;
  std::vector<float> FlashTrackMax;
  std::vector<float> VertexTrackMax;
  std::vector<float> TrackLengthMax;
  
  TLegend* Legend = new TLegend(0.7,0.7,0.9,0.9);
  Legend->SetHeader("Reco Method");
  
  TH1F* XVtxHist;
  TH1F* YVtxHist;
  TH1F* ZVtxHist;
  TH1F* FlashTimeHist;
  TH1F* FlashPEHist;
  
  std::string Version = "v05_06_00";
  
  FileNameVec.push_back("pandoraNuKHit");
  FileNameVec.push_back("pandoraNu");
  FileNameVec.push_back("pandoraCosmic");
  FileNameVec.push_back("pmtrack");
  
//   FileVec.push_back( new TFile(("Hist_pandoraNuKHit_"+Version+".root").c_str(),"READ") );
//   FileVec.push_back( new TFile(("Hist_pandoraNu_"+Version+".root").c_str(),"READ") );
//   FileVec.push_back( new TFile(("Hist_pandoraCosmic_"+Version+".root").c_str(),"READ") );
//   FileVec.push_back( new TFile(("Hist_pmtrack_"+Version+".root").c_str(),"READ") );
  
  
  unsigned int FileNumber = 0;
  
  for(const auto& FileName : FileNameVec)
  {
    TFile* File = new TFile( ("Hist_"+FileName+"_"+Version+".root").c_str(),"READ" ) ;
    File->cd();
    
    if(!FileNumber)
    {
      XVtxHist = (TH1F*)File->Get("X Vertex Distribution");
      YVtxHist = (TH1F*)File->Get("X Vertex Distribution");
      ZVtxHist = (TH1F*)File->Get("X Vertex Distribution");
      FlashTimeHist = (TH1F*)File->Get("Flash Time");
      FlashPEHist = (TH1F*)File->Get("Flash Intensity in PE");
    }
    
//     XVtxHistVec.push_back((TH1F*)File->Get("X Vertex Distribution"));
//     YVtxHistVec.push_back((TH1F*)File->Get("Y Vertex Distribution"));
//     ZVtxHistVec.push_back((TH1F*)File->Get("Z Vertex Distribution"));
    FlashTrackHistVec.push_back((TH1F*)File->Get("Distance from Flash to Track"));
    VertexTrackHistVec.push_back((TH1F*)File->Get("Distance from Vertex to Track"));
    TrackLengthHistVec.push_back((TH1F*)File->Get("Track Length Distribution"));
    TrackMultipHistVec.push_back((TH1F*)File->Get("Track Multiplicity"));
    
    FileNumber++;
  }
  
  TCanvas *Canvas1 = new TCanvas("X Vertex Position", "X Vertex Position", 1400, 1000);
  Canvas1->cd();
  XVtxHist->Draw();
  
  TCanvas *Canvas2 = new TCanvas("Y Vertex Position", "Y Vertex Position", 1400, 1000);
  Canvas2->cd();
  YVtxHist->Draw();
  
  TCanvas *Canvas3 = new TCanvas("Z Vertex Position", "Z Vertex Position", 1400, 1000);
  Canvas3->cd();
  ZVtxHist->Draw();

  TCanvas *Canvas4 = new TCanvas("Flash Time Distribution", "Flash Time Distribution", 1400, 1000);
  Canvas4->cd();
  FlashTimeHist->Draw();

  TCanvas *Canvas5 = new TCanvas("PE Count of Flash", "PE Count of Flash", 1400, 1000);
  Canvas5->cd();
  Canvas5->SetLogy();
  FlashPEHist->Draw();

  TCanvas *Canvas6 = new TCanvas("Flash Track Distance", "Flash Track Distance", 1400, 1000);
  Canvas6->cd();
  Canvas6->SetLogy();
  FlashTrackHistVec.front()->SetMaximum(1.1*GetMaximum(FlashTrackHistVec));
  FlashTrackHistVec.front()->Draw();

  TCanvas *Canvas7 = new TCanvas("Vertex Track Distance", "Vertex Track Distance", 1400, 1000);
  Canvas7->cd();
  Canvas7->SetLogy();
  VertexTrackHistVec.front()->SetMaximum(1.1*GetMaximum(VertexTrackHistVec));
  VertexTrackHistVec.front()->Draw();

  TCanvas *Canvas8 = new TCanvas("Track Length", "Track Length", 1400, 1000);
  Canvas8->cd();
  Canvas8->SetLogy();
  TrackLengthHistVec.front()->SetMaximum(1.1*GetMaximum(TrackLengthHistVec));
  TrackLengthHistVec.front()->Draw();
  
  TCanvas *Canvas9 = new TCanvas("Track Multiplicity", "Track Multiplicity", 1400, 1000);
  Canvas9->cd();
  TrackMultipHistVec.front()->SetMaximum(1.1*GetMaximum(TrackMultipHistVec));
  TrackMultipHistVec.front()->Draw();
  
  std::cout << "Hello" << std::endl;
  
  for(unsigned int no_entries = 1; no_entries < FileNameVec.size(); no_entries++)
  {
//     Canvas1->cd();
//     XVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//     XVtxHistVec.at(no_entries)->Draw("same");
    
//     Canvas2->cd();
//     YVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//     YVtxHistVec.at(no_entries)->Draw("same");
    
//     Canvas3->cd();
//     ZVtxHistVec.at(no_entries)->SetLineColor(no_entries);
//     ZVtxHistVec.at(no_entries)->Draw("same");
    
    Canvas6->cd();
    FlashTrackHistVec.at(no_entries)->SetLineColor(no_entries);
    FlashTrackHistVec.at(no_entries)->Draw("same");
    
    Canvas7->cd();
    VertexTrackHistVec.at(no_entries)->SetLineColor(no_entries);
    VertexTrackHistVec.at(no_entries)->Draw("same");
    
    Canvas8->cd();
    TrackLengthHistVec.at(no_entries)->SetLineColor(no_entries);
    TrackLengthHistVec.at(no_entries)->Draw("same");
    
    Canvas9->cd();
    TrackMultipHistVec.at(no_entries)->SetLineColor(no_entries);
    TrackMultipHistVec.at(no_entries)->Draw("same");
  }
  
  // Fill legend
  for(unsigned int no_entries = 0; no_entries < FileNameVec.size(); no_entries++)
  {
    Legend->AddEntry( FlashTrackHistVec.at(no_entries), (FileNameVec.at(no_entries)).c_str(),"l" );
  }
  
  Canvas1->SaveAs(("XVtxPosition_"+Version+".png").c_str());
  Canvas2->SaveAs(("YVtxPosition_"+Version+".png").c_str());
  Canvas3->SaveAs(("ZVtxPosition_"+Version+".png").c_str());
  Canvas4->SaveAs(("FlashTime_"+Version+".png").c_str());
  Canvas5->SaveAs(("PECount_"+Version+".png").c_str());
  
  Canvas6->cd();
  Legend->Draw();
  Canvas6->SaveAs(("FlashTrackDist_"+Version+".png").c_str());
  
  Canvas7->cd();
  Legend->Draw();
  Canvas7->SaveAs(("VtxTrackDist_"+Version+".png").c_str());
  
  Canvas8->cd();
  Legend->Draw();
  Canvas8->SaveAs(("TrackLength_"+Version+".png").c_str());
  
  Canvas9->cd();
  Legend->Draw();
  Canvas9->SaveAs(("TrackMultiplicity_"+Version+".png").c_str());
  
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
