#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TTree.h>

void ReadBeamErrors()
{
    // Number of columns in file
    unsigned short NumberOfColumns = 5;

    // Initialize data structure
    std::vector <std::vector<float>> BeamSystematics;
    BeamSystematics.resize(NumberOfColumns);
    
    // Line and cell string for ifstream
    std::string FileLine;
    std::string Cell;
    
    std::string FileType = "png";

    // Open field systematic error file
    std::ifstream SysFile ("bnb_sys_error_uboone.txt");

    // check file
    if(SysFile.bad())
    {
        std::cout << "No such file or directory: " << "bnb_sys_error_uboone.txt" << std::endl;
        return;
    }

    // Loop over lines until files end
    while(std::getline(SysFile,FileLine))
    {
        // If not a header line
        if(FileLine[0] != 'E')
        {
            // First column entry gets filled with a factor 1000 for conversion to MeV
            if(SysFile >> Cell)
            {
                BeamSystematics.at(0).push_back(1000*std::stof(Cell));
            }
            
            // Loop over all remaining columns
            for(unsigned column_no = 1; column_no < NumberOfColumns; column_no++)
            {
                // Only read data if data stream works
                if(SysFile >> Cell)
                {
                    // Fill systematic error data into data structure
                    BeamSystematics.at(column_no).push_back(std::stof(Cell));
                }
            } // End of column loop
        } // if not header
    } // line loop
    
    TGraph *NuMuGraph = new TGraph(BeamSystematics.at(0).size(),BeamSystematics.at(0).data(),BeamSystematics.at(1).data());
    NuMuGraph->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
    NuMuGraph->GetYaxis()->SetTitle("Relative Systematic Error");
    TGraph *NuMuBarGraph = new TGraph(BeamSystematics.at(0).size(),BeamSystematics.at(0).data(),BeamSystematics.at(2).data());
    NuMuBarGraph->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
    NuMuBarGraph->GetYaxis()->SetTitle("Relative Systematic Error");
    TGraph *NueGraph = new TGraph(BeamSystematics.at(0).size(),BeamSystematics.at(0).data(),BeamSystematics.at(3).data());
    NueGraph->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
    NueGraph->GetYaxis()->SetTitle("Relative Systematic Error");
    TGraph *NueBarGraph = new TGraph(BeamSystematics.at(0).size(),BeamSystematics.at(0).data(),BeamSystematics.at(4).data());
    NueBarGraph->GetXaxis()->SetTitle("Neutrino Energy [MeV]");
    NueBarGraph->GetYaxis()->SetTitle("Relative Systematic Error");
    
    NuMuGraph->SetTitle("#nu_{#mu} Beam Systematics");
    NuMuBarGraph->SetTitle("#bar{#nu}_{#mu} Beam Systematics");
    NueGraph->SetTitle("#nu_{e} Beam Systematics");
    NueBarGraph->SetTitle("#bar{#nu}_{e} Beam Systematics");
    
    TSpline3 *Spline3NuMu = new TSpline3("Spline3NuMu",NuMuGraph);
    TSpline5 *Spline5NuMu = new TSpline5("Spline5NuMu",NuMuGraph);
    TSpline3 *Spline3NuMuBar = new TSpline3("Spline3NuMuBar",NuMuBarGraph);
    TSpline5 *Spline5NuMuBar = new TSpline5("Spline5NuMuBar",NuMuBarGraph);
    TSpline3 *Spline3Nue = new TSpline3("Spline3Nue",NueGraph);
    TSpline5 *Spline5Nue = new TSpline5("Spline5Nue",NueGraph);
    TSpline3 *Spline3NueBar = new TSpline3("Spline3NueBar",NueBarGraph);
    TSpline5 *Spline5NueBar = new TSpline5("Spline5NueBar",NueBarGraph);
    
    Spline3NuMu->SetLineColor(2);
    Spline5NuMu->SetLineColor(4);
    Spline3NuMuBar->SetLineColor(2);
    Spline5NuMuBar->SetLineColor(4);
    Spline3Nue->SetLineColor(2);
    Spline5Nue->SetLineColor(4);
    Spline3NueBar->SetLineColor(2);
    Spline5NueBar->SetLineColor(4);
    
    TLegend* SplineLegend = new TLegend(0.7,0.8,0.9,0.9);
    SplineLegend->AddEntry( Spline3NuMu, "Spline 3","l" );
    SplineLegend->AddEntry( Spline5NuMu, "Spline 5","l" );
    
    
    TCanvas *Canvas0 = new TCanvas("NuMu", "NuMu", 1400, 1000);
    NuMuGraph->Draw("A*");
    Spline3NuMu->Draw("same");
    Spline5NuMu->Draw("same");
    SplineLegend->Draw("same");
    Canvas0->SaveAs(("NuMuSys."+FileType).c_str());
    
    TCanvas *Canvas1 = new TCanvas("NuMuBar", "NuMuBar", 1400, 1000);
    NuMuBarGraph->Draw("A*");
    Spline3NuMuBar->Draw("same");
    Spline5NuMuBar->Draw("same");
    SplineLegend->Draw("same");
    Canvas1->SaveAs(("NuMuBarSys."+FileType).c_str());
    
    TCanvas *Canvas2 = new TCanvas("Nue", "Nue", 1400, 1000);
    NueGraph->Draw("A*");
    Spline3Nue->Draw("same");
    Spline5Nue->Draw("same");
    SplineLegend->Draw("same");
    Canvas2->SaveAs(("NueSys."+FileType).c_str());
    
    TCanvas *Canvas3 = new TCanvas("NueBar", "NueBar", 1400, 1000);
    NueBarGraph->Draw("A*");
    Spline3NueBar->Draw("same");
    Spline5NueBar->Draw("same");
    SplineLegend->Draw("same");
    Canvas3->SaveAs(("NueBarSys."+FileType).c_str());

} // end main function
