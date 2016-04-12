#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TTree.h>

using namespace std;

//This defines our current settings for the fiducial volume
double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z) {
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    else return false;
}

//This function returns the distance between a flash and a track (in one dimension, here used only for z direction)
double FlashTrackDist(double flash, double start, double end) {
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
    }
}

// Main function
int runOnMCC7_numuCC_QE_wnumuvtx()
{

    string Version = "v05_08_00";
//     string GeneratorName = "prodgenie_bnb_nu_cosmic_";
    string GeneratorName = "prodgenie_bnb_nu_";

    std::vector<string> TrackProdNameVec;

    TrackProdNameVec.push_back("pandoraNuKHit");
    TrackProdNameVec.push_back("pandoraCosmic");
    TrackProdNameVec.push_back("pandoraNu");
    TrackProdNameVec.push_back("pmtrack");
    TrackProdNameVec.push_back("pandoraNuPMA");
    TrackProdNameVec.push_back("trackkalmanhit");

    std::vector<string> VertexProdNameVec;

    VertexProdNameVec.push_back("nuvtx");
    VertexProdNameVec.push_back("pandoraCosmic");
    VertexProdNameVec.push_back("pandoraNu");
    VertexProdNameVec.push_back("pmtrack");

//    string ProductName = "pandoraNuKHit";
//    string ProductName = "pandoraCosmic";
//    string ProductName = "pandoraNu";
//    string ProductName = "pmtrack";

    TChain *treenc = new TChain("analysistree/anatree");
    treenc -> Add( ("/lheppc46/data/uBData/anatrees/" + GeneratorName + Version +"_anatree.root").c_str() );

    // Open FileStream
    ofstream FinalCutFile;
    ofstream EfficiencyFile;
    ofstream PurityFile;

    FinalCutFile.open("FinalCut_"+GeneratorName+Version+".cvs",ios::trunc);
    EfficiencyFile.open("SelectionEfficiency_"+GeneratorName+Version+".cvs",ios::trunc);
    PurityFile.open("SelectionPurity_"+GeneratorName+Version+".cvs",ios::trunc);

    // Set titles
    FinalCutFile << "Number of events which survive all cuts normalized to 20000 events\n";
    EfficiencyFile << "Selection Efficiency\n" ;
    PurityFile << "Selection Purity\n" ;

    // First legend entry
    FinalCutFile << "track\\vertex";
    EfficiencyFile << "track\\vertex";
    PurityFile << "track\\vertex";

    // Fill all other column legends
    for(const auto& VertexingName : VertexProdNameVec)
    {
        FinalCutFile << "," + VertexingName;
        EfficiencyFile << "," + VertexingName;
        PurityFile << "," + VertexingName;
    }
    FinalCutFile << "\n";
    EfficiencyFile << "\n";
    PurityFile << "\n";

    //maximum array sizes
    const int maxentries = 35000;
    const int maxtracks = 5000;
    const int maxvtx = 500;
    const int maxnu = 10;
    const int maxmc = 10;
    const int kMaxFlashes = 2000;

    //Define variables to read from Tree
    Int_t           event;
    Int_t           run;
    Int_t           subrun;
    Int_t           triggerbits; //this is only important for data. 2048 = BNB stream
    Double_t        potbnb; //this is only important for data

    Short_t         ntracks_reco; //number of reco tracks
    Short_t         trkbestplane[maxtracks]; //plane that has most hits for a given track
    Float_t         trklen[maxtracks]; //track length
    Float_t         trkstartx[maxtracks];
    Float_t         trkstarty[maxtracks];
    Float_t         trkstartz[maxtracks];
    Float_t         trkendx[maxtracks];
    Float_t         trkendy[maxtracks];
    Float_t         trkendz[maxtracks];
    Float_t         trktheta[maxtracks];
    Float_t         trkphi[maxtracks];
    Float_t         trkmomrange[maxtracks]; //track momentum calculated from track range
    Short_t         trkId[maxtracks];
    Short_t         trkorigin[maxtracks][3]; //for MC only: which true particle contributes most hits to the reco track: 2 = cosmic. 1 = neutrino
    bool            vertexatstart[maxtracks]; //for analysis: is the vertex at start of the track?
    bool            vertexatend[maxtracks]; //for analysis: ist the vertex at end of track?

    Int_t           no_flashes; //number of flashes in the event
    Float_t         flash_time[kMaxFlashes]; //flash time (in microseconds)
    Float_t         flash_pe[kMaxFlashes]; //total number of photoelectrons corresponding to the flash
    Float_t         flash_zcenter[kMaxFlashes]; //z center of flash (in cm)
    Float_t         flash_ycenter[kMaxFlashes]; //y center of flash (in cm)

    Short_t         nvtx;
    Float_t         vtxx[maxvtx];
    Float_t         vtxy[maxvtx];
    Float_t         vtxz[maxvtx];

    //finding candidate nu interaction vertex in event
    bool            candvertex[maxvtx];
    Short_t         candtrack[maxvtx];
    Float_t         candlength[maxvtx];
    Short_t         numuvertex = -1;
    Short_t         mutrack = -1;
    Float_t         mutracklength = 0;

    //MC truth
    Int_t           mcevts_truth; //neutrino interactions per event
    Float_t         nuvtxx_truth[maxmc]; //true vertex x (in cm)
    Float_t         nuvtxy_truth[maxmc];
    Float_t         nuvtxz_truth[maxmc];
    Int_t           ccnc_truth[maxmc]; //CC = 0, NC = 1
    Int_t           nuPDG_truth[maxmc]; //true neutrino pdg code. numu = 14
    Int_t           mode_truth[maxmc]; //QE = 0, RES = 1, DIS = 2
    Int_t	   PDG_truth[maxtracks];

    Int_t	   NumberOfMCTracks;

    Float_t	   XMCTrackStart[maxtracks];
    Float_t	   YMCTrackStart[maxtracks];
    Float_t	   ZMCTrackStart[maxtracks];

    Float_t	   XMCTrackEnd[maxtracks];
    Float_t	   YMCTrackEnd[maxtracks];
    Float_t	   ZMCTrackEnd[maxtracks];


    //define cut variables
    double flashwidth = 80; //cm. Distance flash-track
    double distcut = 5; //cm. Distance track start/end to vertex
    double lengthcut = 75; //cm. Length of longest track
    double beammin = 3.55; //us. Beam window start
    double beammax = 5.15; //us. Beam window end
    double PEthresh = 50; //PE
    double MCTrackToMCVtxDist = 1; //cm. distance between mc track start and mc vertex
    double TrackToMCDist = 3; //cm. Distance track start/end to mcvertex

    // Loop over all product names
    for(const auto& TrackingName : TrackProdNameVec)
    {
        // Write tracking methode to files first column
        FinalCutFile << TrackingName;
        EfficiencyFile << TrackingName;
        PurityFile << TrackingName;

        for(const auto& VertexingName : VertexProdNameVec)
        {
            treenc -> SetBranchAddress("event", &event);
            treenc -> SetBranchAddress("no_flashes", &no_flashes);
            treenc -> SetBranchAddress("flash_time", flash_time);
            treenc -> SetBranchAddress("flash_pe", flash_pe);
            treenc -> SetBranchAddress("flash_zcenter", flash_zcenter);
            treenc -> SetBranchAddress("flash_ycenter", flash_ycenter);
            treenc -> SetBranchAddress("mcevts_truth", &mcevts_truth);
            treenc -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
            treenc -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
            treenc -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
            treenc -> SetBranchAddress("ccnc_truth", ccnc_truth);
            treenc -> SetBranchAddress("nuPDG_truth", nuPDG_truth);
            treenc -> SetBranchAddress("pdg", PDG_truth);
            treenc -> SetBranchAddress("mode_truth", mode_truth);
            treenc -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
            treenc -> SetBranchAddress("StartPointx", XMCTrackStart);
            treenc -> SetBranchAddress("StartPointy", YMCTrackStart);
            treenc -> SetBranchAddress("StartPointz", ZMCTrackStart);
            treenc -> SetBranchAddress("EndPointx", XMCTrackEnd);
            treenc -> SetBranchAddress("EndPointy", YMCTrackEnd);
            treenc -> SetBranchAddress("EndPointz", ZMCTrackEnd);

            // Product specific stuff
            treenc -> SetBranchAddress(("ntracks_"+TrackingName).c_str(),&ntracks_reco);
            treenc -> SetBranchAddress(("trkstartx_"+TrackingName).c_str(),trkstartx);
            treenc -> SetBranchAddress(("trkstarty_"+TrackingName).c_str(),trkstarty);
            treenc -> SetBranchAddress(("trkstartz_"+TrackingName).c_str(),trkstartz);
            treenc -> SetBranchAddress(("trkendx_"+TrackingName).c_str(),trkendx);
            treenc -> SetBranchAddress(("trkendy_"+TrackingName).c_str(),trkendy);
            treenc -> SetBranchAddress(("trkendz_"+TrackingName).c_str(),trkendz);
            treenc -> SetBranchAddress(("trktheta_"+TrackingName).c_str(),trktheta);
            treenc -> SetBranchAddress(("trkphi_"+TrackingName).c_str(),trkphi);
            treenc -> SetBranchAddress(("trkorigin_"+TrackingName).c_str(),trkorigin);
            // Program hack to apply for non uniform naming of nuvtx
            if(VertexingName != "nuvtx")
            {
                treenc -> SetBranchAddress(("nvtx_"+VertexingName).c_str(), &nvtx);
                treenc -> SetBranchAddress(("vtxx_"+VertexingName).c_str(), vtxx);
                treenc -> SetBranchAddress(("vtxy_"+VertexingName).c_str(), vtxy);
                treenc -> SetBranchAddress(("vtxz_"+VertexingName).c_str(), vtxz);
            }
            else
            {
                treenc -> SetBranchAddress("nnuvtx", &nvtx);
                treenc -> SetBranchAddress("nuvtxx", vtxx);
                treenc -> SetBranchAddress("nuvtxy", vtxy);
                treenc -> SetBranchAddress("nuvtxz", vtxz);
            }

            // Vertex position calues
            TH1F *hXVertexPosition = new TH1F("X Vertex Distribution","X Vertex Distribution",256,0,FVx);
            hXVertexPosition->SetStats(0);
            hXVertexPosition->GetXaxis()->SetTitle("X position [cm]");
            hXVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
            TLine XVertexMinCut(borderx, 0, borderx, 200);
            TLine XVertexMaxCut(FVx-borderx, 0, FVx-borderx, 200);
            XVertexMinCut.SetLineColor(kRed);
            XVertexMaxCut.SetLineColor(kRed);

            TH1F *hYVertexPosition = new TH1F("Y Vertex Distribution","Y Vertex Distribution",333,-FVy/2,(FVy+100)/2);
            hYVertexPosition->SetStats(0);
            hYVertexPosition->GetXaxis()->SetTitle("Y position [cm]");
            hYVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
            TLine YVertexMinCut(-FVy/2+bordery, 0, -FVy/2+bordery, 200);
            TLine YVertexMaxCut(FVy/2-bordery, 0, FVy/2-bordery, 200);
            YVertexMinCut.SetLineColor(kRed);
            YVertexMaxCut.SetLineColor(kRed);

            TH1F *hZVertexPosition = new TH1F("Z Vertex Distribution","Z Vertex Distribution",1000,0,FVz);
            hZVertexPosition->SetStats(0);
            hZVertexPosition->GetXaxis()->SetTitle("Z position [cm]");
            hZVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
            TLine ZVertexMinCut(borderz, 0, borderz, 200);
            TLine ZVertexMaxCut(FVz-borderz, 0, FVz-borderz, 200);
            ZVertexMinCut.SetLineColor(kRed);
            ZVertexMaxCut.SetLineColor(kRed);

            // Flash Time distribution
            TH1F *hFlashTime = new TH1F("Flash Time","Flash Time",200,0,8.2);
            hFlashTime->SetStats(0);
            hFlashTime->GetXaxis()->SetTitle("Time [#mus]");
            hFlashTime->GetYaxis()->SetTitle("Number of Flashes [ ]");
            TLine FlashTimeMinCut(beammin, 0, beammin, 2000);
            TLine FlashTimeMaxCut(beammax, 0, beammax, 2000);
            FlashTimeMinCut.SetLineColor(kRed);
            FlashTimeMaxCut.SetLineColor(kRed);

            // Flash PE distirbution
            TH1F *hIntesityPE = new TH1F("Flash Intensity in PE","Flash Intensity in PE",100,0,1000);
            hIntesityPE->SetStats(0);
            hIntesityPE->GetXaxis()->SetTitle("Number of P.E. [ ]");
            hIntesityPE->GetYaxis()->SetTitle("Number of Flashes [ ]");
            TLine FlashPEMinCut(PEthresh, 0, PEthresh, 20000);
            FlashPEMinCut.SetLineColor(kRed);

            // Flash to track distance
            TH1F *hFlashTrackDist = new TH1F("Distance from Flash to Track","Distance between Flash and Track",500,0,500);
            hFlashTrackDist->SetStats(0);
            hFlashTrackDist->GetXaxis()->SetTitle("Distance [cm]");
            hFlashTrackDist->GetYaxis()->SetTitle("Number of Events [ ]");
            TLine FlashTrackMinCut(flashwidth, 0, flashwidth, 3000);
            FlashTrackMinCut.SetLineColor(kRed);

            // Vertex to Track distance
            TH1F *hVertexTrackDist = new TH1F("Distance from Vertex to Track","Distance between Vertex and Track",200,0,100);
            hVertexTrackDist->SetStats(0);
            hVertexTrackDist->GetXaxis()->SetTitle("Distance [cm]");
            hVertexTrackDist->GetYaxis()->SetTitle("Number of Tracks [ ]");
            TLine VertexTrackMinCut(distcut, 0, distcut, 3000);
            VertexTrackMinCut.SetLineColor(kRed);

            // TrackLength of longest tracks
            TH1I *hTrackLength = new TH1I("Track Length Distribution","Track Length Distribution",1000,0,1000);
            hTrackLength->SetStats(0);
            hTrackLength->GetXaxis()->SetTitle("Track Length [cm]");
            hTrackLength->GetYaxis()->SetTitle("Number of Tracks [ ]");

            // Track Multiplicity at Vertex
            TH1F *hTrackMultip = new TH1F("Track Multiplicity","Track Multiplicity at the Vertex",10,0,10);
            hTrackMultip->SetStats(0);
            hTrackMultip->GetXaxis()->SetTitle("Number of Tracks [ ]");
            hTrackMultip->GetYaxis()->SetTitle("Number of Events [ ]");

            // Track Vertex Distances
            TH1F *hFlashVertexDist = new TH1F("Flash To Vertex Distance","Flash To Vertex Distance",1000,0,1000);
            hFlashVertexDist->SetStats(0);
            hFlashVertexDist->GetXaxis()->SetTitle("Distance [cm]");
            hFlashVertexDist->GetYaxis()->SetTitle("Number of Vertices [ ]");

            // Track True Vertex Distance
            TH1F *hTrackMCVertexDist = new TH1F("Track To MC-Vertex Distance","Track To MC-Vertex Distance",100,0,50);
            hFlashVertexDist->SetStats(0);
            hFlashVertexDist->GetXaxis()->SetTitle("Distance [cm]");
            hFlashVertexDist->GetYaxis()->SetTitle("Number of Tracks [ ]");

            // Track Start and End postion
            TH1F *hXTrackStart = new TH1F("X Track Start Position","X Track Start Position",256,0,FVx);
            hXTrackStart->GetXaxis()->SetTitle("X position [cm]");
            hXTrackStart->GetYaxis()->SetTitle("Number of Tracks [ ]");

            TH1F *hYTrackStart = new TH1F("Y Track Start Position","Y Track Start Position",233,-FVy/2,FVy/2);
            hYTrackStart->GetXaxis()->SetTitle("Y position [cm]");
            hYTrackStart->GetYaxis()->SetTitle("Number of Tracks [ ]");

            TH1F *hZTrackStart = new TH1F("Z Track Start Position","Z Track Start Position",1000,0,FVz);
            hZTrackStart->GetXaxis()->SetTitle("Z position [cm]");
            hZTrackStart->GetYaxis()->SetTitle("Number of Tracks [ ]");

            TH1F *hXTrackEnd = new TH1F("X Track End Position","X Track End Position",256,0,FVx);
            hXTrackEnd->GetXaxis()->SetTitle("X position [cm]");
            hXTrackEnd->GetYaxis()->SetTitle("Number of Tracks [ ]");

            TH1F *hYTrackEnd = new TH1F("Y Track End Position","Y Track End Position",233,-FVy/2,FVy/2);
            hYTrackEnd->GetXaxis()->SetTitle("Y position [cm]");
            hYTrackEnd->GetYaxis()->SetTitle("Number of Tracks [ ]");

            TH1F *hZTrackEnd = new TH1F("Z Track End Position","Z Track End Position",1000,0,FVz);
            hZTrackEnd->GetXaxis()->SetTitle("Z position [cm]");
            hZTrackEnd->GetYaxis()->SetTitle("Number of Tracks [ ]");

            int Size = treenc -> GetEntries();
            if(Size > 20000) Size = 20000;

            cout << "number of events used is: " << Size << endl;

            int theflash = -1;

            double diststart = 0;
            double distend = 0;
            double length = 0;

            int ntrue = 0;

            unsigned int EventsWithFlash = 0;
            unsigned int EventsVtxInFV = 0;
            unsigned int EventsTrackNearVertex = 0;
            unsigned int EventsFlashMatched = 0;
            unsigned int EventsTracksInFV = 0;
            unsigned int EventsLengthCut = 0;
            unsigned int EventsNearVtx = 0;
            unsigned int EventsTrackLong = 0;
            unsigned int EventsTruelyReco = 0;

            unsigned int NumberOfContainedMCTracks = 0;

            //Event Loop
            for(int i = 0; i < Size; i++)
            {
//                 if(i%1000 == 0) cout << "\t... " << i << endl;

                treenc -> GetEntry(i);

                int MCTrackCandidate = -1;
                float MCTrackCandLength = 0;

                // Loop over all MC neutrino vertices
                for(unsigned vertex_no = 0; vertex_no < mcevts_truth; vertex_no++)
                {
                    // Loop over all MC particles
                    for(unsigned track_no = 0; track_no < NumberOfMCTracks; track_no++)
                    {
                        // Calculate MCTrack length
                        float MCTrackLength = sqrt(pow(XMCTrackStart[track_no] - XMCTrackEnd[track_no],2) + pow(YMCTrackStart[track_no] - YMCTrackEnd[track_no],2) + pow(ZMCTrackStart[track_no] - ZMCTrackEnd[track_no],2));

                        // If the a muon is not contained in a singel neutrino event, set mc-track contained flag to false
                        if( abs(PDG_truth[track_no]) == 13 // Track has to be a muon
                                && inFV(nuvtxx_truth[0],nuvtxy_truth[0],nuvtxz_truth[0]) // true vertex has to be in FV
                                && inFV(XMCTrackStart[track_no],YMCTrackStart[track_no],ZMCTrackStart[track_no]) // Track start has to be in FV
                                && inFV(XMCTrackEnd[track_no],YMCTrackEnd[track_no],ZMCTrackEnd[track_no]) // Track end has to be in FV
                                && sqrt(pow(XMCTrackStart[track_no] - nuvtxx_truth[vertex_no],2) + pow(YMCTrackStart[track_no] - nuvtxy_truth[vertex_no],2) + pow(ZMCTrackStart[track_no] - nuvtxz_truth[vertex_no],2)) < MCTrackToMCVtxDist // Track has to start at vertex
                                && MCTrackLength > lengthcut // Track has to be long
                                && MCTrackLength > MCTrackCandLength // If the current candidate length is shorter than the new length
                          )
                        {
                            // Fill new length and candidate index
                            MCTrackCandLength = MCTrackLength;
                            MCTrackCandidate = track_no;
                        }
                    } // MC particle loop
                } // MC vertex loop

                // Count up the number of contained mc-tracks if there are mc candidates
                if(MCTrackCandidate > -1)
                {
                    NumberOfContainedMCTracks++;
                }

                bool flashtag = false;
                float flashmax = 0;

                // Loop over all flashes
                for(int f = 0; f < no_flashes; f++)
                {
                    // Fill flash related histograms
                    hFlashTime->Fill(flash_time[f]);
                    hIntesityPE->Fill(flash_pe[f]);

                    // If the flash is in the beam window and above threshold set flashtag to true
                    if( (flash_time[f] > beammin && flash_time[f] < beammax) && flash_pe[f] > PEthresh )
                    {
                        flashtag = true; //the event does have a flash inside the beam window

                        // If the new flash has more PE than the current maximum, replace the maximum
                        if(flash_pe[f] > flashmax)
                        {
                            theflash = f;
                            flashmax = flash_pe[f];
                        }
                    }
                } // flash loop

                // If the flash tag is ture
                if(flashtag)
                {
                    // Prepare flags
                    bool VertexInFVFlag = true;
                    bool TrackDistanceFlag = true;
                    bool FlashMatchFlag = true;
                    bool TrackContainedFlag = false;

                    // Increase events with flash > 50 PE and within beam window
                    EventsWithFlash++;

                    // Loop over all reconstructed tracks
                    for(int j = 0; j < ntracks_reco; j++)
                    {
                        // Fill track length histogram and the flash track histogram
                        hTrackLength->Fill( pow(trkstartx[j] - trkendx[j],2) + pow(trkstarty[j] - trkendy[j],2) + pow(trkstartz[j] - trkendz[j],2) );
                        hFlashTrackDist->Fill(FlashTrackDist(flash_zcenter[theflash], trkstartz[j], trkendz[j]));
                    } // reco track loop

                    // Initialize Track Candidate properties
                    int TrackCandidate = -1;
                    float TrackCandLength = 0;

                    // Loop over all vertices
                    for(int v = 0; v < nvtx; v++)
                    {
                        // Fill vertex position to histogram
                        hXVertexPosition->Fill(vtxx[v]);
                        hYVertexPosition->Fill(vtxy[v]);
                        hZVertexPosition->Fill(vtxz[v]);

                        // If the vertex is contained
                        if(inFV(vtxx[v], vtxy[v], vtxz[v]))
                        {
                            // Increase count of events with a vertex in FV
                            if(VertexInFVFlag)
                            {
                                EventsVtxInFV++;
                                VertexInFVFlag = false;
                            }

                            unsigned int TrackCountAtVertex = 0;

                            // Loop over all reconstructed tracks
                            for(int j = 0; j < ntracks_reco; j++)
                            {
                                // Calculate distances from track start/end to vertex and calculate track lenth
                                diststart = sqrt((vtxx[v] - trkstartx[j])*(vtxx[v] - trkstartx[j]) + (vtxy[v] - trkstarty[j])*(vtxy[v] - trkstarty[j]) + (vtxz[v] - trkstartz[j])*(vtxz[v] - trkstartz[j]));
                                distend = sqrt((vtxx[v] - trkendx[j])*(vtxx[v] - trkendx[j]) + (vtxy[v] - trkendy[j])*(vtxy[v] - trkendy[j]) + (vtxz[v] - trkendz[j])*(vtxz[v] - trkendz[j]));
                                length = sqrt(pow(trkstartx[j] - trkendx[j],2) + pow(trkstarty[j] - trkendy[j],2) + pow(trkstartz[j] - trkendz[j],2));

                                if(diststart < distend)
                                {
                                    hVertexTrackDist->Fill(diststart);
                                }
                                else
                                {
                                    hVertexTrackDist->Fill(distend);
                                }

                                // If the track vertex distance is within cut, increase track count
                                if(diststart < distcut || distend < distcut)
                                {
                                    // Increase count of events with a distance to vertex < 5 cm
                                    if(TrackDistanceFlag)
                                    {
                                        EventsTrackNearVertex++;
                                        TrackDistanceFlag = false;
                                    }

                                    // Increase track at vertex count
                                    TrackCountAtVertex++;

                                    // Find the longest track from this vertex
                                    if(length > TrackCandLength)
                                    {
                                        TrackCandLength = length;
                                        TrackCandidate = j;
                                    }
                                } // if track vertex distance is within cut distance
                            } // reco track loop

                            // Fill vertex related histograms
                            hTrackMultip->Fill(TrackCountAtVertex);
                            hFlashVertexDist->Fill(fabs(flash_zcenter[theflash]-vtxz[v]));
                        } // if vertex is contained
                    } // vertex loop
                    
                    // If the longest track length is filled
                    if(TrackCandidate > -1)
                    {
                        // If the longest track is flash matched
                        if( FlashTrackDist(flash_zcenter[theflash], trkstartz[TrackCandidate], trkendz[TrackCandidate]) < flashwidth )
                        {
                            if(FlashMatchFlag)
                            {
                                EventsFlashMatched++;
                                FlashMatchFlag = false;
                                // Set track contained flag true, so the other cuts can be applied on this vertex
                                TrackContainedFlag = true;
                            }

                            // If the longest track is fully contained
                            if( inFV(trkstartx[TrackCandidate], trkstarty[TrackCandidate], trkstartz[TrackCandidate])
                                    && inFV(trkendx[TrackCandidate], trkendy[TrackCandidate], trkendz[TrackCandidate])
                                    && TrackContainedFlag)
                            {
                                EventsTracksInFV++;

                                // If longest track is longer than 75 cm
                                if(TrackCandLength > lengthcut)
                                {
                                    EventsTrackLong++;
                                    
                                    double TrkStartMCStartDist = sqrt(pow(XMCTrackStart[MCTrackCandidate] - trkstartx[TrackCandidate],2) + pow(YMCTrackStart[MCTrackCandidate] - trkstarty[TrackCandidate],2) + pow(ZMCTrackStart[MCTrackCandidate] - trkstartz[TrackCandidate],2));
                                    double TrkEndMCEndDist = sqrt(pow(XMCTrackEnd[MCTrackCandidate] - trkendx[TrackCandidate],2) + pow(YMCTrackEnd[MCTrackCandidate] - trkendy[TrackCandidate],2) + pow(ZMCTrackEnd[MCTrackCandidate] - trkendz[TrackCandidate],2));
                                    double TrkStartMCEndDist = sqrt(pow(XMCTrackEnd[MCTrackCandidate] - trkstartx[TrackCandidate],2) + pow(YMCTrackEnd[MCTrackCandidate] - trkstarty[TrackCandidate],2) + pow(ZMCTrackEnd[MCTrackCandidate] - trkstartz[TrackCandidate],2));
                                    double TrkEndMCStartDist = sqrt(pow(XMCTrackStart[MCTrackCandidate] - trkstartx[TrackCandidate],2) + pow(YMCTrackStart[MCTrackCandidate] - trkstarty[TrackCandidate],2) + pow(ZMCTrackStart[MCTrackCandidate] - trkstartz[TrackCandidate],2));

//                                     if(TrkStartMCDist < TrkEndMCDist)
//                                     {
//                                         hTrackMCVertexDist->Fill(TrkStartMCDist);
//                                     }
//                                     else
//                                     {
//                                         hTrackMCVertexDist->Fill(TrkEndMCDist);
//                                     }

                                    // If track end or start are close to montecarlo vertex
                                    if(   (TrkStartMCStartDist < TrackToMCDist && TrkEndMCEndDist < TrackToMCDist)
                                        ||(TrkStartMCEndDist < TrackToMCDist && TrkEndMCStartDist < TrackToMCDist)
                                      )
                                    {
                                        EventsTruelyReco++;
                                    }
                                }
                            }
                            // Set track contained flag false
                            TrackContainedFlag = false;
                        }
                    } // If there is a longest track
                } // if flashtag

                // Increase the neutrino count
                ntrue++;
            }//loop over all events

            TCanvas *Canvas1 = new TCanvas("X Vertex Position", "X Vertex Position", 1400, 1000);
            hXVertexPosition->Draw();
            XVertexMinCut.Draw("same");
            XVertexMaxCut.Draw("same");

            TCanvas *Canvas2 = new TCanvas("Y Vertex Position", "Y Vertex Position", 1400, 1000);
            hYVertexPosition->Draw();
            YVertexMinCut.Draw("same");
            YVertexMaxCut.Draw("same");

            TCanvas *Canvas3 = new TCanvas("Z Vertex Position", "Z Vertex Position", 1400, 1000);
            hZVertexPosition->Draw();
            ZVertexMinCut.Draw("same");
            ZVertexMaxCut.Draw("same");

            TCanvas *Canvas4 = new TCanvas("Flash Time Distribution", "Flash Time Distribution", 1400, 1000);
            hFlashTime->Draw();
            FlashTimeMinCut.Draw("same");
            FlashTimeMaxCut.Draw("same");

            TCanvas *Canvas5 = new TCanvas("PE Count of Flash", "PE Count of Flash", 1400, 1000);
            Canvas5->SetLogy();
            hIntesityPE->Draw();
            FlashPEMinCut.Draw("same");

            TCanvas *Canvas6 = new TCanvas("Flash Track Distance", "Flash Track Distance", 1400, 1000);
            Canvas6->SetLogy();
            hFlashTrackDist->Draw();
            FlashTrackMinCut.Draw("same");

            TCanvas *Canvas7 = new TCanvas("Vertex Track Distance", "Vertex Track Distance", 1400, 1000);
            Canvas7->SetLogy();
            hVertexTrackDist->Draw();
            VertexTrackMinCut.Draw("same");

            TCanvas *Canvas8 = new TCanvas("Track Length", "Track Length", 1400, 1000);
            Canvas8->SetLogy();
            hTrackLength->Draw();

            TCanvas *Canvas9 = new TCanvas("Track Multiplicity", "Track Multiplicity", 1400, 1000);
            hTrackMultip->Draw();

            TCanvas *Canvas10 = new TCanvas("Flash Vertex Distance", "Flash Vertex Distance", 1400, 1000);
            hFlashVertexDist->Draw();

            TCanvas *Canvas11 = new TCanvas("Track MCVertex Distance", "Track MCVertex Distance", 1400, 1000);
            hTrackMCVertexDist->Draw();

            TFile* OutputFile = new TFile(("Hist_Track_"+TrackingName+ "_Vertex_" + VertexingName + "_"+GeneratorName+Version+".root").c_str(),"RECREATE");

            Canvas1->Write();
            Canvas2->Write();
            Canvas3->Write();
            Canvas4->Write();
            Canvas5->Write();
            Canvas6->Write();
            Canvas7->Write();
            Canvas8->Write();
            Canvas10->Write();
            Canvas11->Write();

            hXVertexPosition->Write();
            hYVertexPosition->Write();
            hZVertexPosition->Write();
            hFlashTime->Write();
            hIntesityPE->Write();
            hFlashTrackDist->Write();
            hVertexTrackDist->Write();
            hTrackLength->Write();
            hTrackMultip->Write();
            hFlashVertexDist->Write();
            hTrackMCVertexDist->Write();

            OutputFile->Close();

            Canvas10->SaveAs(("FlashVertexDist"+GeneratorName+Version+".png").c_str());
            Canvas11->SaveAs(("TrackMCVertexDist_Track_"+ TrackingName + "_Vertex_" + VertexingName + "_" +GeneratorName+Version+".png").c_str());
//    Canvas2->SaveAs((TrackingName+"_YVtxPosition_"+GeneratorName+Version+".png").c_str());
//    Canvas3->SaveAs((TrackingName+"_ZVtxPosition_"+GeneratorName+Version+".png").c_str());
//    Canvas4->SaveAs((TrackingName+"_FlashTime_"+GeneratorName+Version+".png").c_str());
//    Canvas5->SaveAs((TrackingName+"_PECount_"+GeneratorName+Version+".png").c_str());
//    Canvas6->SaveAs((TrackingName+"_FlashTrackDist_"+GeneratorName+Version+".png").c_str());
//    Canvas7->SaveAs((TrackingName+"_VtxTrackDist_"+GeneratorName+Version+".png").c_str());
//    Canvas8->SaveAs((TrackingName+"_TrackLength_"+GeneratorName+Version+".png").c_str());

            cout << "--------------------------------------------------------------------------------------------" << endl;
            cout << endl;
            cout << "Track Reco Product Name : " << TrackingName << "  Vertex Reco Product Name : " << VertexingName << endl;
            cout << "number of CC events with vertex in FV : " << ntrue << endl;
            cout << "number of events with flash > 50 PE : " << EventsWithFlash << endl;
            cout << "number of events with vtx in FV : " << EventsVtxInFV << endl;
            cout << "number of events with track start/end within 5cm to vtx : " << EventsTrackNearVertex << endl;
            cout << "number of events with tracks matched within 80cm to flash : " << EventsFlashMatched << endl;
            cout << "number of events with contained tracks : " << EventsTracksInFV << endl;
            cout << "number of events with longest track > 75cm : " << EventsTrackLong << endl;
            cout << "number of events with track start end within 3cm to mc-vtx : " << EventsTruelyReco << endl;
            cout << "number of events with contained MC tracks : " << NumberOfContainedMCTracks << endl;
            cout << "event selection efficiency : " <<  (float)EventsTruelyReco/(float)NumberOfContainedMCTracks << endl;
            cout << "event selection purity : " <<  (float)EventsTruelyReco/(float)EventsTrackLong << endl;
//             cout << "event selection missid rate : " <<  fabs((float)EventsTruelyReco-(float)NumberOfContainedMCTracks)/(float)NumberOfContainedMCTracks << endl;
            cout << endl;
            cout << "--------------------------------------------------------------------------------------------" << endl;

            // Write numbers to cvs File
            FinalCutFile << "," << (float)EventsTrackLong/(float)Size*20000.;
            EfficiencyFile << "," << (float)EventsTruelyReco/(float)NumberOfContainedMCTracks;
            PurityFile << "," << (float)EventsTruelyReco/(float)EventsTrackLong;

            // Garbage collection
            delete Canvas1;
            delete Canvas2;
            delete Canvas3;
            delete Canvas4;
            delete Canvas5;
            delete Canvas6;
            delete Canvas7;
            delete Canvas8;
            delete Canvas9;
            delete Canvas10;
            delete Canvas11;

            delete hXVertexPosition;
            delete hYVertexPosition;
            delete hZVertexPosition;
            delete hFlashTime;
            delete hIntesityPE;
            delete hFlashTrackDist;
            delete hVertexTrackDist;
            delete hTrackLength;
            delete hTrackMultip;
            delete hFlashVertexDist;
            delete hTrackMCVertexDist;
            delete hXTrackStart;
            delete hYTrackStart;
            delete hZTrackStart;
            delete hXTrackEnd;
            delete hYTrackEnd;
            delete hZTrackEnd;

            // Erase all branch addresses for the next iteration
            treenc -> ResetBranchAddresses();

        } // Loop over all vertexing data products

        // Get to next line in file
        FinalCutFile << "\n";
        EfficiencyFile << "\n";
        PurityFile << "\n";
    } // Loop over all tracking data products

    FinalCutFile.close();
    EfficiencyFile.close();
    PurityFile.close();

    return 0;

} // end main function
