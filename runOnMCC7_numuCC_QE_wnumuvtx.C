#include <vector>

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

int runOnMCC7_numuCC_QE_wnumuvtx() {
  
   string Version = "v05_06_00";
   
   string ProductName = "pandoraNuKHit";
//    string ProductName = "pandoraCosmic";
//    string ProductName = "pandoraNu";
//    string ProductName = "pmtrack";
   
   cout << ProductName << endl;

   TChain *treenc = new TChain("analysistree/anatree");
   treenc -> Add( ("/media/christoph/200EFBDA63AA160B/anatrees/prodgenie_bnb_nu_uboone_ana_"+Version+".root").c_str() );

   //maximum array sizes
   const int maxentries = 35000;
   const int maxtracks = 2000;
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

   Int_t           ntracks_reco; //number of reco tracks
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
   
   treenc -> SetBranchAddress("event", &event);
   treenc -> SetBranchAddress(("ntracks_"+ProductName).c_str(),&ntracks_reco);
   treenc -> SetBranchAddress(("trkstartx_"+ProductName).c_str(),trkstartx);
   treenc -> SetBranchAddress(("trkstarty_"+ProductName).c_str(),trkstarty);
   treenc -> SetBranchAddress(("trkstartz_"+ProductName).c_str(),trkstartz);
   treenc -> SetBranchAddress(("trkendx_"+ProductName).c_str(),trkendx);
   treenc -> SetBranchAddress(("trkendy_"+ProductName).c_str(),trkendy);
   treenc -> SetBranchAddress(("trkendz_"+ProductName).c_str(),trkendz);
   treenc -> SetBranchAddress(("trktheta_"+ProductName).c_str(),trktheta);
   treenc -> SetBranchAddress(("trkphi_"+ProductName).c_str(),trkphi);
   treenc -> SetBranchAddress(("trkorigin_"+ProductName).c_str(),trkorigin);
   treenc -> SetBranchAddress("nnuvtx", &nvtx);
   treenc -> SetBranchAddress("nuvtxx", vtxx);
   treenc -> SetBranchAddress("nuvtxy", vtxy);
   treenc -> SetBranchAddress("nuvtxz", vtxz);
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
   
   
   //Angular distribution of events passing the cuts for origin = cosmic
   TH1F *hc_theta = new TH1F("hc_theta", "cos(#theta)", 20, -1, 1);
   TH1F *hc_phi = new TH1F("hc_phi", "#phi [rad]", 20, -3.5, 3.5);
   TH2F *hc_thetaphi = new TH2F("hc_thetaphi", "", 100, -1, 1, 100, -3.5, 3.5);

   //Angular distribution of events passing the cuts for origin = neutrino
   TH1F *hn_theta = new TH1F("hn_theta", "cos(#theta)", 20, -1, 1);
   TH1F *hn_phi = new TH1F("hn_phi", "#phi [rad]", 20, -3.5, 3.5);
   TH2F *hn_thetaphi = new TH2F("hn_thetaphi", "", 100, -1, 1, 100, -3.5, 3.5);

   //Angular distribution of events passing the cuts
   TH1F *hnc_theta = new TH1F("hnc_theta", "cos(#theta)", 20, -1, 1);
   TH1F *hnc_phi = new TH1F("hnc_phi", "#phi [rad]", 20, -3.5, 3.5);
   TH2F *hnc_thetaphi = new TH2F("hnc_thetaphi", "", 100, -1, 1, 100, -3.5, 3.5);
   TH1F *hnc_origin = new TH1F("hnc_origin", "", 6, -1.5, 4.5);

   //Distance between true and reco vertex in R, x, y, and z (in cm)
   TH1F *hdisttruemin = new TH1F("disttruemin", "", 100, 0, 1000);
   TH1F *hdisttrueminx = new TH1F("disttrueminx", "", 1000, -50, 50);
   TH1F *hdisttrueminy = new TH1F("disttrueminy", "", 1000, -50, 50);
   TH1F *hdisttrueminz = new TH1F("disttrueminz", "", 1000, -50, 50);

   //Distance between true and reco vertex
   TH1F *hdisttruemin1 = new TH1F("disttruemin1", "", 100, 0, 1000); //for single track events
   TH1F *hdisttruemin2 = new TH1F("disttruemin2", "", 100, 0, 1000); //for two track events
   TH1F *hdisttruemin2cut = new TH1F("disttruemin2cut", "", 100, 0, 1000); //for two track events with tracks not back-to-back
   
   // Vertex position calues
   TH1F *hXVertexPosition = new TH1F("X Vertex Distribution","X Vertex Distribution",256,0,FVx);
   hXVertexPosition->GetXaxis()->SetTitle("X position [cm]");
   hXVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
   TLine XVertexMinCut(borderx, 0, borderx, 200);
   TLine XVertexMaxCut(FVx-borderx, 0, FVx-borderx, 200);
   XVertexMinCut.SetLineColor(kRed);
   XVertexMaxCut.SetLineColor(kRed);
   
   TH1F *hYVertexPosition = new TH1F("Y Vertex Distribution","Y Vertex Distribution",333,-FVy/2,(FVy+100)/2);
   hYVertexPosition->GetXaxis()->SetTitle("Y position [cm]");
   hYVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
   TLine YVertexMinCut(-FVy/2+bordery, 0, -FVy/2+bordery, 200);
   TLine YVertexMaxCut(FVy/2-bordery, 0, FVy/2-bordery, 200);
   YVertexMinCut.SetLineColor(kRed);
   YVertexMaxCut.SetLineColor(kRed);
   
   TH1F *hZVertexPosition = new TH1F("Z Vertex Distribution","Z Vertex Distribution",1000,0,FVz);
   hZVertexPosition->GetXaxis()->SetTitle("Z position [cm]");
   hZVertexPosition->GetYaxis()->SetTitle("Number of Vertices [ ]");
   TLine ZVertexMinCut(borderz, 0, borderz, 200);
   TLine ZVertexMaxCut(FVz-borderz, 0, FVz-borderz, 200);
   ZVertexMinCut.SetLineColor(kRed);
   ZVertexMaxCut.SetLineColor(kRed);
   
   // Flash Time distribution
   TH1F *hFlashTime = new TH1F("Flash Time","Flash Time",200,0,8.2);
   hFlashTime->GetXaxis()->SetTitle("Time [#mus]");
   hFlashTime->GetYaxis()->SetTitle("Number of Flashes [ ]");
   TLine FlashTimeMinCut(beammin, 0, beammin, 2000);
   TLine FlashTimeMaxCut(beammax, 0, beammax, 2000);
   FlashTimeMinCut.SetLineColor(kRed);
   FlashTimeMaxCut.SetLineColor(kRed);
   
   // Flash PE distirbution
   TH1F *hIntesityPE = new TH1F("Flash Intensity in PE","Flash Intensity in PE",100,0,1000);
   hIntesityPE->GetXaxis()->SetTitle("Number of P.E. [ ]");
   hIntesityPE->GetYaxis()->SetTitle("Number of Flashes [ ]");
   TLine FlashPEMinCut(PEthresh, 0, PEthresh, 20000);
   FlashPEMinCut.SetLineColor(kRed);
   
   // Flash to track distance
   TH1F *hFlashTrackDist = new TH1F("Distance from Flash to Track","Distance between Flash and Track",500,0,500);
   hFlashTrackDist->GetXaxis()->SetTitle("Distance [cm]");
   hFlashTrackDist->GetYaxis()->SetTitle("Number of Events [ ]");
   TLine FlashTrackMinCut(flashwidth, 0, flashwidth, 3000);
   FlashTrackMinCut.SetLineColor(kRed);
   
   // Vertex to Track distance
   TH1F *hVertexTrackDist = new TH1F("Distance from Vertex to Track","Distance between Vertex and Track",200,0,100);
   hVertexTrackDist->GetXaxis()->SetTitle("Distance [cm]");
   hVertexTrackDist->GetYaxis()->SetTitle("Number of Tracks [ ]");
   TLine VertexTrackMinCut(distcut, 0, distcut, 3000);
   VertexTrackMinCut.SetLineColor(kRed);
   
   // TrackLength of longest tracks
   TH1I *hTrackLength = new TH1I("Track Length Distribution","Track Length Distribution",1000,0,1000);
   hTrackLength->GetXaxis()->SetTitle("Track Length [cm]");
   hTrackLength->GetYaxis()->SetTitle("Number of Tracks [ ]");
   
   // Track Multiplicity at Vertex
   TH1F *hTrackMultip = new TH1F("Track Multiplicity","Track Multiplicity at the Vertex",10,0,10);
   hTrackMultip->GetXaxis()->SetTitle("Number of Tracks [ ]");
   hTrackMultip->GetYaxis()->SetTitle("Number of Events [ ]");
   
   // Track Multiplicity at Vertex
   TH1F *hFlashVertexDist = new TH1F("Flash To Vertex Distance","Flash To Vertex Distance",1000,0,1000);
   hFlashVertexDist->GetXaxis()->SetTitle("Distance [cm]");
   hFlashVertexDist->GetYaxis()->SetTitle("Number of Vertices [ ]");
   
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
   
   
   //for two track events, what is the angle between the tracks?
   TH1F *htwotrackangle = new TH1F("twotrackangle", "", 20, -1, 1);

   double ds = 0;
   double de = 0;
   double dx1 = 0;
   double dx2 = 0;
   double dy1 = 0;
   double dy2 = 0;
   double dz1 = 0;
   double dz2 = 0;
   double d1 = 0;
   double d2 = 0;
   double angle = 0;
   int ntracks5 = 0;
   int ntracks10 = 0;

   int Size = treenc -> GetEntries();
   if(Size > 20000) Size = 20000;
   cout << "number of events used is: " << Size << endl;

   bool flashtag = false;
   int theflash = -1;
   double flashmax = 0;

   double sumpot = 0;
   double R = 0;

   double diststart = 0;
   double distend = 0;
   double length = 0;

   int ntrue = 0;
   int nnc = 0;
   int nnue = 0;
   int nnumubar = 0;
   double disttrue[maxmc];
   double disttruemin = 0;
   double disttrueminx = 0;
   double disttrueminy = 0;
   double disttrueminz = 0;
   
   unsigned int EventsWithFlash = 0;
   unsigned int EventsVtxInFV = 0;
   unsigned int EventsTrackNearVertex = 0;
   unsigned int EventsFlashMatched = 0;
   unsigned int EventsTracksInFV = 0;
   unsigned int EventsLengthCut = 0;
   unsigned int EventsNearVtx = 0;
   unsigned int EventsTrackLong = 0;
   
   unsigned int NumberOfContainedMCTracks = 0;
   
//    vector<float> TrackCandidates;

   //Event Loop
   for(int i = 0; i < Size; i++) 
   {
      if(i%1000 == 0) cout << "\t... " << i << endl;

      treenc -> GetEntry(i);
      
      bool ContainedMCTrackFlag = true;
      
      // Loop over all MC particles
      for(int track_no = 0; track_no < NumberOfMCTracks; track_no++)
      {
	// If the a muon is not contained in a singel neutrino event, set mc-track contained flag to false
	if( abs(PDG_truth[track_no])==13 && ContainedMCTrackFlag && mcevts_truth == 1
	  && (!inFV(XMCTrackStart[track_no],YMCTrackStart[track_no],ZMCTrackStart[track_no]) || !inFV(XMCTrackEnd[track_no],YMCTrackEnd[track_no],ZMCTrackEnd[track_no])) )
	{
	  ContainedMCTrackFlag = false;
	}
      } // MC particle loop
      
      // Count up the number of contained mc-tracks if the flag is true
      if(ContainedMCTrackFlag)
      {
	NumberOfContainedMCTracks++;
      }
      
      flashtag = false;
      flashmax = 0;
      
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
	
	// Increase events with flash > 50 PE and within beam window
	EventsWithFlash++;
	
	// Loop over all reconstructed tracks
	for(int j = 0; j < ntracks_reco; j++) 
	{
	  // Fill track length histogram and the flash track histogram
	  hTrackLength->Fill( pow(trkstartx[j] - trkendx[j],2) + pow(trkstarty[j] - trkendy[j],2) + pow(trkstartz[j] - trkendz[j],2) );  
	  hFlashTrackDist->Fill(FlashTrackDist(flash_zcenter[theflash], trkstartz[j], trkendz[j]));
	} // reco track loop
	
	// Loop over all vertices
	for(int v = 0; v < nvtx; v++) 
	{
	  // Fill vertex position to histogram
	  hXVertexPosition->Fill(vtxx[v]);
	  hYVertexPosition->Fill(vtxy[v]);
	  hZVertexPosition->Fill(vtxz[v]);

	  float TrackCandLength = 0;
	  int TrackCandidate = -1;
	  
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
	
	      // Figure out if the end or the start of the track was closer to the vertex
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
	    
	    // If the longest track length is filled
	    if(TrackCandLength)
	    {
	      // If the longest track is flash matched
	      if( FlashTrackDist(flash_zcenter[theflash], trkstartz[TrackCandidate], trkendz[TrackCandidate]) < flashwidth ) 
	      {
		EventsFlashMatched++;
		
		// If the longest track is fully contained
		if( inFV(trkstartx[TrackCandidate], trkstarty[TrackCandidate], trkstartz[TrackCandidate]) && inFV(trkendx[TrackCandidate], trkendy[TrackCandidate], trkendz[TrackCandidate]) )
		{
		  EventsTracksInFV++;
		  
		  // If longest track is longer than 75 cm
		  if(TrackCandLength > lengthcut)
		  {
		    EventsTrackLong++;
		  }
		}
	      }
	    } // If there is a longest track

	  } // if vertex is contained
	} // vertex loop
      } // if flashtag
      
      ntracks5 = 0;
      ntracks10 = 0;

//       if(mcevts_truth == 1) { //to make it easy: require only 1 neutrino interaction per event
//          if(inFV(nuvtxx_truth[0], nuvtxy_truth[0], nuvtxz_truth[0])) { //require that the true vertex is inside the fiducial volume
//             if(nuPDG_truth[0] == 14) { //require that it was a numu primary
//                if(ccnc_truth[0] == 0 && mode_truth[0] == 0) { //require CC QE interaction
                  ntrue++; //cout number of neutrinos that fulfill this requirement

                  diststart = 0;
                  distend = 0;

                  numuvertex = -1;
                  mutrack = -1;
                  mutracklength = 0;

                  if(true) {
                     if(flashtag == true) { //only for events with a flash in the beam window
		       
		       
		       
		       bool TracksInFVFlag = false;
		       bool FlashMatchedFlag = false;
		       
		       
                        for(int v = 0; v < nvtx; v++) { //loop over all reconstructed vertices
                           candvertex[v] = true;
                           candtrack[v] = -1;
                           candlength[v] = 0;
                           if(inFV(vtxx[v], vtxy[v], vtxz[v])) { //check if reco vertex is in fiducial volume
			     if(ntracks_reco)
			     {
			       FlashMatchedFlag = true;
			       TracksInFVFlag = true;
			     }
			     for(int j = 0; j < ntracks_reco; j++) { //loop over all tracks

                                 //calculate distance between track start point and vertex, and track end point and vertex
                                 diststart = sqrt((vtxx[v] - trkstartx[j])*(vtxx[v] - trkstartx[j]) + (vtxy[v] - trkstarty[j])*(vtxy[v] - trkstarty[j]) + (vtxz[v] - trkstartz[j])*(vtxz[v] - trkstartz[j]));
                                 distend = sqrt((vtxx[v] - trkendx[j])*(vtxx[v] - trkendx[j]) + (vtxy[v] - trkendy[j])*(vtxy[v] - trkendy[j]) + (vtxz[v] - trkendz[j])*(vtxz[v] - trkendz[j]));

                                 //find the end which is closest to the vertex, and define this as start of the track to get the track direction
                                 if(diststart < distend) {
                                    vertexatstart[j] = true;
                                    vertexatend[j] = false;
                                 }
                                 else {
                                    vertexatstart[j] = false;
                                    vertexatend[j] = true;
                                 }

                                 if(diststart < distcut || distend < distcut) { //check if the tracks start within distance x of vertex
				   
                                    if(FlashTrackDist(flash_zcenter[theflash], trkstartz[j], trkendz[j]) < flashwidth) { //check if the tracks fulfill flash matching
                                       if(inFV(trkstartx[j], trkstarty[j], trkstartz[j]) && inFV(trkendx[j], trkendy[j], trkendz[j])) { //check if the track is fully contained
                                          length = sqrt(pow(trkstartx[j] - trkendx[j],2) + pow(trkstarty[j] - trkendy[j],2) + pow(trkstartz[j] - trkendz[j],2));
                                          if(length > candlength[v]) { //find the longest track from this vertex
                                             candtrack[v] = j;
                                             candlength[v] = length;
                                          }
                                       }
                                       else {
                                          candvertex[v] = false; //if 1 or more tracks from this vertex are uncontained, it fails
                                          TracksInFVFlag = false;
                                       }
                                    }
                                    else {
                                       candvertex[v] = false; //if 1 or more tracks from this vertex fail flash matching, it fails
                                       FlashMatchedFlag = false;
				       TracksInFVFlag = false;
                                    }
                                 }
                                 else
				 {
				   TracksInFVFlag = false;
				   FlashMatchedFlag = false;
				 }
                              }
                           }
                           else 
			   {
			     candvertex[v] = false; //if vertex is not in fiducial volume
			   }
                           if(candvertex[v] == false) {
                              candtrack[v] = -1; //no candidate muon track
                              candlength[v] = -100; //no candiate muon track length
                           }
                           if(candlength[v] > mutracklength) {
                              mutracklength = candlength[v];
                              mutrack = candtrack[v];
                              numuvertex = v;
                           }
                           if(FlashMatchedFlag) // increase Flash matching in position count
			   {
			     FlashMatchedFlag = false;
			   }
			   
                        }
		     }
                  }
                  if(mutracklength > lengthcut) { //check if candidate muon passes length cut
                     if(vertexatstart[mutrack]) {
                        hnc_theta -> Fill(cos(trktheta[mutrack]));
                        hnc_thetaphi -> Fill(cos(trktheta[mutrack]), trkphi[mutrack]);
                        hnc_phi -> Fill(trkphi[mutrack]);
                        hnc_origin -> Fill(trkorigin[mutrack][2]);
                     }
                     else {
                        hnc_theta -> Fill(-cos(trktheta[mutrack]));
                        hnc_thetaphi -> Fill(-cos(trktheta[mutrack]), trkphi[mutrack]);
                        hnc_phi -> Fill(trkphi[mutrack]);
                        hnc_origin -> Fill(trkorigin[mutrack][2]);
                     }
                     disttruemin = 1000000;
                     for(int m = 0; m < mcevts_truth; m++) {
                        //distmin with vertex instead of track
                        disttrue[m] = sqrt((vtxx[numuvertex] - nuvtxx_truth[m])*(vtxx[numuvertex] - nuvtxx_truth[m]) + (vtxy[numuvertex] - nuvtxy_truth[m])*(vtxy[numuvertex] - nuvtxy_truth[m]) + (vtxz[numuvertex] - nuvtxz_truth[m])*(vtxz[numuvertex] - nuvtxz_truth[m]));
                        if(disttrue[m] < disttruemin) {
                           disttruemin = disttrue[m];
                           disttrueminx = vtxx[numuvertex] - nuvtxx_truth[m];
                           disttrueminy = vtxy[numuvertex] - nuvtxy_truth[m];
                           disttrueminz = vtxz[numuvertex] - nuvtxz_truth[m];
                        }
                     }
                     hdisttruemin -> Fill(disttruemin);
                     hdisttrueminx -> Fill(disttrueminx);
                     hdisttrueminy -> Fill(disttrueminy);
                     hdisttrueminz -> Fill(disttrueminz);

                     //count how many tracks are starting within 10 cm/5cm of vertex 
                     for(int k = 0; k < ntracks_reco; k++) {
                        ds = sqrt((trkstartx[k] - vtxx[numuvertex])*(trkstartx[k] - vtxx[numuvertex]) + (trkstarty[k] - vtxy[numuvertex])*(trkstarty[k] - vtxy[numuvertex]) + (trkstartz[k] - vtxz[numuvertex])*(trkstartz[k] - vtxz[numuvertex]));
                        de = sqrt((trkendx[k] - vtxx[numuvertex])*(trkendx[k] - vtxx[numuvertex]) + (trkendy[k] - vtxy[numuvertex])*(trkendy[k] - vtxy[numuvertex]) + (trkendz[k] - vtxz[numuvertex])*(trkendz[k] - vtxz[numuvertex]));
                        if(ds < 5 || de < 5) ntracks5++;
                        if(ds < 10 || de < 10) ntracks10++;
                     }
                     //make distributions for 1-track and 2-track events
                     if(ntracks10 == 1) hdisttruemin1 -> Fill(disttruemin);
                     else if(ntracks10 == 2) hdisttruemin2 -> Fill(disttruemin);
                     if(ntracks10 == 2) {
                        dx1 = 0;
                        dx2 = 0;
                        dy1 = 0;
                        dy2 = 0;
                        dz1 = 0;
                        dz2 = 0;
                        d1 = 0;
                        d2 = 0;
                        for(int k = 0; k < ntracks_reco; k++) {
                           ds = sqrt((trkstartx[k] - vtxx[numuvertex])*(trkstartx[k] - vtxx[numuvertex]) + (trkstarty[k] - vtxy[numuvertex])*(trkstarty[k] - vtxy[numuvertex]) + (trkstartz[k] - vtxz[numuvertex])*(trkstartz[k] - vtxz[numuvertex]));
                           de = sqrt((trkendx[k] - vtxx[numuvertex])*(trkendx[k] - vtxx[numuvertex]) + (trkendy[k] - vtxy[numuvertex])*(trkendy[k] - vtxy[numuvertex]) + (trkendz[k] - vtxz[numuvertex])*(trkendz[k] - vtxz[numuvertex]));
                           if((ds < 10 || de < 10) && d1 == 0) {
                              dx1 = trkstartx[k] - trkendx[k];
                              dy1 = trkstarty[k] - trkendy[k];
                              dz1 = trkstartz[k] - trkendz[k];
                              d1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
                           }
                           else if((ds < 10 || de < 10) && d1 > 0) {
                              dx2 = trkstartx[k] - trkendx[k];
                              dy2 = trkstarty[k] - trkendy[k];
                              dz2 = trkstartz[k] - trkendz[k];
                              d2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                           }
                        }
                        angle = (dx1*dx2+dy1*dy2+dz1*dz2)/(d1*d2);
                        htwotrackangle -> Fill(angle);
                        if(fabs(angle) < 0.95) hdisttruemin2cut -> Fill(disttruemin);
                     } 
                  }
//                }//interaction is CC QE
//             }//particle is numu
//          }//true vertex in FV
//       }//mcevts == 1
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
   
   TCanvas *Canvas10 = new TCanvas("Track Multiplicity", "Track Multiplicity", 1400, 1000);
   hFlashVertexDist->Draw();
   
   TFile* OutputFile = new TFile(("Hist_"+ProductName+ "_"+Version+".root").c_str(),"RECREATE");
   
   Canvas1->Write();
   Canvas2->Write();
   Canvas3->Write();
   Canvas4->Write();
   Canvas5->Write();
   Canvas6->Write();
   Canvas7->Write();
   Canvas8->Write();
   Canvas10->Write();
   
   hXVertexPosition->SetStats(0);
   hYVertexPosition->SetStats(0);
   hZVertexPosition->SetStats(0);
   hFlashTime->SetStats(0);
   hIntesityPE->SetStats(0);
   hFlashTrackDist->SetStats(0);
   hVertexTrackDist->SetStats(0);
   hTrackLength->SetStats(0);
   hTrackMultip->SetStats(0);
   hFlashVertexDist->SetStats(0);
   
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
   
   OutputFile->Close();
   
   Canvas10->SaveAs(("FlashVertexDist"+Version+".png").c_str());
//    Canvas2->SaveAs((ProductName+"_YVtxPosition_"+Version+".png").c_str());
//    Canvas3->SaveAs((ProductName+"_ZVtxPosition_"+Version+".png").c_str());
//    Canvas4->SaveAs((ProductName+"_FlashTime_"+Version+".png").c_str());
//    Canvas5->SaveAs((ProductName+"_PECount_"+Version+".png").c_str());
//    Canvas6->SaveAs((ProductName+"_FlashTrackDist_"+Version+".png").c_str());
//    Canvas7->SaveAs((ProductName+"_VtxTrackDist_"+Version+".png").c_str());
//    Canvas8->SaveAs((ProductName+"_TrackLength_"+Version+".png").c_str());

//    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
//    hc_theta -> SetLineWidth(2);
//    hc_theta -> SetLineColor(4);
//    hc_theta -> Draw();
//    hn_theta -> SetLineWidth(2);
//    hn_theta -> SetLineColor(8);
//    hn_theta -> Draw("same");
//    hnc_theta -> SetLineWidth(2);
//    hnc_theta -> SetLineColor(2);
//    hnc_theta -> Draw("same");
// //    c1 -> SaveAs("c1_pandoraNuKHit_MC.root");
// 
//    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
//    hc_phi -> SetLineWidth(2);
//    hc_phi -> SetLineColor(4);
//    hc_phi -> Draw();
//    hn_phi -> SetLineWidth(2);
//    hn_phi -> SetLineColor(8);
//    hn_phi -> Draw("same");
//    hnc_phi -> SetLineWidth(2);
//    hnc_phi -> SetLineColor(2);
//    hnc_phi -> Draw("same");
// //    c2 -> SaveAs("c2_pandoraNuKHit_MC.root");
// 
//    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
//    hc_thetaphi -> SetMarkerStyle(20);
//    hc_thetaphi -> SetMarkerColor(4);
//    hc_thetaphi -> Draw("SCAT");
//    hn_thetaphi -> SetMarkerStyle(20);
//    hn_thetaphi -> SetMarkerColor(8);
//    hn_thetaphi -> Draw("SCATSAME");
//    hnc_thetaphi -> SetMarkerStyle(20);
//    hnc_thetaphi -> SetMarkerColor(2);
//    hnc_thetaphi -> Draw("SCATSAME");
// //    c3 -> SaveAs("c3_pandoraNuKHit_MC.root");
// 
//    TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
//    hnc_origin -> SetLineWidth(2);
//    hnc_origin -> Draw();
// //    c4 -> SaveAs("c4_pandoraNuKHit_MC.root");
// 
//    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
//    hdisttruemin -> SetLineWidth(2);
//    hdisttruemin -> Draw();
//    hdisttruemin1 -> SetLineWidth(2);
//    hdisttruemin1 -> SetLineColor(8);
//    hdisttruemin1 -> Draw("same");
//    hdisttruemin2 -> SetLineWidth(2);
//    hdisttruemin2 -> SetLineColor(2);
//    hdisttruemin2 -> Draw("same");
//    hdisttruemin2cut -> SetLineWidth(2);
//    hdisttruemin2cut -> SetLineColor(kOrange+1);
//    hdisttruemin2cut -> Draw("same");
// //    c5 -> SaveAs("c5_pandoraNuKHit_MC.root");
// 
//    TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);
//    hdisttrueminx -> SetLineWidth(2);
//    hdisttrueminx -> Draw();
// //    c6 -> SaveAs("c6_pandoraNuKHit_MC.root");
// 
//    TCanvas *c7 = new TCanvas("c7", "c7", 600, 600);
//    hdisttrueminy -> SetLineWidth(2);
//    hdisttrueminy -> Draw();
// //    c7 -> SaveAs("c7_pandoraNuKHit_MC.root");
// 
//    TCanvas *c8 = new TCanvas("c8", "c8", 600, 600);
//    hdisttrueminz -> SetLineWidth(2);
//    hdisttrueminz -> Draw();
// //    c8 -> SaveAs("c8_pandoraNuKHit_MC.root");
// 
//    TCanvas *c11 = new TCanvas("c11", "c11", 600, 600);
//    htwotrackangle -> SetLineWidth(2);
//    htwotrackangle -> Draw();
   
   cout << "number of CC events with vertex in FV : " << ntrue << endl;
   cout << "number of events with flash > 50 PE : " << EventsWithFlash << endl;
   cout << "number of events with vtx in FV : " << EventsVtxInFV << endl;
   cout << "number of events with track start/end within 5cm to vtx : " << EventsTrackNearVertex << endl;
   cout << "number of events with tracks matched within 80cm to flash : " << EventsFlashMatched << endl;
   cout << "number of events with contained tracks : " << EventsTracksInFV << endl;
   cout << "number of events with longest track > 75cm : " << EventsTrackLong << endl;
   cout << "number of events with contained MC tracks : " << NumberOfContainedMCTracks << endl;
//    cout << "number of selected events: " << hnc_theta -> Integral() << endl;
//    cout << "number of neutrino origin events (backtracker): " << hn_theta -> Integral() << endl;
//    cout << "number of cosmic origin events (backtracker): " << hc_theta -> Integral() << endl;
   cout << "nue: " << nnue << endl;
   cout << "NC: " << nnc << endl;
   cout << "nnumubar: " << nnumubar << endl;



   return 0;

}
