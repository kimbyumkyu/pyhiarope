#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TRandom.h>
#include <TClonesArray.h>
#include <Pythia8/Pythia.h>
#include <TStopwatch.h>
#include <TAxis.h>
#include <THnSparse.h>
#include <THashList.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector2.h>
// Configuration for Toy MC track generation //
const double pi = TMath::Pi();

using namespace std;
using namespace Pythia8;
TFile* fOutput;

THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt) {
  //TString dirname(basename(name)), hname(histname(name));
  //THashList *parent(FindGroup(dirname));
  //if(!parent) parent = CreateHistoGroup(dirname);
  //if(parent->FindObject(hname)){
    //Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
    //return 0;
  //}
  TArrayD xmin(ndim), xmax(ndim);
  TArrayI nbins(ndim);
  for(int idim = 0; idim < ndim; ++idim){
    const TAxis &myaxis = *(axes[idim]);
    nbins[idim] = myaxis.GetNbins();
    xmin[idim] = myaxis.GetXmin();
    xmax[idim] = myaxis.GetXmax();
  }
  THnSparseD *hsparse = new THnSparseD(name, title, ndim, nbins.GetArray(), xmin.GetArray(), xmax.GetArray());
  for(int id = 0; id < ndim; ++id)
    *(hsparse->GetAxis(id)) = *(axes[id]);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    hsparse->Sumw2();
  //parent->Add(hsparse);
  return hsparse;
}



THnSparse * CreateTHnSparse(TString name
    , TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  auto h= CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}

THnSparse * CreateTHnSparse(TString name
    , TString title, TString templ, Option_t * opt){
  auto o = fOutput->FindObject(templ);
  if( !o ) {
    cout<<"ERROR: no "<<templ<<endl;
    exit(1);
  }
  auto ht = dynamic_cast<THnSparse*>( o );
  const TAxis * axises[ht->GetNdimensions()];
  for( int i=0;i<ht->GetNdimensions();i++ ) axises[i]= ht->GetAxis(i);
  auto h= CreateTHnSparse(name, title, ht->GetNdimensions(), axises,opt );
  return h;
}
Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}

Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fOutput->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}

TAxis AxisFix
( TString name, int nbin, Double_t xmin, Double_t xmax ){
  TAxis axis(nbin, xmin, xmax);axis.SetName(name);
  return axis;
}

TAxis AxisStr( TString name, std::vector<TString> bin ){
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax;
}

TAxis AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}

TAxis AxisLog
( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}




int main(int argc, char **argv) {

	//==== Read arguments =====
	if ( argc<4 ) {
		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"+  "<<argv[0]<<" [outputFile] [Seed] [Nevt]  "<<endl;
		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout << endl << endl;
		exit(1);
	}

	char *outFile;
	outFile   = argv[1];
	int randomseed = atoi(argv[2]);
	int Nevt = atoi(argv[3]);

	//char* pythiaconfig  = argv[4];

	//---------------------
	//Pythia initialization
	//---------------------
	Pythia pythia;   // Generator.
	//Event& event      = pythia.event;
	//ParticleData& pdt = pythia.particleData;
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("SoftQCD:nonDiffractive = on");
	pythia.readString("Next:numberShowEvent = 0");
	// Enabling string shoving, setting model parameters.
	// The model is still untuned. These parameter values
	// are choosen for illustrative purposes.
	pythia.readString("Ropewalk:RopeHadronization = on");
	pythia.readString("Ropewalk:doShoving = on");
	pythia.readString("Ropewalk:doFlavour = off");
	pythia.readString("Ropewalk:rCutOff = 10.0");
	pythia.readString("Ropewalk:limitMom = on");
	pythia.readString("Ropewalk:pTcut = 2.0");
	pythia.readString("Ropewalk:r0 = 0.41");
	pythia.readString("Ropewalk:m0 = 0.2");
	pythia.readString("Ropewalk:gAmplitude = 10.0");
	pythia.readString("Ropewalk:gExponent = 1.0");
	pythia.readString("Ropewalk:deltat = 0.1");
	pythia.readString("Ropewalk:tShove = 1.");
	pythia.readString("Ropewalk:deltay = 0.1");
	pythia.readString("Ropewalk:tInit = 1.5");
	// Enabling setting of vertex information.
	pythia.readString("PartonVertex:setVertex = on");
	pythia.readString("PartonVertex:protonRadius = 0.7");
	pythia.readString("PartonVertex:emissionWidth = 0.1");
	pythia.init();




	// Read in commands from external file.
	//pythia.readFile(pythiaconfig);

	// Extract settings to be used in the main program.

	pythia.readString("Random:setSeed = on");
	pythia.readString(Form("Random:seed=%d",randomseed));
	// Initialize. Beam parameters set in .cmnd file.
	pythia.init();


	TStopwatch timer;
	timer.Start();

	typedef std::vector<Double_t> Double1D;
	TH1D *hMultPythia = new TH1D("hMultiPythia","Nch via Pythia",2001,-1,2000);
  Double1D varcentbin = {0,0.001,0.0033,0.01,0.02,0.033,0.05,0.1,0.2,0.5,1,2,5,10,15,20,30,40,50,70,82,100};
	Double1D varmultbin = {0, 35, 80, 105, 1000};
	auto binCent = AxisVar("Cent", varmultbin);
	auto binPhi = AxisFix("phi",32,-0.5*pi-pi/32.0, 1.5*pi-pi/32.0);
	auto binEta = AxisFix("eta",80,-4.0,4.0);
	Double1D ptbin = {0.1,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,14.0,100};
	auto binTPt = AxisVar("TPt",ptbin); //trig
	auto binAPt = AxisVar("APt",ptbin); //associatei
	Double1D ltpttrackbin = {0.2, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 13.0, 20.0};
	auto binLtpt = AxisVar("LPPt",ltpttrackbin);
	auto hRidgeLT = CreateTHnSparse("hRidgeLT","RidgeLT",6,{binCent,binPhi,binEta,binTPt,binAPt,binLtpt},"s");
  auto hRidgeSLT = CreateTHnSparse("hRidgeMixingSLT","RidgeMixingSLT",6,{binCent,binPhi,binEta,binTPt,binAPt,binLtpt},"s");

	auto binNtrig = AxisFix("Ntrig",1,0.5,1.5);
	auto hNtrig = CreateTHnSparse("hNtrig","hNtrig",4,{binCent,binTPt,binNtrig,binLtpt},"s"); 
	//double etaMaxCutForPart = 5.0;
	double fCent = 0;


	typedef std::vector<TLorentzVector*> TRACKS;
	typedef std::deque<TRACKS> EVENTPOOL ;
	typedef std::vector<EVENTPOOL> MIXINGPOOL;


	TRACKS *tracks;

	MIXINGPOOL mixingpool;
	mixingpool.resize(binCent.GetNbins());	


	EVENTPOOL *eventpool;
	UInt_t centbin=0;	
	for(int ievt=0; ievt<Nevt; ievt++){
		if (!pythia.next()) continue;
		if(ievt % 1000 == 0) cout << ievt << endl ;
		// Particle generation ===========================================
		// check the multiploicity
		double mult = 0;
		double MaxPt = 0;
		

	
		for (int i = 0; i < pythia.event.size(); ++i) {//loop over all  particles in the event
			Particle& p = pythia.event[i];
			// Apply simple, particle level, cuts.
			if (!p.isFinal() || !p.isCharged() || abs(p.eta()) > 2.5 || p.pT() < 0.5) continue;
			++mult;
		}
		fCent = mult;
		centbin = binCent.FindBin(fCent) -1;
		eventpool = &mixingpool.at(centbin);
		eventpool->push_back(TRACKS());
		tracks = &(eventpool->back());

	
		for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
			Particle& p = pythia.event[i];
			// Apply simple, particle level, cuts.
			if (!p.isFinal() || !p.isCharged() || abs(p.eta()) > 2.5 ) continue;
			if (p.pT()>MaxPt) MaxPt = p.pT();
			tracks->push_back( new TLorentzVector(p.px(),p.py(),p.pz(),p.e()) );
		}

		
		if (mult==0) eventpool->pop_back();
		if ( eventpool->size() > 11 ){
			for (auto it: eventpool->front()) delete it;
			eventpool->pop_front();
		}

		TRACKS alltracks;
		int nn =0;
		for (auto pool: *eventpool){
			if (nn >10) continue;
			for (auto trk: pool) alltracks.push_back(trk);
			nn++;
		}

	
	
		for (int i = 0; i < pythia.event.size() -1 ; ++i) {//loop over particles
			Particle& p = pythia.event[i];
			if(!p.isFinal() || !p.isCharged() || abs(p.eta()) > 2.5 ) continue;
			TLorentzVector P (p.px(), p.py(), p.pz(), p.e());
			for (int j = i+1; j < pythia.event.size()  ; ++j) {//loop over particles
				Particle& q = pythia.event[j];
				if(!q.isFinal() || !q.isCharged() || abs(q.eta()) > 2.5 ) continue;
				TLorentzVector Q (q.px(), q.py(), q.pz(), q.e());
				double deltaeta = P.Eta() - Q.Eta();
				double deltaphi = P.Phi() - Q.Phi();
				if( P.Pt() < Q.Pt()){ deltaeta *= -1.0; deltaphi *= -1.0;}
				deltaphi = TVector2::Phi_0_2pi(deltaphi);
				if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;
				FillTHnSparse(hRidgeLT,{fCent, deltaphi, deltaeta, max(P.Pt(),Q.Pt()), min(P.Pt(),Q.Pt()),MaxPt}, 1.0 );
				FillTHnSparse(hNtrig,{fCent, max(P.Pt(),Q.Pt()),1.0,MaxPt},1.0);
			}

		}
		
		if (eventpool->size() == 11) {
			for (int i = 0; i < pythia.event.size()  ; ++i) {//loop over particles
				Particle& p = pythia.event[i];
				if(!p.isFinal() || !p.isCharged() || abs(p.eta()) > 2.5 ) continue;
				TLorentzVector P (p.px(), p.py(), p.pz(), p.e());
				for (auto Q :  alltracks ) {//loop over particles
					double deltaeta = P.Eta() - Q->Eta();
					double deltaphi = P.Phi() - Q->Phi();
					if( P.Pt() < Q->Pt()){ deltaeta *= -1.0; deltaphi *= -1.0;}
					deltaphi = TVector2::Phi_0_2pi(deltaphi);
					if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;
					FillTHnSparse(hRidgeSLT,{fCent, deltaphi, deltaeta, max(P.Pt(),Q->Pt()), min(P.Pt(),Q->Pt()),MaxPt}, 1.0 );
				}

			}
		}
		
		hMultPythia -> Fill (mult);
		//EventCounter++;

	} // event loop done.


	fOutput = new TFile(outFile ,"recreate" );

	hMultPythia->Write();
	hRidgeLT->Write();
	hRidgeSLT->Write();
	hNtrig->Write();
	fOutput->Write();
	fOutput->Close();
	timer.Print();
}



