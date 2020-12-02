#define canvassizex 1250
#define canvassizey 800

//_____<GENERIC FUNCTIONS

//Return min value of a vector
double minTime(vector<double> v) {
  double tmin =1e200;;
  for(int i=0; i<v.size(); i++)
    if(v.at(i) < tmin)
      tmin = v.at(i);
  return tmin;
}

//Return max value of a vector
double maxTime(vector<double> v) {
  double tmax = -1e200;;
  for(int i=0; i<v.size(); i++)
    if(v.at(i) > tmax)
      tmax = v.at(i);
  return tmax;
}

//Return min value of a vector of vectors
double minTime2(vector<vector<double>> * v) {
  double tmin=1e200;
  for(int j=0; j<v->size(); j++) {
	  for(int i=0; i<v->at(j).size(); i++)
	    if(v->at(j).at(i) < tmin)
	      tmin = v->at(j).at(i);
  }
  return tmin;
}

//Return max value of a vector of vectors
double maxTime2(vector<vector<double>> * v) {
  double tmax=-1e200;
  int OC=0;
  for(int j=0; j<v->size(); j++) {
	  for(int i=0; i<v->at(j).size(); i++)
	    if(v->at(j).at(i) > tmax){
	      tmax= v->at(j).at(i);
		  OC=j;
	  }
  }
  return tmax;
}

//Return max value between x and y
double fmax(double x, double y) {
  if(y <= x)
    return x;
  else
    return y;
}

//_____>

//---<FUNCTION TO PRINT PARTICLE INFO
//Structure to store event info
struct structInfoTrack {
  int Run, Ev, Tr, PDG, Mother;
	double E, dE;
	string Process;
};
void print_trackInfo(struct structInfoTrack ITr, int file, int iter){
	cout<<"(F, i)=("<<setw(2)<<file<<","<<setw(3)<<iter<<") EvID:"<<setw(3)<<ITr.Ev<<" Energy:"<<setw(10)<<ITr.E;
	cout<<" TrackID:"<<setw(4)<<ITr.Tr<<" MotherID:"<<setw(4)<<ITr.Mother;
	cout <<" PDG:"<<setw(12)<<ITr.PDG<<" dE/dx: "<<setw(8)<<ITr.dE<<" Process: "<<ITr.Process<<endl;
}
//--->

class PDSEventHandle{
	public:
		PDSEventHandle(std::string geofile, int run, int subrun, int event, vector<double>* E, vector<double>* dE,
    vector<vector<double>> *energydepX, vector<vector<double>> *energydepY, vector<vector<double>> *energydepZ,
    vector<vector<double>> *stepX, vector<vector<double>> *stepY, vector<vector<double>> *stepZ, vector<vector<double>> *stepT,
		vector<vector<double>>* Strue, vector<vector<double>>* Sdig, vector<double>* stampTime,
    vector<int> * trackID, vector<int> * motherID, vector<int> * PDGcode, vector<string> * process);

		void make_display(bool draw_energydepositions, bool save, double dEcut, double t1, double t2);

		vector<TH1F> fill_digitizedsignals(vector<vector<double>> *fSdig, vector<double>* fstampTime, int numbins, double min, double max);
		vector<TH1F> fill_discretizetruesignals(vector<vector<double>> *fS, int option);
		TH1F fill_integratedsignals(vector<TH1F> hsig,int numbins, double min, double max, int option);


	private:
    int frun, fsubrun, fevent;
    vector<vector<double>> *fstepX;
    vector<vector<double>> *fstepY;
    vector<vector<double>> *fstepZ;
    vector<vector<double>> *fstepT;
    vector<vector<double>> *fstepPositions= new vector<vector<double>>;
    vector<double> *fstepTimes = new vector<double>;
    vector<vector<double>> *fStrue;
    vector<vector<double>> *fSdig;
    vector<double> *fdE;
    vector<vector<double>> *fenergydepX;
    vector<vector<double>> *fenergydepY;
    vector<vector<double>> *fenergydepZ;
    vector<vector<double>> *fenergydepPositions= new vector<vector<double>>;
    vector<double> *fE;
    vector<double> *fstampTime;
    vector<int>  *ftrackID;
    vector<int>  *fmotherID;
    vector<int>  *fPDGcode;
    vector<string>  *fprocess;

    void stepXtostepPositions(double dEcut, double tmin, double tmax);
    void energydepotoenergydepPositions(double tmin, double tmax);
		vector<TMarker> fill_markersPMTs (vector< vector <double> > * Signals);
		double GetColor(double x,double time_int, double minT);
		double timeinterval(vector<double> *fstepTimes, vector<vector<double>> *fstepPositions);
		double meantrack(vector<vector<double>> *fstepPositions, int l);
		double minimumtime(vector<double> *fstepTimes, vector<vector<double>> *fstepPositions);
		vector<double> GetSignalAfterPMTQE(vector<double> v, const double PMT_QE);


		const std::vector<double> *fPMTX;
		const std::vector<double> *fPMTY;
		const std::vector<double> *fPMTZ;
		const std::vector<std::string> *fPDtype;
		const std::vector<unsigned int> *fPDSmapTPC;
		const TParameter<double> *fXCryo;
		const TParameter<double> *fYCryo;
		const TParameter<double> *fZ1Cryo;
		const TParameter<double> *fZ2Cryo;
		const TParameter<double> *fXActive;
		const TParameter<double> *fYActive;
		const TParameter<double> *fZ1Active;
		const TParameter<double> *fZ2Active;
		const TParameter<int> *fNOpCh;
    const double fExtralength_Steps=50;

    const double TimeWindowMin=-10000;
    const double TimeWindowMax=20000;
    const double dEcut_track=0.5;

		//Path Particles Colors
		int minCol=50, maxCol=100;
    const double fBaselineArapuca=-1500;
    const double fDAQSamplingTime=5.;

};



PDSEventHandle::PDSEventHandle(std::string geofile, int run, int subrun, int event,vector<double>* E, vector<double>* dE,
vector<vector<double>> *energydepX, vector<vector<double>> *energydepY, vector<vector<double>> *energydepZ,
vector<vector<double>> *stepX, vector<vector<double>> *stepY, vector<vector<double>> *stepZ, vector<vector<double>> *stepT,
vector<vector<double>>* Strue, vector<vector<double>>* Sdig, vector<double>* stampTime,
vector<int> * trackID, vector<int> * motherID, vector<int> * PDGcode, vector<string> * process){
	TFile* filegeo = TFile::Open(geofile.c_str(), "READ");
  filegeo->GetObject("PMTX", fPMTX);
	filegeo->GetObject("PMTY", fPMTY);
	filegeo->GetObject("PMTZ", fPMTZ);
	filegeo->GetObject("PDtype", fPDtype);
	filegeo->GetObject("PDSmapTPC", fPDSmapTPC);
	filegeo->GetObject("XActive", fXActive);
	filegeo->GetObject("YActive", fYActive);
	filegeo->GetObject("Z1Active", fZ1Active);
	filegeo->GetObject("Z2Active", fZ2Active);
	filegeo->GetObject("nOpCh", fNOpCh);

  /*double promz=0, promy=0;
  int nsteps;
  for(int j=0; j<fPMTZ->size(); j++){
    if(strcmp(PDSmap::PDSopch_type[j], "pmt_coated") == 0 || strcmp(PDSmap::PDSopch_type[j], "pmt_uncoated") == 0 ){
      promz+=fPMTZ->at(j);
      promy+=fPMTY->at(j);
      nsteps++;
    }
  }
  cout<<nsteps<<"  Z="<<promz/nsteps<<"   Y="<<promy/nsteps<<endl;*/

  frun=run;
  fsubrun=subrun;
  fevent=event;
  fstepX=stepX;
  fstepY=stepY;
  fstepZ=stepZ;
  fstepT=stepT;
  fE=E;
  fdE=dE;
  fenergydepX=energydepX;
  fenergydepY=energydepY;
  fenergydepZ=energydepZ;
  fmotherID=motherID;
  ftrackID=trackID;
  fPDGcode=PDGcode;
  fprocess=process;
  fStrue=Strue;
  fSdig=Sdig;
  fstampTime=stampTime;

  filegeo->Close();
}

void PDSEventHandle::stepXtostepPositions(double dEcut, double tmin, double tmax)
{
  cout<<"NParticles   StepX size: "<<" "<<fstepX->size()<<" PDGsize: "<<fPDGcode->size()<<endl;
	for(int k=0; k<fstepX->size(); k++){
    if( (fdE->at(k)>dEcut || fmotherID->at(k)==0) && (fstepT->at(k).at(0)>tmin && fstepT->at(k).at(0)<tmax) ) {
      cout<<k<<" Mother="<<fmotherID->at(k)<<" PDG="<<fPDGcode->at(k)<<" ID="<<ftrackID->at(k)<<" E="<<fE->at(k)<<" dE="<<fdE->at(k)<<" Proc:"<<fprocess->at(k);
      cout<<" X="<<fstepX->at(k).at(0)<<" Y="<<fstepY->at(k).at(0)<<" Z="<<fstepZ->at(k).at(0)<<" T="<<fstepT->at(k).at(0)<<endl;
    }
		if(fdE->at(k)<dEcut) continue;
    for(int j=0; j<fstepX->at(k).size(); j++){
      if(fstepT->at(k).at(j)<tmin || fstepT->at(k).at(j)>tmax) continue;


      fstepPositions->push_back({});
      fstepPositions->at(fstepPositions->size()-1).push_back(fstepX->at(k).at(j));
      fstepPositions->at(fstepPositions->size()-1).push_back(fstepY->at(k).at(j));
      fstepPositions->at(fstepPositions->size()-1).push_back(fstepZ->at(k).at(j));
      fstepTimes->push_back(fstepT->at(k).at(j));
    }
  }
}

void PDSEventHandle::energydepotoenergydepPositions(double tmin, double tmax)
{
  cout<<"NEnergyDepo   StepX size: "<<" "<<fenergydepX->size()<<" PDGsize: "<<fPDGcode->size()<<endl;
	for(int k=0; k<fenergydepX->size(); k++){
    if(fstepT->at(k).at(0)<tmin || fstepT->at(k).at(0)>tmax) continue;
      for(int j=0; j<fenergydepX->at(k).size(); j++){
      //if(j==0) cout<<k<<" X="<<fstepX->at(k).at(j)<<" Y="<<fstepY->at(k).at(j)<<" Z="<<fstepZ->at(k).at(j)<<" T="<<fstepT->at(k).at(j)<<" PDG="<<fPDGcode->at(k)<<" dE="<<fdE->at(k)<<endl;
      fenergydepPositions->push_back({});
      fenergydepPositions->at(fenergydepPositions->size()-1).push_back(fenergydepX->at(k).at(j));
      fenergydepPositions->at(fenergydepPositions->size()-1).push_back(fenergydepY->at(k).at(j));
      fenergydepPositions->at(fenergydepPositions->size()-1).push_back(fenergydepZ->at(k).at(j));
    }
  }
}

//_____<FILL SIGNALS FUNCTIONS
//Function to get Light Signals after the Quantum Efficiency
vector<double> PDSEventHandle::GetSignalAfterPMTQE(vector<double> v, const double PMT_QE){
	vector<double> signal_QE;
	signal_QE.clear();
	for(size_t i=0; i<v.size(); i++){
		if(gRandom->Uniform(0.,1.) <= PMT_QE)
			signal_QE.push_back(v.at(i));
	}
	return signal_QE;
}

/*/Returns a vector of TH1 with digitized signal in each PMT
vector<TH1F> fill_digitizedsignals(vector<vector<double>> *fSdig, vector<double>* fstampTime, vector<double> fPMTX, int numbins, double min, double max){
	vector<TH1F> hsig; hsig.clear();
	//Vector initialization
	for(int OC=0; OC<nOpChannels; OC++){
		TH1F h("","",numbins,min,max);
		h.GetXaxis()->SetTitle("Bins");
		h.GetYaxis()->SetTitle("ADC");
		hsig.push_back(h);
	}
	//Filling
	for(int OC=0; OC<nOpChannels; OC++){
		double stamptime=fstampTime->at(OC)/2.;
		for(int j=0; j<fSdig->at(OC).size(); j++){
			hsig.at(OC).Fill(j+stamptime, -fSdig->at(OC).at(j));
			hsig.at(OC).Fill(j+stamptime, 8000);
		}
	}

	return hsig;
}*/
//Returns a vector of TH1 with digitized signal in each PMT
vector<TH1F> PDSEventHandle::fill_digitizedsignals(vector<vector<double>> *fSdig, vector<double>* fstampTime, int numbins, double min, double max){
  vector<TH1F> hsig; hsig.clear();
  //Vector initialization
  for(int OC=0; OC<fSdig->size(); OC++){
    TH1F h("","",numbins,min,max);
    h.GetXaxis()->SetTitle("Bins");
    h.GetYaxis()->SetTitle("ADC");
    hsig.push_back(h);
  }
  //Fillingcout<<"Size "<<fSdig->size()<<endl;
  for(int OC=0; OC<fSdig->size(); OC++){
    double stamptime=fstampTime->at(OC)/fDAQSamplingTime;
    cout<<OC<<" "<<fSdig->at(OC).size()<<endl;
    for(int j=0; j<fSdig->at(OC).size(); j++){
      hsig.at(OC).Fill(stamptime+j, -fSdig->at(OC).at(j)+fBaselineArapuca);
      //cout<<j<<"  "<<-fSdig->at(OC).at(j)<<" ";
      //hsig.at(OC).Fill(j+stamptime, 8000);
    }
  }
  return hsig;
}

//Returns a vector of TH1 with true signal in each PMT
vector<TH1F> PDSEventHandle::fill_discretizetruesignals(vector<vector<double>> *fS, int option){
	vector<TH1F> hsig; hsig.clear();
  double min=minTime2(fS);
  double max=maxTime2(fS);
  int numbins=(int)(max-min);
  cout<<min<<" "<<max<<endl;
	//Vector initialization
	for(int OC=0; OC<fS->size(); OC++){
		TH1F h("",";time_slot;# photons",numbins,min,max);
    h.SetTitle(("OpCh "+to_string(OC)).c_str());
		hsig.push_back(h);
	}
  for(int OC=0; OC<fS->size(); OC++){
    for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
  }
	/*/Filling
	if(option==0){
		for(int OC=0; OC<nOpChannels; OC++){
			if(abs(fPMTX.at(OC))!=xpos_PMT) continue;
			for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
		}
	}
	else if(option==1){
		for(int OC=0; OC<nOpChannels; OC++){
			if(abs(fPMTX.at(OC))!=xpos_PMT) continue;
			if(fPMTY.at(OC)==130 || fPMTY.at(OC)==0 || fPMTY.at(OC)==-130){
				for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
			}
		}
	}
	else if(option==2){
		for(int OC=0; OC<nOpChannels; OC++){
			if(abs(fPMTX.at(OC))!=xpos_PMT) continue;
			if(fPMTY.at(OC)!=130 && fPMTY.at(OC)!=0 && fPMTY.at(OC)!=-130){
				for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
			}
		}
	}
	return hsig;
}

//Returns a vector of TH1 with digitized signal in each PMT
vector<TH1F> fill_digitizedsignals(vector<vector<double>> *fSdig, vector<double>* fstampTime, int numbins, double min, double max){
	vector<TH1F> hsig; hsig.clear();
	//Vector initialization
	for(int OC=0; OC<nOpChannels; OC++){
		TH1F h("","",numbins,min,max);
		h.GetXaxis()->SetTitle("Bins");
		h.GetYaxis()->SetTitle("ADC");
		hsig.push_back(h);
	}
	//Filling
	for(int OC=0; OC<nOpChannels; OC++){
		double stamptime=fstampTime->at(OC)/2.;
		for(int j=0; j<fSdig->at(OC).size(); j++){
			hsig.at(OC).Fill(j+stamptime, -fSdig->at(OC).at(j));
      //cout<<j<<"  "<<-fSdig->at(OC).at(j)<<" ";
			hsig.at(OC).Fill(j+stamptime, 8000);
		}
	}*/

	return hsig;
}
/*vector<TH1F> fill_discretizetruesignals(vector<vector<double>> *fS, vector<double> fPMTX, vector<double> fPMTY, int numbins, double min, double max, int option){
	vector<TH1F> hsig; hsig.clear();
	//Vector initialization
	for(int OC=0; OC<nOpChannels; OC++){
		TH1F h("","",numbins,min,max);
		h.GetXaxis()->SetTitle("t [ns]");
		h.GetYaxis()->SetTitle("PE");
		hsig.push_back(h);
	}
	//Filling
	if(option==0){
		for(int OC=0; OC<nOpChannels; OC++){
			for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
		}
	}
	else if(option==1){
		for(int OC=0; OC<nOpChannels; OC++){
			if(fPMTY.at(OC)==130 || fPMTY.at(OC)==0 || fPMTY.at(OC)==-130){
				for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
			}
		}
	}
	else if(option==2){
		for(int OC=0; OC<nOpChannels; OC++){
			if(fPMTY.at(OC)!=130 && fPMTY.at(OC)!=0 && fPMTY.at(OC)!=-130){
				for(int j=0; j<fS->at(OC).size();j++) hsig.at(OC).Fill(fS->at(OC).at(j));
			}
		}
	}
	return hsig;
}
*/
//Returns a TH1 with integrated signal to all PMTs
//Option 0: Consider both Corner and Central PMTs
//Option 1: Consider Central PMTs
//Option 2: Consider Corner PMTs
TH1F PDSEventHandle::fill_integratedsignals(vector<TH1F> hsig,int numbins, double min, double max, int option){
	TH1F h("","",numbins,min,max);
	//Filling
	if(option==0){
		for(int OC=0; OC<hsig.size(); OC++){
			h.Add(&hsig.at(OC), 1);
		}
	}
	if(option==1){
		for(int OC=0; OC<hsig.size(); OC++){
			if(fPDtype->at(OC)!="pmt_coated" || fPDtype->at(OC)!="pmt_uncoated"){
				h.Add(&hsig.at(OC), 1);
			}
		}
	}
	else if(option==2){
		for(int OC=0; OC<hsig.size(); OC++){
			if(fPDtype->at(OC)!="pmt_coated"){
				h.Add(&hsig.at(OC), 1);
			}
		}
	}
	else if(option==3){
		for(int OC=0; OC<hsig.size(); OC++){
			if(fPDtype->at(OC)!="pmt_uncoated"){
				h.Add(&hsig.at(OC), 1);
			}
		}
	}

	else{cout << "Invalid option in function fill_Intdigitizedsignals" << endl;}

	return h;
}
//_____>

//---<FUNCTION meantrack(): computes mean value of fstepPositions (l=0 mean x, l=1 mean y, l=2 mean z)
double PDSEventHandle::meantrack(vector<vector<double>> *fstepPositions, int l){
	int nSteps=fstepPositions->size();
	double sum=0;
	int num_steps=0;
	if(l==0 || l==1 || l==2){
		for(int j=0; j<nSteps; j++) {
			//Bypassing tracks outside the detector
		  	if(fstepPositions->at(j).at(0) < -fXActive->GetVal() || fstepPositions->at(j).at(0) > +fXActive->GetVal()
				|| fstepPositions->at(j).at(1) < -fYActive->GetVal() || fstepPositions->at(j).at(1) > +fYActive->GetVal()
					|| fstepPositions->at(j).at(2) < fZ1Active->GetVal() || fstepPositions->at(j).at(2) > fZ2Active->GetVal()) continue;

			sum+=fstepPositions->at(j).at(l);
			num_steps++;
		}
		return sum/num_steps;
	}
	else {
		cout << "Invalid option in function meantrack" << endl;
		return 0;
	}
}
//--->

//---<FUNCTION timeinterval(): computes time spent by the particle inside the detector
double PDSEventHandle::timeinterval(vector<double> *fstepTimes, vector<vector<double>> *fstepPositions){
	int nSteps=fstepTimes->size();
	double max=-1e30, min=1e30;

	for(int j=0; j<nSteps; j++) {
    //Bypassing tracks outside the detector
    if(fstepPositions->at(j).at(0) < -fXActive->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(0) > +fXActive->GetVal()+fExtralength_Steps
    || fstepPositions->at(j).at(1) < -fYActive->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(1) > +fYActive->GetVal()+fExtralength_Steps
      || fstepPositions->at(j).at(2) < fZ1Active->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(2) > fZ2Active->GetVal()+fExtralength_Steps) continue;

		if(fstepTimes->at(j)>max) max=fstepTimes->at(j);
		if(fstepTimes->at(j)<min) min=fstepTimes->at(j);
	}
	return max-min;
}
//--->

//---<FUNCTION minimumtime(): computes minimum step time of the particle inside the detector
double PDSEventHandle::minimumtime(vector<double> *fstepTimes, vector<vector<double>> *fstepPositions){
	int nSteps=fstepTimes->size();
	double min=1e30;

	for(int j=0; j<nSteps; j++) {
	  //Bypassing tracks outside the detector
	  	if(fstepPositions->at(j).at(0) < -fXActive->GetVal() || fstepPositions->at(j).at(0) > +fXActive->GetVal()
			|| fstepPositions->at(j).at(1) < -fYActive->GetVal() || fstepPositions->at(j).at(1) > +fYActive->GetVal()
				|| fstepPositions->at(j).at(2) < fZ1Active->GetVal() || fstepPositions->at(j).at(2) > fZ2Active->GetVal()) continue;


		if(fstepTimes->at(j)<min) min=fstepTimes->at(j);
	}
	return min;
}
//--->

double PDSEventHandle::GetColor(double x,double time_int, double minT){
  double color = minCol + (maxCol- minCol)*(x-minT)/time_int;
  return color;
}

//Rellena los TMarkers de los PMT's con tamaño según la cantidad de señal recibida
vector<TMarker> PDSEventHandle::fill_markersPMTs (vector< vector <double> > * Signals){
	vector<TMarker> Markers;


  double maxsizePMT=0, maxsizeAR=0;
  /*for(int OC=0; OC<fNOpCh->GetVal(); OC++) {
    if(fPDtype->at(OC)=="pmt_uncoated" || fPDtype->at(OC)=="xarapuca_vis"){
      //if( Signals->at(OC).size() > maxsizePMT) maxsizePMT=(int)Signals->at(OC).size();
			maxsizePMT+= Signals->at(OC).size();
    }
    else{
      //if( Signals->at(OC).size() > maxsizeAR) maxsizeAR=(int)Signals->at(OC).size();
			maxsizeAR+= Signals->at(OC).size();
    }
  }*/
	double fPMTInDAQSize, factor;

	//Loop over the Optical Channels
	for(int OC=0; OC<fNOpCh->GetVal(); OC++) {
    //PMT Marker
    TMarker marker(fPMTZ->at(OC),fPMTY->at(OC),24);
    fPMTInDAQSize = 1;
    factor=maxsizePMT;
    //ARAPUCA Marker
    if(fPDtype->at(OC)!="pmt_uncoated" && fPDtype->at(OC)!="pmt_coated"){
    	marker.SetMarkerStyle(21);
      fPMTInDAQSize=0.7;
      factor=maxsizeAR;
    }
    //Selecting Color
    int color;
    if(fPDtype->at(OC)=="pmt_uncoated" || fPDtype->at(OC)=="xarapuca_vis") color=kRed;
		else if( fPDtype->at(OC)=="xarapuca_vuv" ) color=kBlue;
    else color=kBlack;
    marker.SetMarkerColorAlpha(color, 0.65);

    //Selecting Size
		double size=0.05;// = fPMTInDAQSize;
	  if(Signals->size()==fNOpCh->GetVal()){
      if(Signals->at(OC).size() > 0 ) size = 0.05+fPMTInDAQSize*fmax(0.,log10( Signals->at(OC).size() ) );
    }
	  marker.SetMarkerSize(size);
		//cout<<OC<<" "<<Signals->at(OC).size()<<" "<<size<<endl;

	  Markers.push_back(marker);
	}
	return Markers;
}


//ARGUMENTS:
//Event info: structInfoTrack ITr
//Track Info: StepTimes (t0, t1, t2...) and StepPositions [(x(t0), y(t0), z(t0)), (x(t1), y(t1), z(t1))...]
//True light signal: Strue (OpChan1, OpChan2, OpChan3...)
//Digitized light signal: Sdig (OpChan1, OpChan2, OpChan3...) and stamp time
/*void PDSEventHandle::make_display(struct structInfoTrack ITr, vector<double>* fstepTimes, vector<vector<double>>* fstepPositions,
vector<vector<double>>* fStrue, vector<vector<double>>* fSdig, vector<double>* fstampTime,
vector<int> * ftrackID, vector<int> * fmotherID, vector<int> * fPDGcode, bool save){*/

void PDSEventHandle::make_display(bool draw_energydepositions, bool save, double dEcut, double t1, double t2){
	cout << "-*-*-*PDS DISPLAY-*-*-*"<<endl;
  cout<<fstepX->size()<<" "<<fstepY->size()<<" "<<fstepZ->size()<<" "<<fstepT->size()<<" "<<fPDGcode->size()<<" "<<fdE->size()<<" "<<fmotherID->size()<<endl;

  stepXtostepPositions(dEcut_track, TimeWindowMin, TimeWindowMax);
  energydepotoenergydepPositions(TimeWindowMin, TimeWindowMax);

  delete gROOT->FindObject("cdisplay");
	TCanvas * cdisplay = new TCanvas("cdisplay","Event Display",0,0,canvassizex,canvassizey);
	TPad *pad0 = new TPad("pad0","This is pad1",0,0,1,1);
	TPad *pad1 = new TPad("pad1","This is pad1",0.005,0.5,0.51,0.99);
  TPad *pad2 = new TPad("pad2","This is pad2",0.005,0.01,0.51,0.5);
	TPad *pad3 = new TPad("pad3","This is pad3",0.5,0.5,0.995,0.99);
  TPad *pad4 = new TPad("pad4","This is pad4",0.5,0.01,0.995,0.5);
	pad0->SetFillColor(12);pad0->Draw();pad0->Draw();pad1->Draw();pad2->Draw();pad3->Draw();pad4->Draw();
	pad1->SetFillColor(19);pad2->SetFillColor(19);
	pad3->SetFillColor(19);pad4->SetFillColor(19);
	//COMPUTING MEAN X Y AND Z
	double promx=meantrack(fstepPositions, 0);
	double promy=meantrack(fstepPositions, 1);
	double promz=meantrack(fstepPositions, 2);
	int nSteps = fstepPositions->size();
	TLatex latex;

	//---<
  //Vector with TMarker objects portraying the OpCh
	vector<TMarker> Markers; Markers.clear();
	//Filling Markers
	Markers=fill_markersPMTs(fStrue);
	//Filling the path of the particles
	double time_int = timeinterval(fstepTimes, fstepPositions);
	double minT=minimumtime(fstepTimes, fstepPositions);

	//Vector with TMarker objects portraying particle's tracks
	vector<TMarker> pathZY0, pathZY1, pathXY; pathZY0.clear();pathZY1.clear(); pathXY.clear();
	for(int j=0; j<nSteps; j++) {
		//Bypassing tracks outside the detector
	  	if(fstepPositions->at(j).at(0) < -fXActive->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(0) > +fXActive->GetVal()+fExtralength_Steps
			|| fstepPositions->at(j).at(1) < -fYActive->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(1) > +fYActive->GetVal()+fExtralength_Steps
				|| fstepPositions->at(j).at(2) < fZ1Active->GetVal()-fExtralength_Steps || fstepPositions->at(j).at(2) > fZ2Active->GetVal()+fExtralength_Steps) continue;

      if(fstepTimes->at(j)<t1 || fstepTimes->at(j)>t2) continue;


		double color = GetColor(fstepTimes->at(j), time_int, minT);

		TMarker* trackStepsZY = new TMarker(fstepPositions->at(j).at(2),fstepPositions->at(j).at(1),20);
		trackStepsZY->SetMarkerSize(.35);
		trackStepsZY->SetMarkerColor(color);
		if(fstepPositions->at(j).at(0)<0){
			pathZY0.push_back(*trackStepsZY);
		}
		else{
			pathZY1.push_back(*trackStepsZY);
		}

		TMarker* trackStepsXY = new TMarker(fstepPositions->at(j).at(0),fstepPositions->at(j).at(1),20);
		trackStepsXY->SetMarkerSize(.35);
		trackStepsXY->SetMarkerColor(color);
		pathXY.push_back(*trackStepsXY);
	}

  vector<TMarker> pathZY0_E, pathZY1_E, pathXY_E; pathZY0_E.clear();pathZY1_E.clear(); pathXY_E.clear();
  if(draw_energydepositions){
    for(int j=0; j< fenergydepPositions->size(); j++) {
      TMarker* trackStepsZY = new TMarker(fenergydepPositions->at(j).at(2),fenergydepPositions->at(j).at(1),20);
      trackStepsZY->SetMarkerSize(.05);
      trackStepsZY->SetMarkerColor(kBlack);
      if(fenergydepPositions->at(j).at(0)<0){
        pathZY0_E.push_back(*trackStepsZY);
      }
      else{
        pathZY1_E.push_back(*trackStepsZY);
      }

      TMarker* trackStepsXY = new TMarker(fenergydepPositions->at(j).at(0),fenergydepPositions->at(j).at(1),20);
      trackStepsXY->SetMarkerColor(kBlack);
      trackStepsXY->SetMarkerSize(.05);
      pathXY_E.push_back(*trackStepsXY);
    }
  }
  //--->

	//---<Subcanvas 1-> ZY plane (TPC0)
	pad1->cd();
  const double Wd=25.; //extra size to H2 histograms
  TH2D *H2ZY = new TH2D("","",500,fZ1Active->GetVal()-Wd,fZ2Active->GetVal()+Wd,400,-fYActive->GetVal()-Wd,fYActive->GetVal()+Wd);
  H2ZY->SetStats(0); //H2ZY->SetTitle("ZY plane");
  H2ZY->GetYaxis()->SetTitle("y [cm]"); H2ZY->GetXaxis()->SetTitle("z [cm]");
	//Drawing
  H2ZY->Draw("col");
  TBox activevolumeZY;
  activevolumeZY.SetFillColorAlpha(kBlue-10, 1);
  activevolumeZY.DrawBox(fZ1Active->GetVal(),-fYActive->GetVal(),fZ2Active->GetVal(),fYActive->GetVal());
	latex.SetTextSize(0.05);
	latex.DrawLatexNDC(.45,.92,"(x<0)");
	for(int jj=0; jj<Markers.size(); jj++) {
		if(fPMTX->at(jj)>0) continue;
		Markers.at(jj).Draw("same");
	}
	for(size_t ii = 0; ii < pathZY0.size(); ++ii) {
		pathZY0[ii].Draw("same");
	}
  if(draw_energydepositions){
    for(size_t ii = 0; ii < pathZY0_E.size(); ++ii) {
  		pathZY0_E[ii].Draw("same");
  	}
  }
  //--->

	//---<Subcanvas 3-> ZY plane (TPC1)
	pad3->cd();
	//Drawing
	H2ZY->Draw("col");
	activevolumeZY.SetFillColorAlpha(kBlue-10, 1);
	activevolumeZY.DrawBox(fZ1Active->GetVal(),-fYActive->GetVal(),fZ2Active->GetVal(),fYActive->GetVal());
	latex.DrawLatexNDC(.45,.92,"(x>0)");
	for(int jj=0; jj<Markers.size(); jj++) {
		if(fPMTX->at(jj)<0) continue;
		Markers.at(jj).Draw("same");
	}
	for(size_t ii = 0; ii < pathZY1.size(); ++ii) {
		pathZY1[ii].Draw("same");
	}
  if(draw_energydepositions){
    for(size_t ii = 0; ii < pathZY1_E.size(); ++ii) {
      pathZY1_E[ii].Draw("same");
    }
  }
	//--->

	//---<Subcanvas 2-> XY plane
	pad2->cd();
  TH2D *H2XY = new TH2D("","",400,-fXActive->GetVal()-Wd-50,fXActive->GetVal()+Wd+50,400,-fYActive->GetVal()-Wd,fYActive->GetVal()+Wd);
  H2XY->SetStats(0); //H2XY->SetTitle("XY plane");
  H2XY->GetYaxis()->SetTitle("y [cm]");
  H2XY->GetXaxis()->SetTitle("x [cm]");
  H2XY->Draw("colz");
  TBox activevolumeXY;
	activevolumeXY.SetFillColorAlpha(kBlue-10, 1);
	activevolumeXY.DrawBox(-fXActive->GetVal(),-fYActive->GetVal(),fXActive->GetVal(),fYActive->GetVal());
	for(size_t ii = 0; ii < pathXY.size(); ++ii) {
		pathXY[ii].Draw("same");
	}
  if(draw_energydepositions){
    for(size_t ii = 0; ii < pathXY_E.size(); ++ii) {
      pathXY_E[ii].Draw("same");
    }
  }
	//Drawing detector's cathode and anhode
	TBox cathode;
	cathode.SetFillColor(kRed);
	cathode.DrawBox(-2.5,-fYActive->GetVal(),2.5,fYActive->GetVal());
	/*TBox photocathode;
	photocathode.SetFillColor(kOrange);
	photocathode.DrawBox(212,-200,214,200);
	TBox photocathode2;
	photocathode2.SetFillColor(kOrange);
	photocathode2.DrawBox(-212,-200,-214,200);*/
  //--->


  //_____<SUBCANVAS 4: PRACTICAL INFO
	pad4->cd();
	//Drawing time interval scale
	int Ncolors=50;
		TBox ColorCode;
		double DY=0.3/(Ncolors+1);
		for(int col=0; col<=Ncolors; col++){
			ColorCode.SetFillColor(minCol+col*(maxCol-minCol)/(Ncolors));
			ColorCode.DrawBox(0.05+(col-0.5)*DY, 0.7, 0.05+(col+0.5)*DY ,0.73);
	}
	//Drawing PMTs legend
	latex.SetTextSize(0.05);
	latex.DrawLatex(.06,.64,( "#leftarrow  #Delta t="+to_string((int)time_int)+" [ns]  #rightarrow" ).c_str());
	latex.SetTextSize(0.04);
	latex.DrawLatex(.2,.48,"PMT_{coated}");
	TMarker PMTco(0.15,0.5,24); PMTco.Draw("same");
	PMTco.SetMarkerColor(kRed);PMTco.SetMarkerSize(3);
	latex.DrawLatex(.2,.38,"PMT");
	TMarker PMT(0.15,0.4,24); PMT.Draw("same");
	PMT.SetMarkerColor(kBlack);PMT.SetMarkerSize(3);
	latex.DrawLatex(.2,.23,"XARAPUCA_VUV");
	TMarker ARVUV(0.15,0.25,21); ARVUV.Draw("same");
	ARVUV.SetMarkerColor(kBlue);ARVUV.SetMarkerSize(2.5);
	latex.DrawLatex(.2,.13,"XARAPUCA_VIS");
	TMarker ARVIS(0.15,0.15,21); ARVIS.Draw("same");
	ARVIS.SetMarkerColor(kRed);ARVIS.SetMarkerSize(2.5);
	//Printing event info
  TPaveText *t=new TPaveText(.4,.1,.95,.8, "TL");
  t->AddText( ( "Run: "+to_string(frun)+", Subrun: "+to_string(fsubrun)+", Event: "+to_string(fevent) ).c_str());
  //+", E: "+to_string(ITr.E)+" [GeV]"
  t->AddText("Primary Particles: ");
  for(int k=0; k<fmotherID->size(); k++){
    //if(fmotherID->at(k)==0)
      //cout<<fPDGcode->at(k)<<endl;
        //t->AddText( ( "PDG: "+to_string(fPDGcode->at(k))).c_str());
  }
  //t->AddText( ( "Deposited Energy: "+to_string(ITr.dE)+" [MeV]" ).c_str());
  t->AddText( ( "<x>="+to_string(promx)+" [cm], <y>="+to_string(promy)+" [cm], <z>="+to_string(promz) ).c_str());
  t->Draw("same");
   latex.SetTextSize(0.1);
	 latex.SetTextColor(13);
   //latex.SetTextAlign(13);  //align at top
   latex.DrawLatex(.3,.85,"PDS Event Display");
  //_____>

	cdisplay->GetCanvas()->SetFillColor(18);
  cdisplay->Update(); cdisplay->Modified(); cdisplay->WaitPrimitive();
	if(save==true){
		pad1->SaveAs("./displayoutput/YZ.pdf", "pdf");
		pad2->SaveAs("./displayoutput/XY.pdf", "pdf");
		pad3->SaveAs("./displayoutput/signal.pdf", "pdf");
		pad4->SaveAs("./displayoutput/info.pdf", "pdf");
		cdisplay->SaveAs("./displayoutput/display.pdf", "pdf");
	}
  delete gROOT->FindObject("cdisplay");
  delete H2ZY, H2XY;pathZY0.clear();pathZY1.clear(); pathXY.clear();Markers.clear();

	return;
}


//_____<SUBCANVAS 3: Digitized signals (only PMTs)
/*pad3->cd();
double tlow, tup; int numbins_signals=1100; tlow=0; tup=numbins_signals;

vector<TH1F> h_digsignal=fill_digitizedsignals(fSdig,fstampTime,numbins_signals,tlow,tup);

TH1F h_digcentral("h_digcentral","Integrated digitized signal",numbins_signals,tlow,tup);
h_digcentral=fill_integratedsignals(h_digsignal, fPMTX, fPMTY,numbins_signals,tlow,tup, 1);
h_digcentral.SetLineColor(kRed);h_digcentral.GetXaxis()->SetTitle("Bin Number"); h_digcentral.GetYaxis()->SetTitle("ADCs");
h_digcentral.SetTitle("Integrated digitized signals (PMTs)");

TH1F h_digcorner("h_digcorner","Integrated digitized signal",numbins_signals,tlow,tup);
h_digcorner=fill_integratedsignals(h_digsignal, fPMTX, fPMTY,numbins_signals,tlow,tup, 2);
h_digcorner.SetLineColor(kBlack);h_digcorner.GetXaxis()->SetTitle("Bin Number"); h_digcorner.GetYaxis()->SetTitle("ADCs");
h_digcorner.SetTitle("Integrated digitized signals (PMTs)");

if(h_digcentral.GetMaximum()>h_digcorner.GetMaximum()){
	h_digcentral.SetStats(0);h_digcentral.Draw("hist");h_digcorner.Draw("hist same");
}
else{
	h_digcorner.SetStats(0);h_digcorner.Draw("hist");h_digcentral.Draw("hist same");
}

auto leg3 = new TLegend(0.77, 0.7, 1, 1);
leg3->AddEntry(&h_digcentral,"Central Signal (x4)", "lpf");
leg3->AddEntry(&h_digcorner,"Corner Signal", "lpf");
leg3->SetHeader("Integrated Digitized Signals");
leg3->Draw("same");
*/
//_____>
