#include "../include/PDSmapping.hpp"

void make_DUNEgeometry(){
  TFile* fileout = TFile::Open("DUNEgeometry_PDSdisplay.root", "RECREATE");

  //_____<OpCh POSITIONS{
    //Almacenamos las posiciones de los PMT del fichero y definimos QE
      string positions = "DUNEgeometry_PDSdisplay.txt";
      //Getting the pmt positions (y and z)
      ifstream Traks_file1(positions.c_str());
      if(!Traks_file1) cerr << "WARNING:  Failed to open file OpChPositions File"<< endl;
      Traks_file1.seekg(0);
      vector<double> fPMTX;
      vector<double> fPMTY;
      vector<double> fPMTZ;
      vector<string> fPDtype;
      vector<unsigned int> fPDSmapTPC;
      string t;
      double id, x, y, z;
      while(!(Traks_file1.eof())) {
        Traks_file1 >> id >> x >> y >> z >> t;
        cout<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<t<<endl;
        fPMTX.push_back(x);
        fPMTY.push_back(y);
        fPMTZ.push_back(z);
        fPDtype.push_back(t);
        fPDSmapTPC.push_back(x/abs(x));
      }
  //_____>
  fileout->WriteObject(&fPMTX, "PMTX");
  fileout->WriteObject(&fPMTY, "PMTY");
  fileout->WriteObject(&fPMTZ, "PMTZ");
  fileout->WriteObject(&fPDtype, "PDtype");
  fileout->WriteObject(&fPDSmapTPC, "PDSmapTPC");

  //Detector size
  TParameter<double> xdet_size("xdet_size", 400);
  TParameter<double> ydet_size("ydet_size", 600);
  TParameter<double> zmindet_size("zmindet_size", -100);
  TParameter<double> zmaxdet_size("zmaxdet_size", 1500);
  TParameter<int> nOpCh("nOpCh", 480);

  fileout->WriteObject(&xdet_size, "XActive");
  fileout->WriteObject(&ydet_size, "YActive");
  fileout->WriteObject(&zmindet_size, "Z1Active");
  fileout->WriteObject(&zmaxdet_size, "Z2Active");
  fileout->WriteObject(&nOpCh, "nOpCh");

  fileout->Close();

  TFile* filein = TFile::Open("SBNDgeometry_PDSdisplay.root", "READ");
  std::vector<double> *vin;
  filein->GetObject("PMTX", vin);
  //cout<<vin.size()<<" "<<vin.at(45);
  TParameter<double> *pp;
  filein->GetObject("xdet_size", pp);
  cout<<" cvdsfs  "<<pp->GetVal();
  //pp=( (TParameter<double>*) filein->Get("xpos") );//->GetVal();

  return 0;
}
