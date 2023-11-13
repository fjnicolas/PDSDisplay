#include "include/root_mygeneralfunctionsdefinition.h"
#include "include/PDSmapping.hpp"
#include "include/OpDisplay_library.hpp"


#define frun 1
#define fsubrun 1
#define fevent 27


void MacroPDSDisplay(std::string path, bool display_daqsignals=false, bool display_ophits=false, bool option_userunsubruneventlabels=false){


  //_____<READING TREE
  TFile* fh= new TFile(path.c_str());
  TTree* tree = (TTree*)fh->Get("opanatree/OpAnaTree");
  //TTree* tree = (TTree*)fh->Get("truthphotonstree/TruthPhotonsTree");
  tree->Print();

  unsigned int run; tree->SetBranchAddress("runID",&run);
  unsigned int subrun; tree->SetBranchAddress("subrunID",&subrun);
  unsigned int event; tree->SetBranchAddress("eventID",&event);
  vector<vector<double> >* fSimPhotonsLiteVIS=new vector<vector<double> >();
  vector<vector<double> >* fSimPhotonsLiteVUV=new vector<vector<double> >();
  tree->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS);
  tree->SetBranchAddress("SimPhotonsLiteVUV",&fSimPhotonsLiteVUV);
  vector<double>* fSimPhotonsperOpChVUV=new vector<double>();
  vector<double>* fSimPhotonsperOpChVIS=new vector<double>();
  tree->SetBranchAddress("SimPhotonsperOpChVUV",&fSimPhotonsperOpChVUV);
  tree->SetBranchAddress("SimPhotonsperOpChVIS",&fSimPhotonsperOpChVUV);
  vector<vector<double>>* fSignalsDigi=new vector<vector<double>>;
  tree->SetBranchAddress("SignalsDigi",&fSignalsDigi);
  vector<double>* fStampTime=new vector<double>;
  tree->SetBranchAddress("StampTime",&fStampTime);
  vector<int>* fOpChDigi=new vector<int>;
  tree->SetBranchAddress("OpChDigi",&fOpChDigi);
  vector<vector<double>> * fstepX=new vector<vector<double>>;
  vector<vector<double>> * fstepY=new vector<vector<double>>;
  vector<vector<double>> * fstepZ=new vector<vector<double>>;
  vector<vector<double>> * fstepT=new vector<vector<double>>;
  tree->SetBranchAddress("stepX",&fstepX);
  tree->SetBranchAddress("stepY",&fstepY);
  tree->SetBranchAddress("stepZ",&fstepZ);
  tree->SetBranchAddress("stepT",&fstepT);
  vector<vector<double>> * fenergydepX=new vector<vector<double>>;
  vector<vector<double>> * fenergydepY=new vector<vector<double>>;
  vector<vector<double>> * fenergydepZ=new vector<vector<double>>;
  tree->SetBranchAddress("energydepX",&fenergydepX);
  tree->SetBranchAddress("energydepY",&fenergydepY);
  tree->SetBranchAddress("energydepZ",&fenergydepZ);
  vector<double> * fE=new vector<double>; tree->SetBranchAddress("E",&fE);
  vector<double> * fdE=new vector<double>; tree->SetBranchAddress("dE",&fdE);
  vector<int> * ftrackID=new vector<int>; tree->SetBranchAddress("trackID",&ftrackID);
  vector<int> * fmotherID=new vector<int>; tree->SetBranchAddress("motherID",&fmotherID);
  vector<int> * fPDGcode=new vector<int>; tree->SetBranchAddress("PDGcode",&fPDGcode);
  vector<string> * fprocess=new vector<string>; tree->SetBranchAddress("process",&fprocess);

  vector<double> * fophit_opch=new vector<double>;
  tree->SetBranchAddress("ophit_opch",&fophit_opch);
  vector<double> * fophit_pe=new vector<double>;
  tree->SetBranchAddress("ophit_pe",&fophit_pe);
  vector<double> * fophit_peakT=new vector<double>;
  tree->SetBranchAddress("ophit_peakT",&fophit_peakT);
  vector<double> * fophit_area=new vector<double>;
  tree->SetBranchAddress("ophit_area",&fophit_area);
  vector<double> * fophit_width=new vector<double>;
  tree->SetBranchAddress("ophit_width",&fophit_width);

  int nopflash; tree->SetBranchAddress("nopflash",&nopflash);
  vector<double>* flash_time=new vector<double>; tree->SetBranchAddress("flash_time",&flash_time);
  vector<double>* flash_total_pe=new vector<double>; tree->SetBranchAddress("flash_total_pe",&flash_total_pe);
  vector<double>* flash_y=new vector<double>; tree->SetBranchAddress("flash_y",&flash_y);
  vector<double>* flash_z=new vector<double>; tree->SetBranchAddress("flash_z",&flash_z);
  vector<double>* flash_yerr=new vector<double>; tree->SetBranchAddress("flash_yerr",&flash_yerr);
  vector<double>* flash_zerr=new vector<double>; tree->SetBranchAddress("flash_zerr",&flash_zerr);
  vector<int>* flash_tpc=new vector<int>; tree->SetBranchAddress("flash_tpc",&flash_tpc);
  vector<vector<double>>* flash_pe_v = new vector<vector<double>> ; tree->SetBranchAddress("flash_pe_v",&flash_pe_v);
  vector<vector<double>> * flash_ophit_opch=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_ch",&flash_ophit_opch);
  vector<vector<double>> * flash_ophit_pe=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_pe",&flash_ophit_pe);
  vector<vector<double>> * flash_ophit_peakT=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_time",&flash_ophit_peakT);
  vector<vector<double>> * flash_ophit_area=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_area",&flash_ophit_area);

  /*int nmcflash; tree->SetBranchAddress("nmcflash",&nmcflash);
  vector<double>* mcflash_time=new vector<double>; tree->SetBranchAddress("mcflash_time",&mcflash_time);
  vector<double>* mcflash_abstime=new vector<double>; tree->SetBranchAddress("mcflash_abstime",&mcflash_abstime);
  vector<double>* mcflash_y=new vector<double>; tree->SetBranchAddress("mcflash_y",&mcflash_y);
  vector<double>* mcflash_z=new vector<double>; tree->SetBranchAddress("mcflash_z",&mcflash_z);
  vector<double>* mcflash_yerr=new vector<double>; tree->SetBranchAddress("mcflash_yerr",&mcflash_yerr);
  vector<double>* mcflash_zerr=new vector<double>; tree->SetBranchAddress("mcflash_zerr",&mcflash_zerr);
  vector<double>* mcflash_total_pe=new vector<double>; tree->SetBranchAddress("mcflash_total_pe",&mcflash_total_pe);
  vector<int>* mcflash_tpc=new vector<int>; tree->SetBranchAddress("mcflash_tpc",&mcflash_tpc);
  vector<vector<double>>* mcflash_pe_v=new vector<vector<double>>; tree->SetBranchAddress("mcflash_pe_v",&mcflash_pe_v);*/


  vector<double>* nuvX=new vector<double>; tree->SetBranchAddress("nuvX",&nuvX);
  vector<double>* nuvY=new vector<double>; tree->SetBranchAddress("nuvY",&nuvY);
  vector<double>* nuvZ=new vector<double>; tree->SetBranchAddress("nuvZ",&nuvZ);
  vector<double>* nuvT=new vector<double>; tree->SetBranchAddress("nuvT",&nuvT);
  vector<double>* nuvE=new vector<double>; tree->SetBranchAddress("nuvE",&nuvE);

  unsigned int InTimeCosmics; tree->SetBranchAddress("InTimeCosmics",&InTimeCosmics);

  string sbndgeofile;
  sbndgeofile="include/SBNDgeometry_PDSdisplay.root";
  //sbndgeofile="include/sbndPDSMap_v02_00.root";
  string dunegeofile="include/DUNEgeometry_PDSdisplay.root";
  //_____>

  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    std::cout<<"i="<<i<<std::endl;
    if(option_userunsubruneventlabels && run!=frun && fsubrun!=subrun && fevent!=event){continue;}


    PDSEventHandle pds_evt(sbndgeofile, run, subrun, event, fE, fdE, fenergydepX, fenergydepY, fenergydepZ, fstepX, fstepY, fstepZ, fstepT, fSimPhotonsLiteVUV, fSignalsDigi, fStampTime, fOpChDigi, ftrackID, fmotherID, fPDGcode, fprocess);

    cout<<"Run="<<run<<" Subrun="<<subrun<<" Event="<<event<<endl;

    pds_evt.make_display(1, 0, 0, -1e10, 1e10, false);

    if(display_ophits){
      cout<<" nFlash="<<nopflash<<endl;
      for(int l=0; l<nopflash; l++){
        //displayOpHits_dune(*fophit_peakT, *fophit_pe, *fophit_opch, flash_ophit_peakT->at(l), flash_ophit_pe->at(l), flash_ophit_opch->at(l), 0.01, 320);
      }
    }

  



  }






  return 0;
}
