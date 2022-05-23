#include "include/root_mygeneralfunctionsdefinition.h"
#include "include/PDSmapping.hpp"
#include "include/OpDisplay_library.hpp"

#define frun 11
#define fsubrun 80
#define fevent 15



void macro_OpAnalyzer(){
  bool display_daqsignals=false; display_daqsignals=true;
  bool option_userunsubruneventlabels=true;
  string path;

  path="../../FlashRecoTest/Test2/OpAna_nusample.root";
  path="OpAna_1mu1p.root";
  path="OpAna_nusample.root";
  path="Study/OpAnaTree_DebugSaturation.root";

  path="OpAna_039765.root";
  path="OpAna_041621.root";

  path="../LightAna_09_09_01/Test2/FlashTimeResStudy/TimeRes10ns/flashfinder_tree_2500.root";

  path="../LightAna_09_09_01/Test2/NueSample/OpAna_Tree_15-20k.root";
  //_____<READING TREE
  TFile* fh= new TFile(path.c_str());
  TTree* tree = (TTree*)fh->Get("opanatree/OpAnaTree");
  tree->Print();

  unsigned int run; tree->SetBranchAddress("runID",&run);
  unsigned int subrun; tree->SetBranchAddress("subrunID",&subrun);
  unsigned int event; tree->SetBranchAddress("eventID",&event);
  vector<vector<double> >* fSimPhotonsLiteVIS=new vector<vector<double> >();
  vector<vector<double> >* fSimPhotonsLiteVUV=new vector<vector<double> >();
  tree->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS);
  tree->SetBranchAddress("SimPhotonsLiteVUV",&fSimPhotonsLiteVUV);
  vector<vector<double>>* fSignalsDigi=new vector<vector<double>>;
  tree->SetBranchAddress("SignalsDigi",&fSignalsDigi);
  vector<double>* fStampTime=new vector<double>;
  tree->SetBranchAddress("StampTime",&fStampTime);
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

  vector<int> * flash_tpc=new vector<int>;
  tree->SetBranchAddress("flash_tpc",&flash_tpc);
  vector<vector<double>> * flash_ophit_opch=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_ch",&flash_ophit_opch);
  vector<vector<double>> * flash_ophit_pe=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_pe",&flash_ophit_pe);
  vector<vector<double>> * flash_ophit_peakT=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_time",&flash_ophit_peakT);
  vector<vector<double>> * flash_ophit_area=new vector<vector<double>>;
  tree->SetBranchAddress("flash_ophit_area",&flash_ophit_area);

  string sbndgeofile="include/SBNDgeometry_PDSdisplay.root";

  std::cout<<" ---- Reading...done!"<<std::endl;
  //_____>

  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);

    if(option_userunsubruneventlabels){
      if(run!=frun || fsubrun!=subrun || fevent!=event) continue;
    }


    //PDSEventHandle pds_evt(sbndgeofile, run, subrun, event, fE, fdE, fstepX, fstepY, fstepZ, fstepT, fSimPhotonsLiteVUV,
    //  fSignalsDigi, fStampTime, ftrackID, fmotherID, fPDGcode);
    PDSEventHandle pds_evt(sbndgeofile, run, subrun, event, fE, fdE, fenergydepX, fenergydepY, fenergydepZ, fstepX, fstepY, fstepZ, fstepT, fSimPhotonsLiteVUV,
      fSignalsDigi, fStampTime, ftrackID, fmotherID, fPDGcode, fprocess);
    cout<<"Run="<<run<<" Subrun="<<subrun<<" Event="<<event<<endl;


    double TimeRes=0.01;
    cout<<"--- MinTime [us]  Flash:"<<minTime2(flash_ophit_peakT)<<" Tot:"<<minTime(*fophit_peakT)<<endl;
    cout<<"--- MaxTime [us]  Flash:"<<maxTime2(flash_ophit_peakT)<<" Tot:"<<maxTime(*fophit_peakT)<<endl;
    double mint=minTime(*fophit_peakT);
    double maxt=maxTime(*fophit_peakT);
    mint-=10*TimeRes;
    size_t nbins=1000;
    maxt=mint+nbins*TimeRes;

    //maxt+=10*TimeRes;
    //size_t nbins = (size_t)((maxt - mint) / TimeRes) ;

    TH2F h_OpHits("h_OpHits_Totales", "OpHits (all);Time [#us];OpCh",nbins, mint, maxt, 320, 0, 320);
    TH1F h_OpHitsPE("h_OpHitsPE", ";Time [#us];#PE",nbins, mint, maxt);
    TH1F h_OpHitsNOpCh("h_OpHitsNOpCh", ";Time [#us];OpCh in time coincidence",nbins, mint, maxt);
    for(int l=0; l<fophit_pe->size(); l++){
      //int OC=fophit_opch->at(l);
      //if( fPMTX->at(OC)>0 && inTPC==0) continue;
      //if( fPMTX->at(OC)<0 && inTPC==1) continue;
      //cout<<l<<" OC="<<fophit_opch->at(l)<<" PE="<<fophit_pe->at(l)<<" T="<<fophit_peakT->at(l)<<" Area="<<fophit_area->at(l)<<endl;
      //size_t index=(size_t)((fophit_peakT->at(l) - mint) / TimeRes);
      //h_OpHits.Fill( index , fophit_opch->at(l), fophit_pe->at(l));
      h_OpHits.Fill(  fophit_peakT->at(l) , fophit_opch->at(l), fophit_pe->at(l));
      h_OpHitsPE.Fill( fophit_peakT->at(l) , fophit_pe->at(l));
      h_OpHitsNOpCh.Fill( fophit_peakT->at(l) );
    }
    TH1F h_OpHitsPEmin3("h_OpHitsPEmin3", ";BinNumber (10ns);#PE",nbins, mint, maxt);
    TimeRes=0.02;maxt=mint+nbins*TimeRes;
    TH1F h_OpHitsPEnewTimeRes("h_OpHitsPEnewTimeRes", ";BinNumber (10ns);#PE",nbins, mint, maxt);
    for(int j=1; j<=h_OpHitsPE.GetNbinsX(); j++){
      if(h_OpHitsNOpCh.GetBinContent(j)>=3) h_OpHitsPEmin3.SetBinContent(j, h_OpHitsPE.GetBinContent(j));
      h_OpHitsPEnewTimeRes.Fill(h_OpHitsPE.GetBinCenter(j), h_OpHitsPE.GetBinContent(j));
    }

    TH2F h_OpHits2("h_OpHits2", "OpHits (from OpFlash);Time [#us];OpCh",nbins, mint, maxt, 320, 0, 320);
    TH1F h_OpHits2PE("h_OpHits2PE", "OpHits (from OpFlash);Time [#us];#PE",nbins, mint, maxt);
    TH1F h_OpHits2NOpCh("h1_OpHits2NOpCh", "OpHits (from OpFlash);Time [#us];NOpCh in time coincidence",nbins, mint, maxt);
    for(int l=0; l<flash_ophit_pe->at(0).size(); l++){
      //cout<<l<<" OC="<<flash_ophit_opch->at(0).at(l)<<" PE="<<flash_ophit_pe->at(0).at(l)<<" T="<<flash_ophit_peakT->at(0).at(l)<<" Area="<<flash_ophit_area->at(0).at(l)<<endl;
      //size_t index=(size_t)((flash_ophit_peakT->at(0).at(l) - mint) / TimeRes);
      h_OpHits2.Fill(flash_ophit_peakT->at(0).at(l), flash_ophit_opch->at(0).at(l), flash_ophit_pe->at(0).at(l));
      h_OpHits2PE.Fill(flash_ophit_peakT->at(0).at(l),flash_ophit_pe->at(0).at(l));
      h_OpHits2NOpCh.Fill( flash_ophit_peakT->at(0).at(l) );
    }
    TH1F h_OpHits2PEmin3("h_OpHits2PEmin3", "OpHits (from OpFlash);BinNumber (10ns);#PE",nbins, mint, maxt);
    for(int j=1; j<=h_OpHitsPE.GetNbinsX(); j++){
      if(h_OpHits2NOpCh.GetBinContent(j)>=3) h_OpHits2PEmin3.SetBinContent(j, h_OpHits2PE.GetBinContent(j));
    }

    int maxbin=h_OpHitsPE.GetMaximumBin();
    cout<<"MaxBin="<<maxbin<<endl;

    TCanvas cOH("cOH","",100,100,1200,650);vector<TPad*>  TOH=CustomDisplays::buildpadcanvas(2,2);
    gStyle->SetOptStat(0);

    TOH[1]->cd(); gPad->SetGrid();
    h_OpHits.Draw("colz");
    h_OpHits.GetXaxis()->SetRangeUser(0, 4);
    TOH[3]->cd(); gPad->SetGrid();
    h_OpHitsPE.Draw("hist");h_OpHitsPEmin3.Draw("hist same");h_OpHitsPEnewTimeRes.Draw("hist same");
    h_OpHitsPEmin3.SetLineColor(kRed);h_OpHitsPEnewTimeRes.SetLineColor(kGreen);
    h_OpHitsPE.GetXaxis()->SetRangeUser(0, 4);
    auto leg= new TLegend(0.75, 0.7, 0.95, 0.97);
    leg->AddEntry(&h_OpHitsPE,"All Bins", "l");
    leg->AddEntry(&h_OpHitsPEmin3,"#splitline{>3OpCh in time}{coincidence}", "l");
    leg->AddEntry(&h_OpHitsPEnewTimeRes,"All Bins (TimeRes=20ns)", "l");
    leg->Draw("same");

    TOH[2]->cd(); gPad->SetGrid();
    h_OpHits2.Draw("colz");
    h_OpHits2.GetXaxis()->SetRangeUser(0, 4);

    TOH[4]->cd(); gPad->SetGrid();
    h_OpHits2PE.Draw("hist");h_OpHits2PEmin3.Draw("hist same");
    h_OpHits2PEmin3.SetLineColor(kRed);h_OpHits2PEmin3.SetLineStyle(kDashed);
    h_OpHits2PE.GetXaxis()->SetRangeUser(0, 4);
    leg->Draw("same");
    cOH.cd();cOH.Update(); cOH.WaitPrimitive();


    //for(int j=0; j<fSimPhotonsLiteVIS->size(); j++)
      //cout<<j<<" "<<fSimPhotonsLiteVIS->at(j).size()<<" VUV:"<<fSimPhotonsLiteVUV->at(j).size()<<endl;

    pds_evt.make_display(1,0,.5, -100, 10000);

    if(display_daqsignals){
      int min=-100, max=2000;
      cout<<"SIZE  "<<fSignalsDigi->size()<<endl;
      vector<TH1F> digisig=pds_evt.fill_digitizedsignals(fSignalsDigi, fStampTime, max-min, min,max);
      TCanvas c("c","",100,100,1200,650);vector<TPad*>  T=CustomDisplays::buildpadcanvas(1,1);
      int OC=-1;
      cout<<endl<<"Select OpCh: ";cin>>OC;
      while(OC>0){
        T[1]->cd();gPad->SetGrid();
        digisig.at(OC).Draw("hist");
        c.cd();c.Update(); c.WaitPrimitive();
        cout<<endl<<"Select OpCh: ";
        cin>>OC;
      }
      T[1]->cd();gPad->SetGrid();
      //TH1F digisig_allOC=pds_evt.fill_integratedsignals(digisig,max-min, min,max, 0);digisig_allOC.Draw("hist");
      c.cd();c.Update(); c.WaitPrimitive();
    }

  }






  return 0;
}
