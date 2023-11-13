#include "/Applications/root_v6.16.00/include/root_mygeneralfunctionsdefinition.h"
#include "../include/PDSmapping.hpp"
#include "../include/OpDisplay_library.hpp"

void display_opflashshape(vector<vector<double>> v, int option, string title){
  //delete gROOT->FindObject("c");TCanvas * c = new TCanvas("c","Op Flash",100,100,1600,850);

  vector<double> x, y;
  x.clear();y.clear();
  for(size_t OC=0; OC<v.size(); OC++){
    x.push_back(OC);
    y.push_back(v.at(OC).size());
    //cout << OC << " "<<v.at(OC)<<endl;
  }
  TGraph *g = new TGraph(x.size(), &x.at(0), &y.at(0));
  g->SetMarkerStyle(20); g->SetMarkerSize(0.7);g->SetMarkerColor(kBlue);g->SetLineColor(kBlue);
  g->GetYaxis()->SetRangeUser(0,g->GetMaximum());
  g->GetYaxis()->SetTitle("# Photons");g->GetXaxis()->SetTitle("OpCh");
  g->SetTitle(title.c_str());
  g->Draw("alp");

  return;
}


void macro_check(){
  TFile* fh1 = new TFile("OpAna_Tree_lite.root");
  TTree* tree1 = (TTree*)fh1->Get("opanatree/OpAnaTree");
  TFile* fh2 = new TFile("OpAna_Tree_nolite.root");
  TTree* tree2 = (TTree*)fh2->Get("opanatree/OpAnaTree");
  tree1->Print();tree2->Print();

  vector<vector<double> >* fSimPhotonsLiteVIS_Semi=new vector<vector<double> >();
  vector<vector<double> >* fSimPhotonsLiteVUV_Semi=new vector<vector<double> >();
  vector<vector<double> >* fSimPhotonsLiteVIS_OpLib=new vector<vector<double> >();
  vector<vector<double> >* fSimPhotonsLiteVUV_OpLib=new vector<vector<double> >();
  tree1->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS_Semi);
  tree1->SetBranchAddress("SimPhotonsLiteVUV",&fSimPhotonsLiteVUV_Semi);
  tree2->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS_OpLib);
  tree2->SetBranchAddress("SimPhotonsLiteVUV",&fSimPhotonsLiteVUV_OpLib);

  tree1->GetEntry(3);
  tree2->GetEntry(3);
  cout<<fSimPhotonsLiteVUV_OpLib->size()<<endl;
  cout<<fSimPhotonsLiteVIS_OpLib->size()<<endl;
  cout<<fSimPhotonsLiteVUV_Semi->size()<<endl;
  cout<<fSimPhotonsLiteVIS_Semi->size()<<endl;


  TH1F h_VUV_OpLib("h_VUV_OpLib", " SimPhotonsLite VUV (OpLib);Time Slot;# PE", 6e6, -3e6, 3e6);
  TH1F h_VIS_OpLib("h_VIS_OpLib", " SimPhotonsLite VIS (OpLib);Time Slot;# PE", 6e6, -3e6, 3e6);
  TH1F h_VUV_Semi("h_VUV_Semi", " SimPhotonsLite VUV (Semianalytic);Time Slot;# PE", 6e6, -3e6, 3e6);
  TH1F h_VIS_Semi("h_VIS_Semi", " SimPhotonsLite VIS (Semianalytic);Time Slot;# PE", 6e6, -3e6, 3e6);
  for(int i=0; i<fSimPhotonsLiteVUV_OpLib->size(); i++){
    for(int j=0; j<fSimPhotonsLiteVUV_OpLib->at(i).size(); j++){
      h_VUV_OpLib.Fill(fSimPhotonsLiteVUV_OpLib->at(i).at(j));
    }
  }
  for(int i=0; i<fSimPhotonsLiteVIS_OpLib->size(); i++){
    for(int j=0; j<fSimPhotonsLiteVIS_OpLib->at(i).size(); j++){
      h_VIS_OpLib.Fill(fSimPhotonsLiteVIS_OpLib->at(i).at(j));
    }
  }
  for(int i=0; i<fSimPhotonsLiteVUV_Semi->size(); i++){
    for(int j=0; j<fSimPhotonsLiteVUV_Semi->at(i).size(); j++){
      h_VUV_Semi.Fill(fSimPhotonsLiteVUV_Semi->at(i).at(j));
    }
  }
  for(int i=0; i<fSimPhotonsLiteVIS_Semi->size(); i++){
    for(int j=0; j<fSimPhotonsLiteVIS_Semi->at(i).size(); j++){
      h_VIS_Semi.Fill(fSimPhotonsLiteVIS_Semi->at(i).at(j));
    }
  }



  TCanvas c("c","",100,100,1200,650);
    vector<TPad*>  T=CustomDisplays::buildpadcanvas(2,2);

    T[1]->cd();gPad->SetGrid();
    display_opflashshape(*fSimPhotonsLiteVUV_Semi, 0, "VUV Semi");

    T[2]->cd();gPad->SetGrid();
    display_opflashshape(*fSimPhotonsLiteVIS_Semi, 0, "VIS Semi");

    T[3]->cd();gPad->SetGrid();
    display_opflashshape(*fSimPhotonsLiteVUV_OpLib, 0, "VUV OpLib");

    T[4]->cd();gPad->SetGrid();
    display_opflashshape(*fSimPhotonsLiteVIS_OpLib, 0, "VIS OpLib");

    c.Update(); c.WaitPrimitive();


    TCanvas c1("c1","",100,100,650,650);
    vector<TPad*>  T1=CustomDisplays::buildpadcanvas(1,2);

    T1[1]->cd();gPad->SetGrid();
    h_VUV_OpLib.SetStats(0);h_VUV_Semi.SetStats(0);h_VUV_OpLib.Draw();h_VUV_Semi.Draw("same");
    h_VUV_Semi.SetLineColor(kRed);
    auto leg1 = new TLegend(0.8, 0.8, 1, 1);
    leg1->AddEntry(&h_VUV_OpLib, "VUV OpLib", "l");
    leg1->AddEntry(&h_VUV_Semi, "VUV Semi", "l");
      leg1->SetTextSize(0.035);leg1->Draw("same");
    h_VUV_OpLib.SetTitle("");

    T1[2]->cd();gPad->SetGrid();
    h_VIS_OpLib.SetStats(0);h_VIS_Semi.SetStats(0);h_VIS_OpLib.Draw();h_VIS_Semi.Draw("same");
    h_VIS_Semi.SetLineColor(kRed);
    auto leg2 = new TLegend(0.8, 0.8, 1, 1);
    leg2->AddEntry(&h_VIS_OpLib, "VIS OpLib", "l");
    leg2->AddEntry(&h_VIS_Semi, "VIS Semi", "l");
      leg2->SetTextSize(0.035);leg2->Draw("same");
    h_VIS_OpLib.SetTitle("");

    T1[0]->cd();c1.Update(); c1.WaitPrimitive();





  vector<TH1F> truesignals=fill_discretizetruesignals(fSimPhotonsLiteVUV_Semi, 0);
  int OC; cout<<"Select OpCh: ";cin>>OC;
  while(OC>=0 && OC<PDSmap::nOpChannels){
    TCanvas cOC("cOC","",100,100,1200,650);
    cOC.cd();
    truesignals.at(OC).Draw("hist");
    cOC.Update(); cOC.WaitPrimitive();
    cout<<"Select OpCh: ";cin>>OC;
  }



  /*T[1]->cd();gPad->SetGrid();
  h_VUV_OpLib.Draw();
  T[3]->cd();gPad->SetGrid();
  h_VIS_OpLib.Draw();
  T[2]->cd();gPad->SetGrid();
  h_VUV_Semi.Draw();
  T[4]->cd();gPad->SetGrid();
  h_VIS_Semi.Draw();
  T[0]->cd();c.Update(); c.WaitPrimitive();
  */







  return 0;
}
