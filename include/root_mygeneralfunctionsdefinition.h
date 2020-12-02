#ifndef MYGERNEALFUNCTIONS
#define MYGERNEALFUNCTIONS

namespace CustomDisplays{
  vector<TPad*> buildpadcanvas(int nx, int ny){
    vector<TPad*> Tp;
    double x=0, y=1, dx=1./nx, dy=1./ny;
    TPad *pad = new TPad("","", 0, 0, 1, 1, -1, -1, -1);
    Tp.push_back(pad);
    for(int i=1; i<=nx; i++){
      y=1;
      for(int j=1; j<=ny; j++){
        //TPad *pad = new TPad("a", "a",x, y, x+dx,y+dy);
        //TPad pad("a", "a",x, y, x+dx,y+dy,1, 1, 2);
        //cout<<x<<" "<<y<<endl;
        TPad *pad = new TPad("","", x, y-dy, x+dx, y, -1, -1, -1);
        Tp.push_back(pad);
        y-=dy;
      }
      x+=dx;
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
    }
    return Tp;
  }

  void DrawBiasStdDev_fromTProfile(TProfile h){
    TMultiGraph * mg = new TMultiGraph("","");
    vector<double> X, Y1, Y2;
    for(int k=1; k<=h.GetNbinsX(); k++){
      X.push_back(h.GetBinCenter(k));
      Y1.push_back(h.GetBinContent(k));
      Y2.push_back(h.GetBinError(k));
      //cout<<k<<" "<<h.GetBinCenter(k)<<" "<<h.GetBinContent(k)<<endl;
    }
    TGraph *g1=new TGraph(X.size(), &X.at(0), &Y1.at(0));
    g1->SetMarkerColor(kRed); g1->SetMarkerStyle(20); g1->SetMarkerSize(0.6);
    g1->SetTitle("Bias");
    TGraph *g2=new TGraph(X.size(), &X.at(0), &Y2.at(0));
    g2->SetMarkerColor(kBlue); g2->SetMarkerStyle(21); g2->SetMarkerSize(0.6);
    g2->SetTitle("StdDev");

    mg->GetXaxis()->SetTitle( h.GetXaxis()->GetTitle() );
    mg->GetYaxis()->SetTitle( h.GetYaxis()->GetTitle() );
    mg->Add(g1);mg->Add(g2);
    mg->Draw("alp");
  }
}



#endif
