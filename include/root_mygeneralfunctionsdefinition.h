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

  void DrawBiasStdDev_fromTProfile(TProfile h, int minentries_perbin, double miny, double maxy){
    TMultiGraph * mg = new TMultiGraph("","");
    vector<double> X, Y1, Y2;
    for(int k=1; k<=h.GetNbinsX(); k++){
      if(h.GetBinEntries(k)<minentries_perbin) continue;
      X.push_back(h.GetBinCenter(k));
      Y1.push_back(h.GetBinContent(k));
      Y2.push_back(h.GetBinError(k));
      //cout<<k<<" "<<h.GetBinCenter(k)<<" "<<h.GetBinContent(k)<<" "<<h.GetBinEntries(k)<<endl;
    }
    TGraph *g1=new TGraph(X.size(), &X.at(0), &Y1.at(0));
    g1->SetMarkerColor(kRed); g1->SetMarkerStyle(20); g1->SetMarkerSize(0.6);
    g1->SetTitle("Bias");
    TGraph *g2=new TGraph(X.size(), &X.at(0), &Y2.at(0));
    g2->SetMarkerColor(kBlue); g2->SetMarkerStyle(21); g2->SetMarkerSize(0.6);
    g2->SetTitle("StdDev");

    mg->SetMinimum(miny);mg->SetMaximum(maxy);
    mg->GetXaxis()->SetTitle( h.GetXaxis()->GetTitle() );
    mg->GetYaxis()->SetTitle( h.GetYaxis()->GetTitle() );
    mg->SetTitle( h.GetTitle() );
    mg->Add(g1);mg->Add(g2);
    mg->Draw("alp");

  }
}

//_______Function to interpolate two values (used in xreconstrution)
double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )
{
	int size = xData.size();
	int i = 0;                                          // find left end of interval for interpolation
	if ( x >= xData[size - 2] )                         // special case: beyond right end
	{
	i = size - 2;
	}
	else
	{
		while ( x > xData[i+1] ) i++;
	}
	double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
	if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
	{
	if ( x < xL ) yR = yL;
	if ( x > xR ) yL = yR;
	}
	double dydx = ( yR - yL ) / ( xR - xL );            // gradient
	return yL + dydx * ( x - xL );                      // linear interpolation
}



#endif
