void plot(char* input_filename){
  gStyle->SetPalette(1);
  gStyle->SetOptStat(11111);
//     TFile *file=new TFile("output.root");
//   TFile *file=new TFile("output_neutron.root");
//     TFile *file=new TFile("output_proton.root");
//     TFile *file=new TFile("output_deuteron.root");
//       TFile *file=new TFile("output_3he.root");
//       TFile *file=new TFile("output_3he_5.7GeV.root");
  
  TFile *file=new TFile(input_filename);
  
      if (file->IsZombie()) {
       cout << "Error opening file" << input_filename << endl;
       exit(-1);
    }
    else cout << "open file " << input_filename << endl;

      
TTree *T = (TTree*) file->Get("T");

TCanvas *c_rate_Q2x = new TCanvas("rate_Q2x","rate_Q2x",1200,900);
// c_rate->Divide(1,2);
// c_rate->cd(1);
gPad->SetLogz();
// T->Draw("Q2:x>>h1(100,0,1,120,0,12)","rate","colz");
// // T->Draw("Q2:x>>h1(100,0,1,120,0,12)","rate*(9./180.*3.14159<theta && theta<24./180.*3.14159)","colz");
// T->Draw("Q2:x>>h1(11,0.225,0.775,10,3,12)","rate*(9./180.*3.14159<theta && theta<24./180.*3.14159)","TEXT");
// T->Draw("Q2:x>>h1(11,0.225,0.775,10,3,12)","rate*( 9./180.*3.14159<theta && theta<24./180.*3.14159)","TEXT");
// T->Draw("Q2:x>>h1(11,0.225,0.775,10,3,12)","rate*(W>2 && Q2 >1 && 22./180.*3.14159<theta && theta<35./180.*3.14159)","TEXT");
T->Draw("Q2:x>>h1(11,0.225,0.775,10,3,13)","rate*(W>2 && Q2 >1 && 9./180.*3.14159<theta && theta<24./180.*3.14159)","TEXT");

//       dimension AXSOLMIN(11),AXSOLMAX(11),ARATESOL(5,11),AQMU2SOL(5,11) ! in Hz
//       dimension NSOL(11)
//       DATA AXSOLMIN /0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 
//      >     0.55, 0.60, 0.65, 0.70, 0.75/
//       DATA AXSOLMAX /0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 
//      >     0.60, 0.65, 0.70, 0.75, 0.80/
//       DATA NSOL   /  2,    3,    3,    4,    3,    4,
//      >      5,    4,    5,    4,    2/ !3/
//       DATA AQMU2SOL / 
//      >     3.5, 4.5, 0.0, 0.0, 0.0, ! 0.25
//      >     3.5, 4.5, 5.5, 0.0, 0.0, ! 0.30
//      >     4.5, 5.5, 6.5, 0.0, 0.0, ! 0.35
//      >     4.5, 5.5, 6.5, 7.5, 0.0, ! 0.40
//      >     5.5, 6.5, 7.5, 0.0, 0.0, ! 0.45
//      >     5.5, 6.5, 7.5, 8.5, 0.0, ! 0.50
//      >     5.5, 6.5, 7.5, 8.5, 9.5, ! 0.55
//      >     6.5, 7.5, 8.5, 9.5, 0.0, ! 0.60
//      >     6.5, 7.5, 8.5, 9.5,10.5, ! 0.65
//      >     7.5, 8.5, 9.5,10.5, 0.0, ! 0.70
//      >     9.5,10.5,11.5, 0.0, 0.0/ ! 0.75
// !          qmu2=0.0
// !          rate=0.0
// !          do iq2=1,NSOL(ix)
// !             print *,'x=',xbj,' iq2=',iq2,' rate=',ARATESOL(iq2,ix)
// !             qmu2 = qmu2 + AQMU2SOL(iq2,ix)*ARATESOL(iq2,ix)
// !             rate = rate + ARATESOL(iq2,ix)
// !          enddo
// !          qmu2 = qmu2/rate

const int xbin=11;
const int Q2bin=10; 
 double x_low=0.225,x_high=0.775;
  double Q2_low=3.,Q2_high=13.;
  double Q2width=(Q2_high-Q2_low)/Q2bin;
 double Q2[xbin],rate[xbin];
for (int i=0;i<xbin;i++){
  Q2[i]=0.;
  rate[i]=0.;
  for (int j=0;j<Q2bin;j++){
  Q2[i] = Q2[i] + (Q2_low+Q2width*(j+1./2))*h1->GetBinContent(i+1,j+1);
  rate[i] = rate[i] + h1->GetBinContent(i+1,j+1);
  }
  Q2[i] = Q2[i]/rate[i];  
}
for (int i=0;i<xbin;i++) cout << Q2[i] << ",";
cout << endl;
for (int i=0;i<xbin;i++) cout << rate[i] << ",";
cout << endl;

TFile *outputfile=new TFile("output.root", "recreate");

gStyle->SetPaintTextFormat("2.0f");

TCanvas *c_rate_thetaEf = new TCanvas("rate_ThetaP","rate_ThetaP",1200,900);
gPad->SetLogz();
// T->Draw("(g2_theta*180/3.1415926):g2_p>>h2(11000,0,11,180,0,180)","rate","colz");
// T->Draw("pf:(theta*180/3.1415926)>>h2(20,20,40,22,0,11)","","colz text");
T->Draw("pf:(theta*180/3.1415926)>>h2(20,20,40,22,0,11)","rate","colz text");
// T->Draw("theta:Ef>>h2(11,0,11,10,0,3.14)","rate","TEXT");
// T->Draw("theta:Ef>>h2(11,0,11,10,0,3.14)","rate*(W>2 && Q2 >1)","TEXT");
// T->Draw("(theta*180/3.1415926):pf>>h2(11000,0,11,180,0,180)","rate*(W>2 && Q2 >1)","colz");
// T->Draw("(theta*180/3.1415926):pf>>h2(11000,0,11,180,0,180)","rate","colz");
// T->Draw("(theta*180/3.1415926):pf>>h2(11000,0,11,180,0,180)","rate*(9<=theta*180/3.1415926 && theta*180/3.1415926<15)","colz");
// T->Draw("(theta*180/3.1415926):pf>>h2(11000,0,11,180,0,180)","rate*(15<=theta*180/3.1415926 && theta*180/3.1415926<25)","colz");
// cout << h2->GetSum() << endl;
// T->Draw("theta:Ef>>h2(1,1.337,1.463,1,0.74,0.83)","rate*(0<phi&&phi<0.092)","TEXT"); 
// 42.4/180*3.1415926=0.74
// 47.6/180*3.1415926=0.83
// ./main 1 1000. 5.7 1 1.337 1.463 42.4 47.6 SoLID

outputfile->Write();
outputfile->Flush();

}