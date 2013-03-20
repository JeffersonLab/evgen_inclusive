{  
//   gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

const int m=10;

// char* input_filename[m]={
// "rate_solid_SIDIS_3he_eDIS_1e6.root",
// "rate_solid_SIDIS_3he_pip_1e6.root",
// "rate_solid_SIDIS_3he_pim_1e6.root",
// "rate_solid_SIDIS_3he_pi0_1e6.root",
// "rate_solid_SIDIS_3he_eES_1e6.root",
// "rate_solid_SIDIS_3he_Kp_1e6.root",
// "rate_solid_SIDIS_3he_Km_1e6.root",
// "rate_solid_SIDIS_3he_Ks_1e6.root",
// "rate_solid_SIDIS_3he_Kl_1e6.root",
// "rate_solid_SIDIS_3he_p_1e6.root"
// };
// double xmin=0,xmax=30;

char* input_filename[m]={
"rate_solid_PVDIS_LD2_eDIS_1e6.root",
"rate_solid_PVDIS_LD2_pip_1e6.root",
"rate_solid_PVDIS_LD2_pim_1e6.root",
"rate_solid_PVDIS_LD2_pi0_1e6.root",
"rate_solid_PVDIS_LD2_eES_1e6.root",
"rate_solid_PVDIS_LD2_Kp_1e6.root",
"rate_solid_PVDIS_LD2_Km_1e6.root",
"rate_solid_PVDIS_LD2_Ks_1e6.root",
"rate_solid_PVDIS_LD2_Kl_1e6.root",
"rate_solid_PVDIS_LD2_p_1e6.root"
};
double xmin=20,xmax=40;

int color[m]={1,2,3,4,5,6,7,8,9,46};
// int pid[m]={1,2,2,5,6,1,2,7,8,5,6,1,4};

char *label[m]={"e(DIS)","#pi^{+}","#pi^{-}","#pi^{0}","e(ES)","K^{+}","K^{-}","K_{s}","K_{l}","p"};

TCanvas *c_ThetaP = new TCanvas("ThetaP","ThetaP",1800,900);
c_ThetaP->Divide(5,2);

TCanvas *c = new TCanvas("c","c",1800,900);
c->Divide(2,1);

TH1F *hTheta[m],*hP[m];
TH2F *hThetaP[m];

for(int i=0;i<m;i++){

  TFile *file=new TFile(input_filename[i]);
  if (file->IsZombie()) {
       cout << "Error opening file" << input_filename[i] << endl;
       exit(-1);
    }
    else cout << "open file " << input_filename[i] << endl;
      

// TFile *outputfile=new TFile("output.root", "recreate");

TTree *T= (TTree*) file->Get("T"); 

char hstname[100];
// sprintf(hstname,"Theta_%i",i);
// hTheta[i]=new TH1F(hstname,hstname,30,0,30);
// hTheta[i]->SetLineColor(i+1);
// sprintf(hstname,"P_%i",i);
// hP[i]=new TH1F(hstname,hstname,220,0,11);
// hP[i]->SetLineColor(i+1);
// sprintf(hstname,"ThetaP_%s",label[i]);
// hThetaP[i]=new TH2F(hstname,hstname,30,0,30,220,0,11);

char content[100];

c->cd(1);
gPad->SetLogy();
sprintf(content,"(theta*180/3.1415926)>>htemp1(180,0,180)");
// T->Draw(content,"rate","goff");
T->Draw(content,"rate*(pf>0.5)","goff");
hTheta[i] = (TH1F*)htemp1->Clone();
hTheta[i]->SetLineColor(color[i]);
cout << hTheta[i]->Integral(8+1,24) << endl;
// cout << hTheta[i]->Integral(8+1,15) << endl;
if (i==0) hTheta[i]->Draw();
else  hTheta[i]->Draw("same");

c->cd(2);
gPad->SetLogy();
sprintf(content,"pf>>htemp2(220,0,11)");
// T->Draw(content,"rate","goff");
// T->Draw(content,"rate*(8 <theta*180/3.14159 && theta*180/3.14159 < 24)","goff");
T->Draw(content,"rate*(22 <theta*180/3.14159 && theta*180/3.14159 < 35)","goff");
hP[i] = (TH1F*)htemp2->Clone();
hP[i]->SetLineColor(color[i]);
if (i==0) hP[i]->Draw();
else  hP[i]->Draw("same");

c_ThetaP->cd(i+1);
gPad->SetRightMargin(2);
gPad->SetLogz();
sprintf(content,"pf:(theta*180/3.1415926)>>htemp3(360,0,180,220,0,11)");
// sprintf(content,"pf:(theta*180/3.1415926)>>htemp3(180,0,180,22,0,11)");
T->Draw(content,"rate","goff");
hThetaP[i] = (TH2F*)htemp3->Clone();
sprintf(hstname,"ThetaP %s",label[i]);
hThetaP[i]->SetTitle(hstname);
hThetaP[i]->Draw("colz");
}

c->cd(1);
TPaveText *pt1 = new TPaveText(0.8,0.6,0.95,0.9,"brNDC");
pt1->SetFillColor(17);
pt1->SetTextAlign(12);
pt1->Draw();
for(int i=0;i<m;i++){
TText *text=pt1->AddText(label[i]);
text->SetTextColor(color[i]);
text->SetTextSize(0.024);
}

TH2F *he=hThetaP[0].Clone();
he->Rebin2D(2,10);  
TH2F *hpip=hThetaP[1]->Clone();
hpip->Rebin2D(2,10);  
TH2F *hpim=hThetaP[2]->Clone();
hpim->Rebin2D(2,10);
TH2F *hpi=hpip->Clone();
hpi->Add(hpip,hpim);
TH2F *hpie=hpi->Clone();
hpie->Divide(he);
TH2F *hpime=hpim->Clone();
hpime->Divide(he);

gStyle->SetPaintTextFormat("2.0f");

TCanvas *c_ThetaP_e = new TCanvas("ThetaP_e","ThetaP_e",1800,900);
gPad->SetLogz();
he->SetAxisRange(xmin,xmax,"X");
he->Draw("colz text");
TCanvas *c_ThetaP_pip = new TCanvas("ThetaP_pip","ThetaP_pip",1800,900);
gPad->SetLogz();
hpip->SetAxisRange(xmin,xmax,"X");
hpip->Draw("colz text");
TCanvas *c_ThetaP_pim = new TCanvas("ThetaP_pim","ThetaP_pim",1800,900);
gPad->SetLogz();
hpim->SetAxisRange(xmin,xmax,"X");
hpim->Draw("colz text");
TCanvas *c_ThetaP_pi = new TCanvas("ThetaP_pi","ThetaP_pi",1800,900);
gPad->SetLogz();
hpi->SetTitle("ThetaP #pi^{+} + #pi^{-}");
hpi->SetAxisRange(xmin,xmax,"X");
hpi->Draw("colz text");
TCanvas *c_ThetaP_pie = new TCanvas("ThetaP_pie","ThetaP_pie",1800,900);
gPad->SetLogz();
hpie->SetTitle("ThetaP #pi^{+} + #pi^{-} / e(DIS) ratio");
hpie->SetAxisRange(xmin,xmax,"X");
hpie->Draw("colz text");
TCanvas *c_ThetaP_pime = new TCanvas("ThetaP_pime","ThetaP_pime",1800,900);
gPad->SetLogz();
hpime->SetTitle("ThetaP #pi^{-} / e(DIS) ratio");
hpime->SetAxisRange(xmin,xmax,"X");
hpime->Draw("colz text");

// outputfile->Write();
// outputfile->Flush();

}