{  
//   gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

const int m=10;

// // char* input_filename[m]={
// // "rate_solid_SIDIS_3he_eDIS_1e6.root",
// // "rate_solid_SIDIS_3he_eES_1e6.root",
// // "rate_solid_SIDIS_3he_pip_1e6.root",
// // "rate_solid_SIDIS_3he_pim_1e6.root",
// // "rate_solid_SIDIS_3he_pi0_1e6.root",
// // "rate_solid_SIDIS_3he_Kp_1e6.root",
// // "rate_solid_SIDIS_3he_Km_1e6.root",
// // "rate_solid_SIDIS_3he_Ks_1e6.root",
// // "rate_solid_SIDIS_3he_Kl_1e6.root",
// // "rate_solid_SIDIS_3he_p_1e6.root"
// // };
// // double xmin=5,xmax=30;
// // double thetamin=7.5,thetamax=24;

char* input_filename[m]={
"rate_solid_PVDIS_LD2_eDIS_1e6.root",
"rate_solid_PVDIS_LD2_eES_1e6.root",
"rate_solid_PVDIS_LD2_pip_1e6.root",
"rate_solid_PVDIS_LD2_pim_1e6.root",
"rate_solid_PVDIS_LD2_pi0_1e6.root",
"rate_solid_PVDIS_LD2_Kp_1e6.root",
"rate_solid_PVDIS_LD2_Km_1e6.root",
"rate_solid_PVDIS_LD2_Ks_1e6.root",
"rate_solid_PVDIS_LD2_Kl_1e6.root",
"rate_solid_PVDIS_LD2_p_1e6.root"
};
// char* input_filename[m]={
// "rate_solid_PVDIS_LD2_eDIS_5e7.root",
// "rate_solid_PVDIS_LD2_eDIS_5e7.root",
// "rate_solid_PVDIS_LD2_pip_5e7.root",
// "rate_solid_PVDIS_LD2_pim_5e7.root",
// "rate_solid_PVDIS_LD2_pi0_5e7.root",
// "rate_solid_PVDIS_LD2_Kp_5e7.root",
// "rate_solid_PVDIS_LD2_Km_5e7.root",
// "rate_solid_PVDIS_LD2_Ks_5e7.root",
// "rate_solid_PVDIS_LD2_Kl_5e7.root",
// "rate_solid_PVDIS_LD2_p_5e7.root"
// };
double xmin=20,xmax=40;
double thetamin=21,thetamax=36;

int color[m]={1,2,3,4,5,6,7,8,9,46};
// int pid[m]={1,2,2,5,6,1,2,7,8,5,6,1,4};

char *label[m]={"e(DIS)","e(ES)","#pi^{+}","#pi^{-}","#pi^{0}","K^{+}","K^{-}","K_{s}","K_{l}","p"};

TCanvas *c_ThetaP = new TCanvas("ThetaP","ThetaP",1800,900);
c_ThetaP->Divide(5,2);

TCanvas *c_ThetaP_1D_nocut = new TCanvas("ThetaP_1D_nocut","ThetaP_1D_nocut",1800,900);
c_ThetaP_1D_nocut->Divide(2,1);

TCanvas *c_ThetaP_1D_cut = new TCanvas("ThetaP_1D_cut","ThetaP_1D_cut",1800,900);
c_ThetaP_1D_cut->Divide(2,1);

TH1F *hTheta[2][m],*hP[2][m];
TH2F *hThetaP[m];

double rate_total[2][m][2];

for(int i=0;i<m;i++){

  TFile *file=new TFile(input_filename[i]);
  if (file->IsZombie()) {
       cout << "Error opening file" << input_filename[i] << endl;
       exit(-1);
    }
    else cout << "open file " << input_filename[i] << endl;
      

TFile *outputfile=new TFile("output.root", "recreate");

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
char weight[100];

c_ThetaP_1D_nocut->cd(1);
gPad->SetLogy();
sprintf(content,"(theta*180/3.1415926)>>htemp1(180,0,180)");
if (i==0) sprintf(weight,"rate*(W>2)");
else sprintf(weight,"rate");
T->Draw(content,weight,"goff");
hTheta[0][i] = (TH1F*)htemp1->Clone();
hTheta[0][i]->SetLineColor(color[i]);
rate_total[0][i][0]= hTheta[0][i]->Integral();
if (i==0) hTheta[0][i]->Draw();
else  hTheta[0][i]->Draw("same");

c_ThetaP_1D_nocut->cd(2);
gPad->SetLogy();
sprintf(content,"pf>>htemp2(220,0,11)");
if (i==0) sprintf(weight,"rate*(W>2)");
else sprintf(weight,"rate");
T->Draw(content,weight,"goff");
hP[0][i] = (TH1F*)htemp2->Clone();
hP[0][i]->SetLineColor(color[i]);
rate_total[0][i][1]= hP[0][i]->Integral();
if (i==0) hP[0][i]->Draw();
else  hP[0][i]->Draw("same");

c_ThetaP_1D_cut->cd(1);
gPad->SetLogy();
sprintf(content,"(theta*180/3.1415926)>>htemp3(180,0,180)");
if (m==0) sprintf(weight,"rate*(W>2)*(pf>2)");
else sprintf(weight,"rate*(pf>2)");
T->Draw(content,weight,"goff");
hTheta[1][i] = (TH1F*)htemp3->Clone();
hTheta[1][i]->SetLineColor(color[i]);
rate_total[1][i][0]= hTheta[1][i]->Integral();
if (i==0) hTheta[1][i]->Draw();
else  hTheta[1][i]->Draw("same");

c_ThetaP_1D_cut->cd(2);
gPad->SetLogy();
sprintf(content,"pf>>htemp4(220,0,11)");
if (m==0) sprintf(weight,"rate*(W>2 && %i<theta*180/3.14159 && theta*180/3.14159<%i)",thetamin,thetamax);
else sprintf(weight,"rate*(%i<theta*180/3.14159 && theta*180/3.14159<%i)",thetamin,thetamax);
T->Draw(content,weight,"goff");
hP[1][i] = (TH1F*)htemp4->Clone();
hP[1][i]->SetLineColor(color[i]);
rate_total[1][i][1]= hP[1][i]->Integral();
if (i==0) hP[1][i]->Draw();
else  hP[1][i]->Draw("same");

c_ThetaP->cd(i+1);
gPad->SetRightMargin(2);
gPad->SetLogz();
sprintf(content,"pf:(theta*180/3.1415926)>>htemp5(360,0,180,220,0,11)");
// sprintf(content,"pf:(theta*180/3.1415926)>>htemp3(180,0,180,22,0,11)");
if (i==0) {sprintf(weight,"rate*(W>2)");}
else {sprintf(weight,"rate");}
T->Draw(content,weight,"goff");
hThetaP[i] = (TH2F*)htemp5->Clone();
sprintf(hstname,"ThetaP %s",label[i]);
hThetaP[i]->SetTitle(hstname);
hThetaP[i]->Draw("colz");
}

gStyle->SetPaintTextFormat("2.0f");

c_ThetaP_1D_nocut->cd(1);
TPaveText *pt1 = new TPaveText(0.75,0.6,0.95,0.9,"brNDC");
pt1->SetFillColor(17);
pt1->SetTextAlign(12);
pt1->Draw();
for(int i=0;i<m;i++){
char textcontent[100];
sprintf(textcontent,"%s \t\t %.2g Hz",label[i],rate_total[0][i][0]);
TText *text1=pt1->AddText(textcontent);
text1->SetTextColor(color[i]);
text1->SetTextSize(0.024);
}

c_ThetaP_1D_nocut->cd(2);
TPaveText *pt2 = new TPaveText(0.75,0.6,0.95,0.9,"brNDC");
pt2->SetFillColor(17);
pt2->SetTextAlign(12);
pt2->Draw();
for(int i=0;i<m;i++){
char textcontent[100];
sprintf(textcontent,"%s \t\t %.2g Hz",label[i],rate_total[0][i][1]);
TText *text2=pt2->AddText(textcontent);
text2->SetTextColor(color[i]);
text2->SetTextSize(0.024);
}

c_ThetaP_1D_cut->cd(1);
TPaveText *pt3 = new TPaveText(0.75,0.6,0.95,0.9,"brNDC");
pt3->SetFillColor(17);
pt3->SetTextAlign(12);
pt3->Draw();
for(int i=0;i<m;i++){
char textcontent[100];
sprintf(textcontent,"%s \t\t %.2g",label[i],rate_total[1][i][0]);
TText *text3=pt3->AddText(textcontent);
text3->SetTextColor(color[i]);
text3->SetTextSize(0.024);
}

c_ThetaP_1D_cut->cd(2);
TPaveText *pt4 = new TPaveText(0.75,0.6,0.95,0.9,"brNDC");
pt4->SetFillColor(17);
pt4->SetTextAlign(12);
pt4->Draw();
for(int i=0;i<m;i++){
char textcontent[100];
sprintf(textcontent,"%s \t\t %.2g",label[i],rate_total[1][i][1]);
TText *text4=pt4->AddText(textcontent);
text4->SetTextColor(color[i]);
text4->SetTextSize(0.024);
}

TH2F *he=hThetaP[0].Clone();
he->Rebin2D(2,10);
TH2F *hpip=hThetaP[2]->Clone();
hpip->Rebin2D(2,10);
TH2F *hpim=hThetaP[3]->Clone();
hpim->Rebin2D(2,10);
TH2F *hpi=hpip->Clone();
hpi->Add(hpip,hpim);
TH2F *hpie=hpi->Clone();
hpie->Divide(he);
TH2F *hpime=hpim->Clone();
hpime->Divide(he);

TCanvas *c_ThetaP_eDIS = new TCanvas("ThetaP_eDIS","ThetaP_eDIS",1800,900);
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
TCanvas *c_ThetaP_pie_ratio = new TCanvas("ThetaP_pie_ratio","ThetaP_pie_ratio",1800,900);
gPad->SetLogz();
hpie->SetTitle("ThetaP #pi^{+} + #pi^{-} / e(DIS) ratio");
hpie->SetAxisRange(xmin,xmax,"X");
hpie->Draw("colz text");
TCanvas *c_ThetaP_pime_ratio = new TCanvas("ThetaP_pime_ratio","ThetaP_pime_ratio",1800,900);
gPad->SetLogz();
hpime->SetTitle("ThetaP #pi^{-} / e(DIS) ratio");
hpime->SetAxisRange(xmin,xmax,"X");
hpime->Draw("colz text");

outputfile->Write();
outputfile->Flush();

}