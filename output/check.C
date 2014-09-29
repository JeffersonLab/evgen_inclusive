#include<iostream>
#include<fstream>
using namespace std;


void check(){

  gStyle->SetOptStat(11111111);
//	TFile *file_1=new TFile("../zhiwen_target/PVDIS_LD2/rate_solid_PVDIS_LD2_pip_5e7.root"); //location of zhiwen's rate tree from target
// TFile *file_1=new TFile("../../eicRate_20101102/output/PVDIS_LD2/rate_solid_PVDIS_LD2_pim_1e6.root");
// TFile *file_1=new TFile("../../eicRate_20101102/output/SIDIS_He3/rate_solid_SIDIS_He3_Km_1e6.root");
  //TFile *file_1=new TFile("../../eicRate_20101102/output/output.root");
// 	TFile *file_1=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/JPsi_LH2/rate_solid_JPsi_LH2_pi0_1e6.root");  
	TFile *file_1=new TFile("solid_SIDIS_He3_pip.root");  	
// 	TFile *file_1=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/b.root");  			
// 	TFile *file_1=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/SIDIS_He3/rate_solid_SIDIS_He3_Kl_1e6.root");
// 	TFile *file_1=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/SIDIS_He3/rate_solid_SIDIS_He3_window_upstream_Kl_1e6.root");	

	
//       TFile *file_2=new TFile("PVDIS_LD2/rate_solid_PVDIS_LD2_pim_1e6.root");
// 	TFile *file_2=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/a.root");  		
	TFile *file_2=new TFile("solid_SIDIS_He3_pip_new.root");  		
	
// 	TFile *file_2=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_20101102/output/JPsi_LH2_noweight/rate_solid_JPsi_LH2_pi0_1e6.root");  	
// 	      TFile *file_2=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_new/output/SIDIS_He3/rate_solid_SIDIS_He3_Kl_1e6.root");
//       TFile *file_2=new TFile("/home/zwzhao/work_halla/solid/solid_svn/solid/evgen/eicRate_new/output/SIDIS_He3/rate_solid_SIDIS_He3_window_upstream_Kl_1e6.root");
// 	      TFile *file_2=new TFile("SIDIS_He3/output.root");
	
	TTree *T_1=(TTree*)file_1->Get("T");	
	TTree *T_2=(TTree*)file_2->Get("T");

	double rate_1;
	double vz_1;
	double pf_1;
	double theta_1;

	T_1->SetBranchAddress("rate",&rate_1);
	T_1->SetBranchAddress("vz",&vz_1);
	T_1->SetBranchAddress("pf",&pf_1);
	T_1->SetBranchAddress("theta",&theta_1);

	double rate_2;
	double vz_2;
	double pf_2;
	double theta_2;

	T_2->SetBranchAddress("rate",&rate_2);	
	T_2->SetBranchAddress("vz",&vz_2);
	T_2->SetBranchAddress("pf",&pf_2);
	T_2->SetBranchAddress("theta",&theta_2);

	const int bin_num=100;
// 	const double min=-3.5;
// 	const double max=-3.0;	
// 	const double min=-0.4;
// 	const double max=0.4;
	const double min=-3.8;
	const double max=-3.2;
// 	const double min=-3.30010;
// 	const double max=-3.29990;
// 	const double min=-3.70010;
// 	const double max=-3.69990;	
	
// 	TH1F *h_1=new TH1F("h_1","h_1",bin_num,min,max);
// 	TH1F *h_2=new TH1F("h_2","h_2",bin_num,min,max);
// 	TH1F *h_ratio=new TH1F("h_ratio","h_ratio of h_1/h_2",bin_num,min,max);
// 	
// 	T_1->Project("h_1","vz","rate");
// 	T_2->Project("h_2","vz","rate");	
// 	h_ratio->Divide(h_1,h_2);
// 
// 	TCanvas *can=new TCanvas("can","can",1200,1000);
// 	can->Divide(1,3);
// 	can->cd(1);
// 	h_1->Draw();
// 	can->cd(2);
// 	h_2->Draw();
// 	can->cd(3);
// 	h_ratio->Draw();

	TH2F *h_1=new TH2F("h_1","h_1",110,0,11,35,0,3.5);
	TH2F *h_2=new TH2F("h_2","h_2",110,0,11,35,0,3.5);
	TH2F *h_ratio=new TH2F("h_ratio","h_ratio of h_1/h_2",110,0,11,35,0,3.5);
	T_1->Project("h_1","theta:pf","rate");
	T_2->Project("h_2","theta:pf","rate");	
// 	T_1->Project("h_1","theta:pf");
// 	T_2->Project("h_2","theta:pf");		
	h_ratio->Divide(h_1,h_2);	
// 	h_ratio->Add(h_1,h_2,1,-1);		

	TCanvas *can=new TCanvas("can","can",1200,1000);
	can->Divide(1,3);
	can->cd(1);
	h_1->Draw("colz");
	can->cd(2);
	h_2->Draw("colz");
	can->cd(3);
	h_ratio->Draw("colz");

}



