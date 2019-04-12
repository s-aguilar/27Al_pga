#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TError.h"

// Fitting routines
#include "calibration/fitFunctions.h"


// if 1, use local, if 0 use CRC
int loc = 0;


void peakFitter(const char *fileName, const char *detector,
	int detLoop){


	// Reset global variables
	gROOT->Reset();

	// Get root file
	TFile *fyield = new TFile(fileName);

	// Get integrated charge from root file
	TVectorD *qcharge = static_cast<TVectorD*>(fyield->Get("pulses"));
	double charge = (*qcharge)[0];

	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));

	gStyle->SetOptFit(1111);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,2);
	c0->Update();

	c0->cd(1);

	// // Prepare histograms and name them
	// TH1D *h1 = new TH1D(Form("h1 - det%i",detLoop),Form("h1 - det%i",detLoop),4096,0,4096);
	// TH1D *h2 = new TH1D(Form("h2 - det%i",detLoop),Form("h2 - ECAL - det%i",detLoop),4096,0,4096);

	double area;
	double area_err;
	double chi2NDF;
	double sig1;
	double sig2;
	bool isValid;
	int status;

	// Draw results
	hyield->Draw();
	hyield->GetXaxis()->SetRangeUser(600,1350);
	hyield->SetStats(kFALSE);
	c0->cd(2);
	hyield->GetXaxis()->SetRangeUser(600,1350);
	hyield->Draw();
	hyield->SetStats(kFALSE);


	// ATTEMPT TO FIT
	try {
		// .Get() returns the contained pointer to TFitResult. Dereference it with "*""
		TFitResult fitResults = static_cast<TFitResult>(*single_gauss_area_p2(hyield).Get());
		// fitResults.Print();

		area = fitResults.Parameter(3);
		area_err = fitResults.Error(3);
		chi2NDF = fitResults.Chi2()/fitResults.Ndf();
		sig1 = fitResults.Parameter(5);
		isValid = fitResults.IsValid();
		status = fitResults.Status();


		// The charge is integrated charge of proton which is 1
		double yield = area/(charge);
		double yield_err = area_err/(charge);

		string runNum = fileName;
		if (loc==1) runNum = runNum.substr(8,4);
		else runNum = runNum.substr(70,4);

		string detNum = detector;

		c0->SaveAs(Form("Yields/P2/run_%s/%s_Fit.png",runNum.c_str(),detNum.c_str()));
		c0->SaveAs(Form("Yields/P2/det-%i/run_%s_Fit.png",detLoop,runNum.c_str()));

		ofstream myfile;
		myfile.open ("Yields/P2/_P2.csv",std::ios::app);
		myfile<<Form("run_%s",runNum.c_str())<<","<< Form("%s",detNum.c_str())
			<<","<<yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<sig1
			<<","<<chi2NDF<<","<<isValid<<","<<status<<","<<charge<<"\n";
		myfile.close();
	}catch(...){}

	c0->Clear();
	fyield->Close();
	delete c0;
}



/*==============================MAIN=========================================*/
void p2Yields_pa(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/27Al_pa/";
	chdir(path);

	const char *detect;
	const char *files;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("Yields/P2/_P2.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","
		<<"Area"<<","<<"Area err"<<","<<"sig1"<<","<<"X2NDF"<<","<<"IsValid"<<","
		<<"Status"<<","<<"Q_int"<<"\n";
	myfile.close();

	// try {
	// 	gSystem->Exec(Form("mkdir Yields"));
	// 	gSystem->Exec(Form("mkdir Yields/A1"));
	// 	gSystem->Exec(Form("mkdir Yields/P1"));
	// 	gSystem->Exec(Form("mkdir Yields/P2"));
	// }catch(...){}
	// try {
	// 	gSystem->Exec(Form("mkdir TFitResult"));
	// 	gSystem->Exec(Form("mkdir TFitResult/A1"));
	// 	gSystem->Exec(Form("mkdir TFitResult/P1"));
	// 	gSystem->Exec(Form("mkdir TFitResult/P2"));
	// }catch(...){}
	// // Make directory to visually inspect the fits
	// for(int ii = 0; ii<13; ii++){
	// 	try {
	// 		gSystem->Exec(Form("mkdir Yields/P2/det-%i",ii));
	// 		gSystem->Exec(Form("mkdir TFitResult/P2/det-%i",ii));
	// 	}catch(...){}
	// }


	// Loop through runs: 97-107
	cout << "\nBEGINNING PEAK FITTING:" << endl;

	int fileNum = 1;
	int runStart = 97 ;
	int upToRun;

	if (loc==1) upToRun = 108;	// 108
	else upToRun = 1177;		// 1177

	for(int i=runStart;i<upToRun;i++){

		// Skip bad runs
		if(i==98) continue;
		else if(i==102) continue;
		else if(i==123) continue;
		else if(i==124) continue;
		else if(i==146) continue;
		else if(i==152) continue;
		else if(i==160) continue;
		else if(i>=167 && i<=179) continue;
		else if(i==221) continue;
		else if(i==222) continue;
		else if(i==223) continue;
		else if(i==299) continue;
		else if(i==308) continue;
		else if(i==331) continue;
		else if(i==360) continue;
		else if(i==385) continue;
		else if(i==462) continue;
		else if(i==479) continue;
		else if(i==488) continue;
		else if(i==506) continue;
		else if(i==576) continue;
		else if(i==590) continue;
		else if(i==599) continue;
		else if(i==609) continue;
		else if(i==615) continue;
		else if(i==628) continue;
		else if(i==633) continue;
		else if(i==636) continue;
		else if(i==640) continue;
		else if(i==645) continue;
		else if(i==649) continue;
		else if(i==670) continue;
		else if(i==672) continue;
		else if(i==679) continue;
		else if(i==683) continue;
		else if(i==685) continue;
		else if(i==686) continue;
		else if(i==690) continue;
		else if(i==727) continue;
		// else if(i==863) continue; // BG
		// else if(i==887) continue;
		// else if(i==922) continue;
		// else if(i==939) continue; // BG
		// else if(i==942) continue;
		else if(i>=808 && i<=952) continue;	// TUNE PROBLEM IN P1 channel
		else if(i>=953 && i<=959) continue;
		else if(i==980) continue;
		else if(i==984) continue;
		else if(i==997) continue;	// Rescan of tune issue, 1 point jumps up
		else if(i==1003) continue;
		else if(i==1148) continue;
		else if(i==1149) continue;
		else if(i==1151) continue;
		else if(i==1152) continue;
		else if(i==1153) continue;
		else if(i==1164) continue;

		TString runNum_TString;

		if(i < 100) runNum_TString = "00";
		else if((i >= 100) && (i < 1000)) runNum_TString = "0";
		else runNum_TString = "";

		runNum_TString += i;	// Should be format 0001 -> 9999
		const char *runNum_String = (const char*)runNum_TString;

		// try {
		// 	gSystem->Exec(Form("mkdir Yields/P2/run_%s",runNum_String));
		// }catch(...){}

		// // Want to open TFILE in subdirectory and save fit results evetually.
		// TFile *fitResults = new TFile(Form(""));

		// Loop through detectors on board
		for(int j=0;j<13;j++){ // 13

			if (loc==1) files = Form("run/run_%s.root",runNum_String);
			else files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/27Al_p_March_2019/run/run_%s.root",runNum_String);

			detect = Form("det%d",j);

			// Perform peak fitting
			peakFitter(files,detect,j);
		}

		cout << Form("Fitting run_%s complete",runNum_String) << endl;
		fileNum+=1;
	}
}
