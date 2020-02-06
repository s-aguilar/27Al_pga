#include <iostream>
using std::cout;
using std::endl;

#include <chrono>  // for high_resolution_clock

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


// if true use local, false use CRC
bool loc = false;

// if true save plots, false don't save
bool plot = true;


void peakFitter(TFile *TFitOut, const char *fileName, const char *detector,
	int detLoop){

	// Get root file
	TFile *fyield = new TFile(fileName);

	// Change output directory to TFitOut
	TFitOut->cd();

	// Get integrated charge and run time from root file
	TVectorD *qcharge = static_cast<TVectorD*>(fyield->Get("pulses"));
	double charge = (*qcharge)[0];
	TVectorD *ttimes = static_cast<TVectorD*>(fyield->Get("times"));
	double times = (*ttimes)[0];

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
	// hyield->Draw();
	hyield->GetXaxis()->SetRangeUser(1300,1700);
	hyield->SetStats(kFALSE);
	c0->cd(2);
	hyield->GetXaxis()->SetRangeUser(1400,1550);
	hyield->SetStats(kFALSE);

	// ATTEMPT TO FIT
	try {
		string runNum = fileName;
		if (loc) runNum = runNum.substr(8,4);
		else runNum = runNum.substr(70,4);

		// .Get() returns the contained pointer to TFitResult. Dereference it with "*""
		TFitResult fitResults = static_cast<TFitResult>(*single_gauss_area_Thick(hyield).Get());

		// Write out fitResults to current TFile (TFitOut)
		fitResults.Write(Form("det-%i",detLoop));

		area = fitResults.Parameter(3);
		area_err = fitResults.Error(3);
		chi2NDF = fitResults.Chi2()/fitResults.Ndf();
		sig1 = fitResults.Parameter(5);
		isValid = fitResults.IsValid();	// Check if fit was succesful, minimum
										// was found, return type bool

		status = fitResults.Status();	// Return status code of error minimization
										// This is minimizer dependent!


		// The charge is integrated charge of proton which is 1
		double yield = area/(charge);
		double yield_err = area_err/(charge);

		string detNum = detector;

		if(plot){
			c0->SaveAs(Form("Yields/Thick/run_%s/%s_Fit.pdf",runNum.c_str(),detNum.c_str()));
			c0->SaveAs(Form("Yields/Thick/det-%i/run_%s_Fit.pdf",detLoop,runNum.c_str()));
		}

		ofstream myfile;
		myfile.open ("Yields/Thick/Thick.csv",std::ios::app);
		myfile<<Form("run_%s",runNum.c_str())<<","<< Form("%s",detNum.c_str())
			<<","<<yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<sig1
			<<","<<chi2NDF<<","<<isValid<<","<<status<<","<<charge<<","<<times<<"\n";
		myfile.close();
	}catch(...){}

	c0->Clear();
	fyield->Close();

	// Delete objects on the heap
	delete c0;
	c0 = nullptr;
	delete fyield;
	fyield = nullptr;
}



/*==============================MAIN=========================================*/
void thicknessAnalysis(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/27Al_pa/";
	chdir(path);

	const char *detect;
	const char *files;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("Yields/Thick/Thick.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","
		<<"Area"<<","<<"Area err"<<","<<"sig1"<<","<<"X2NDF"<<","<<"IsValid"<<","
		<<"Status"<<","<<"Q int"<<","<<"Time"<<"\n";
	myfile.close();

	try {
		// gSystem->Exec(Form("mkdir Yields"));
		// gSystem->Exec(Form("mkdir Yields/A1"));
		// gSystem->Exec(Form("mkdir Yields/P1"));
		// gSystem->Exec(Form("mkdir Yields/P2"));
		gSystem->Exec(Form("mkdir Yields/Thick"));
	}catch(...){}
	try {
		// gSystem->Exec(Form("mkdir TFitResult"));
		// gSystem->Exec(Form("mkdir TFitResult/A1"));
		// gSystem->Exec(Form("mkdir TFitResult/P1"));
		// gSystem->Exec(Form("mkdir TFitResult/P2"));
		gSystem->Exec(Form("mkdir TFitResult/Thick"));
	}catch(...){}
	// Make directory to visually inspect the fits
	for(int ii = 0; ii<13; ii++){
		try {
			gSystem->Exec(Form("mkdir Yields/Thick/det-%i",ii));
		}catch(...){}
	}


	// Loop through runs: 97-107
	cout << "\nBEGINNING PEAK FITTING:" << endl;

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	int fileNum = 1;
	int runStart = 77 ;
	int upToRun = 97;

	for(int i=runStart;i<upToRun;i++){

		TString runNum_TString;

		if(i < 100) runNum_TString = "00";
		else if((i >= 100) && (i < 1000)) runNum_TString = "0";
		else runNum_TString = "";

		runNum_TString += i;	// Should be format 0001 -> 9999
		const char *runNum_String = (const char*)runNum_TString;

		try {
			gSystem->Exec(Form("mkdir Yields/Thick/run_%s",runNum_String));
		}catch(...){}

		// Save TFitResult results.
		TFile *TFitFiles = new TFile(Form("TFitResult/Thick/run_%s.root",runNum_String),"RECREATE");

		// Loop through detectors on board
		for(int j=0;j<13;j++){

			files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/27Al_p_March_2019/run/run_%s.root",runNum_String);

			detect = Form("det%d",j);

			// Perform peak fitting
			peakFitter(TFitFiles,files,detect,j);
		}

		// Write all TObjects in memory (TFitResult) to TFile
		TFitFiles->Write();
		TFitFiles->Close();

		// Delete objects on the heap
		delete TFitFiles;
		TFitFiles = nullptr;

		cout << Form("Fitting run_%s complete",runNum_String) << endl;
		fileNum+=1;
	}

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	cout << "Elapsed time: " << elapsed.count() << " s\n";
}
