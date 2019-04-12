#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TLine.h"

// Fitting routines
#include "fitFunctions.h"


void calibration(const char *fileName,const char *fileBack,const char *detector,
	int fileNameLoop, int detLoop){


	gROOT->Reset();

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile(fileBack);    // bg spectra

	double runTime0;
	double backTime0 = 58136;	// Length of run_1196 (BG)
	// double backTime0 = 19682;	// Length of run_1197 (BG)

	// if(fileNameLoop <= 1) runTime0 = 4735;	// Length of run_1188 (60Co)
	if(fileNameLoop <= 1) runTime0 = 2525;	// Length of run_1189 (60Co)
	else runTime0 = 68766;					// Length of run_1194 (137Cs)

	// Scale background spectra to ratio of run time
	double scale0 = runTime0/backTime0;

	// Get histograms from root file, prepare bg subtracted histogram
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	TH1D *HBACK = static_cast<TH1D*>(fbackground->Get(detector));

	HBACK->Scale(scale0);		// Scale the background
	HBACK->SetDirectory(0);
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);


	// Draw results
	hyield->Draw();
	HBACK->Draw("SAME");
	HBACK->SetLineColor(kGreen);
	hyield->SetStats(kFALSE);

	// gPad->SetLogy();

	if(fileNameLoop==0) hyield->GetXaxis()->SetRangeUser(800,1300);
	else if(fileNameLoop==1) hyield->GetXaxis()->SetRangeUser(1000,1300);
	else if(fileNameLoop==2) hyield->GetXaxis()->SetRangeUser(400,800);


	// Background subtracted spectrum
	c0->cd(2);


    // Recalculate errors manually
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = HBACK->GetBinContent(i);
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = HBACK->GetBinError(i);

		ysubtracted->SetBinContent(i,yval-yval2);
		ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
	}



	// Calculate the peak area and error
	string detNum = detector;


	ysubtracted->Draw();

	// .Get() returns the contained pointer to TFitResult. Dereference it with "*""
	TFitResult fitResults;

	if(fileNameLoop==0) {
		fitResults = static_cast<TFitResult>(*single_gauss_area_co1(ysubtracted).Get());
		ysubtracted->GetXaxis()->SetRangeUser(800,1300);
	}
	else if(fileNameLoop==1) {
		fitResults = static_cast<TFitResult>(*single_gauss_area_co2(ysubtracted).Get());
		ysubtracted->GetXaxis()->SetRangeUser(800,1300);
	}
	else if(fileNameLoop==2) {
		fitResults = static_cast<TFitResult>(*single_gauss_area_cs(ysubtracted).Get());
		ysubtracted->GetXaxis()->SetRangeUser(400,800);
	}

	double mean = fitResults.Parameter(4);
	double sigma = fitResults.Parameter(5);

	// Peak area using fit parameter
	double A = fitResults.Parameter(3);
	double A_err = fitResults.Error(3);

	ysubtracted->SetStats(kFALSE);
	// gPad->SetLogy();

	// Save fits and calibration information
	if (fileNameLoop==0) {
		c0->SaveAs(Form("calPlots/60Co_1173peak/%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1173cal.csv",std::ios::app);
		myfile<< detNum.c_str()<<","<<mean<<","<<sigma<<","<<A
				<<","<<A_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==1) {
		c0->SaveAs(Form("calPlots/60Co_1332peak/%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1332cal.csv",std::ios::app);
		myfile<< detNum.c_str()<<","<<mean<<","<<sigma<<","<<A
				<<","<<A_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==2) {
		c0->SaveAs(Form("calPlots/137Cs_661peak/%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/137Cs_661cal.csv",std::ios::app);
		myfile<< detNum.c_str()<<","<<mean<<","<<sigma<<","<<A
				<<","<<A_err<<","<<runTime0<<"\n";
		myfile.close();
	}

	c0->Clear();

	fyield->Close();
	fbackground->Close();

	delete fyield;
	delete fbackground;
	delete c0;

	delete HBACK;

}


/*============================START OF MAIN==================================*/
void effCalibration(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;


	// // Create file directory for output incase its not made
	// try {
	// 	cout << "\nATTEMPTING TO CREATE OUTPUT FILE DIRECTORIES" << endl;
	// 	gSystem->Exec(Form("mkdir csv"));
	// 	gSystem->Exec(Form("mkdir calPlots"));
	// 	gSystem->Exec(Form("mkdir calPlots/60Co_1173peak"));
	// 	gSystem->Exec(Form("mkdir calPlots/60Co_1332peak"));
	// 	gSystem->Exec(Form("mkdir calPlots/137Cs_661peak"));
	// }catch(...){}


	// List of relative root file directories and detector names\
	   These will then be fed into actual calibration function
	const char *fileLoc[] = {"60Co/run_1189/run_1189.root","60Co/run_1189/run_1189.root","137Cs/run_1194/run_1194.root"}; //1188
	const char *fileBackground = "background/run_1196/run_1196.root";
	const char *detect[] = {"det0","det1","det2","det3","det4","det5","det6","det7","det8","det9","det10","det11","det12"};

		// Outer loop: Runs over the calibration spectra
	cout << "\nBEGINNING CALIBRATION FITTING:" << endl;
	for (int i=0;i<=2;i++){

			// 1173 keV gamma
		if (i == 0){
			ofstream myfile;
			myfile.open ("csv/60Co_1173cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"A"<<','<<"A err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 1332 keV gamma
		else if (i == 1){
			ofstream myfile;
			myfile.open ("csv/60Co_1332cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"A"<<','<<"A err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 661 keV gamma
		else if (i == 2){
			ofstream myfile;
			myfile.open ("csv/137Cs_661cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"A"<<','<<"A err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// Inner loop: Runs over all 13 detectors
		for (int j=0;j<=12;j++){
			calibration(fileLoc[i],fileBackground,detect[j],i,j);
		}
		cout << Form("Files %i/3 complete",i+1) << endl;
	}

	cout << "\nDONE!" << endl;
}
