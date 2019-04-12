#ifndef FITFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FITFUNCTIONS_H

# include <iostream>

#include <vector>
using std::vector;

#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLine.h"


double background_func(double *x, double *par){
	/*
	COEFFICIENTS:
	===============
	constant = par[0]
	linear = par[1]
	quadratic = par[2]
	*/

	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}


double single_gauss_func(double *x, double *par){
	/*
	norm = par[0]
	mean = par[1]
	sigma = par[2]
	*/

	return par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
}


double double_gauss_func(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma1 = par[2]
	norm2 = par[3]
	mean2 = par[4]
	sigma2 = par[5]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[4],par[5],kTRUE);
	return g1+g2;
}


double double_gauss_same_width_func(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma = par[2]
	norm2 = par[3]
	mean2 = par[4]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[4],par[2],kTRUE); // par[5]->par[2]
	return g1+g2;
}


double double_gauss_same_width_func_p1(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma = par[2]
	norm2 = par[3]
	mean2 = par[4]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[1]+24,par[2],kTRUE); // par[4]->par[1]+24 , par[5]->par[2]
	return g1+g2;
}


double triple_gauss_same_width_func(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma = par[2]

	norm2 = par[3]

	norm3 = par[4]
	mean3 = par[5]
	sigma3 = par[6]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[5]-25,par[6],kTRUE);
	double g3 = par[4]*TMath::Gaus(x[0],par[5],par[6],kTRUE);
	return g1+g2+g3;
}


double fit_single_gauss_func(double *x, double *par){
    // Fit a Gaussian with a Quadratic bakground

	return background_func(x,par) + single_gauss_func(x,&par[3]);
}


double fit_double_gauss_func(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_func(x,&par[3]);
}


double fit_double_gauss_same_width_func(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_same_width_func(x,&par[3]);
}


double fit_double_gauss_same_width_func_p1(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_same_width_func_p1(x,&par[3]);
}


double fit_triple_gauss_same_width_func(double *x, double *par) {
    // Fit a triple Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + triple_gauss_same_width_func(x,&par[3]);
}


/*===========================================================================*/

TFitResultPtr triple_gauss_area(TH1D *H0){

	double low = 1090;
	double high = 1300;
	double a1PeakEst = 1135;

	TF1 *ffit1 = new TF1("ffit1",fit_triple_gauss_same_width_func,low,high,10);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","norm3","mean3","sigma3");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,a1PeakEst,10,2000,400,a1PeakEst+80,10);

	// Turn off linear background terms, fix peak positions
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(4,a1PeakEst-15,a1PeakEst+15);
	ffit1->SetParLimits(5,7,17);	// Std dev range
	ffit1->SetParLimits(6,0,1e6);	// Normalization
	ffit1->SetParLimits(7,0,1e6);	// Normalization
	ffit1->SetParLimits(8,a1PeakEst+75,a1PeakEst+95);
	ffit1->SetParLimits(9,7,20);	// Std dev range
	TFitResultPtr r = H0->Fit("ffit1","SQR");
	return r;

}

TFitResultPtr single_gauss_area_p1(TH1D *H0){

	double low = 663; //660
	double high = 750; //760

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetParameters(1,0,0,2000,703,12);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	TFitResultPtr r = H0->Fit("ffit1","SQR");
	return r;

}

TFitResultPtr single_gauss_area_p2(TH1D *H0){

	double low = 785;	// 780
	double high = 920;	// 900

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
    ffit1->SetParameters(1,0,0,2000,850,12); //840
    ffit1->FixParameter(2,0);		// Makes it a linear background

	TFitResultPtr r = H0->Fit("ffit1","SQR");
	return r;

}

TFitResultPtr single_gauss_area_co1(TH1D *H0){

	double low = 930;
	double high = 1035;

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetParameters(1,0,0,2000,965,12);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	TFitResultPtr r = H0->Fit("ffit1","SQR");
	TLine *line = new TLine(r->Parameter(4),H0->GetMinimum(),
							r->Parameter(4),1.1*H0->GetMaximum());
	line->SetLineColor(2);
	line->Draw("SAME");
	return r;

}

TFitResultPtr single_gauss_area_co2(TH1D *H0){

	double low = 1030;
	double high = 1165;

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetParameters(1,0,0,2000,1100,12);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	TFitResultPtr r = H0->Fit("ffit1","SQR");
	TLine *line = new TLine(r->Parameter(4),H0->GetMinimum(),
							r->Parameter(4),1.1*H0->GetMaximum());
	line->SetLineColor(2);
	line->Draw("SAME");
	return r;

}

TFitResultPtr single_gauss_area_cs(TH1D *H0){

	double low = 500;
	double high = 600;

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetParameters(1,0,0,2000,550,12);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	TFitResultPtr r = H0->Fit("ffit1","SQR");
	TLine *line = new TLine(r->Parameter(4),H0->GetMinimum(),
							r->Parameter(4),1.1*H0->GetMaximum());
	line->SetLineColor(2);
	line->Draw("SAME");
	return r;

}


/*===========================================================================*/
#endif
