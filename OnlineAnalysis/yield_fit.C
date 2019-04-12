#include <iostream>
using namespace std;

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
#include "fitFunctions.h"

Int_t yield_fit(){

  TString data_path("/afs/crc/group/nsl/activetarget/data/27Al_p_March_2019/root/");
  Int_t runNum;


  cout << "Enter run number: ";
  cin >> runNum;

  TString filename;
  filename += data_path;
  filename += "run_";
  if(runNum < 10){
    filename += "00";
  }
  else if((runNum >= 10) && (runNum < 100)){
    filename += "0";
  }

  filename += runNum;
  filename += ".root";

  //cout << filename << endl;

  TFile *f = new TFile(filename.Data());
  if(f->IsOpen()){
    cout << "File " << filename.Data() << " is open" << endl;
  }
  else{
    cout << "Could not open file " << filename.Data() << endl;
    return -1;
  }

  // histogram definitions
  TH1D *hist[13];

  for(Int_t det_i = 0; det_i < 13; det_i++){
    TString hist_name[13];
    hist_name[det_i] += TString("det");
    hist_name[det_i] += to_string(det_i);
    hist[det_i] = new TH1D(hist_name[det_i].Data(),"ADC spectrum",4096,0.0,4095.0);
  }

  // tree setup
  TTreeReader reader("evtTree",f);
  TTreeReaderArray<UShort_t> adc(reader, "event.l");

  TTreeReader scaler_reader("scalerTree",f);
  TTreeReaderArray<UInt_t> scaler(scaler_reader, "scaler.sc");
  
  Int_t iEvt = 0;
  while(reader.Next()){
    if((iEvt%10000) == 0){
      cout << "Event " << iEvt << endl;
    }
    for(Int_t det_i=0; det_i<13; det_i++){
      //cout << "detector " << iDet << ": " << adc[det_adc[iDet]] << "\t" << tdc[det_tdc[iDet]] << "\t" << endl;
      hist[det_i]->Fill(adc[det_i]);
    }
    iEvt++;
  }

  // fit histogram

  Int_t detNum;
  detNum=0;
  // cout << "Enter detector no: ";
  // cin >> detNum;
  if((detNum < 0) || (detNum > 12)){
    cout << "Detector number must be between 0 and 12!" << endl;
    return -1;
  }

  TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,1400,1525,6);
  ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
  ffit1->SetParameters(1,0,0,2000,1475,12);
  ffit1->FixParameter(2,0);		// Makes it a linear background
  
  hist[detNum]->GetXaxis()->SetRangeUser(1400,1525);
  hist[detNum]->Fit("ffit1","R");

  cout << "Area = " << ffit1->GetParameter(3) << " +- " << ffit1->GetParError(3) << endl;

  Double_t qlive = 0.0;
  while(scaler_reader.Next()){
    qlive = scaler[3];
    //cout << "qlive: " << qlive << endl;
  }

  cout << endl << "Q.live = " << qlive << endl;

  return 0;
}
