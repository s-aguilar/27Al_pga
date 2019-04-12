#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TVectorD.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"


void extractSpectra(){

    int local = 0;

    cout << "\nBeginning Processing:" << endl;
    // Loop through runs, creating the histograms
    for (int runNum=1; runNum<=10; runNum++){

        try {
            TString data_path, out_path;

            if (local == 1) {
                data_path = "/Users/sebastian/Desktop/27Al_pga/rawRun/";
                out_path = "/Users/sebastian/Desktop/27Al_pga/run/";
            }
            else {
                data_path = "/afs/crc/group/nsl/activetarget/data/27Al_p_March_2019/root/";
                out_path = "/afs/crc/group/nsl/activetarget/data/27Al_p_March_2019/run/";
            }

            TString filename, outname;
            filename += data_path;
            filename += "run_";
            outname += "run_";
            cout << runNum << endl;
            if(runNum < 10){
                filename += "00";
                outname += "000"; // "00"
            }
            else if((runNum >= 10) && (runNum < 100)){
                filename += "0";
                outname += "00";
            }
            else if((runNum >= 100) && (runNum < 1000)){
                outname += "0";
            }

            cout << outname << endl;

            filename += runNum;
            outname += runNum;
            filename += ".root";
            outname += ".root";



            cout << outname << " is being processed." << endl;

            // Input file
            TFile *f = new TFile(filename.Data());

            // Output file
            TFile *f1 = new TFile(out_path+outname,"RECREATE");

            // TTree *T = new TTree("T","Data");

            // histogram definitions
            TH1D *hist[13];

            for(Int_t det_i = 0; det_i < 13; det_i++){
                TString hist_name[13];
                hist_name[det_i] += TString("det");
                hist_name[det_i] += to_string(det_i);
                hist[det_i] = new TH1D(hist_name[det_i].Data(),"ADC spectrum",4096,0.0,4095.0);
                // hist[det_i]->Write();
                // T->Branch(TString::Format("hist[%i]",det_i),"TH1D",&hist[det_i],32000,0);
            }


            // tree setup
            TTreeReader reader("evtTree",f);
            TTreeReaderArray<UShort_t> adc(reader, "event.l");
            TTreeReaderArray<UShort_t> eventTime(reader, "event.t");

            TTreeReader scaler_reader("scalerTree",f);
            TTreeReaderArray<UInt_t> scaler(scaler_reader, "scaler.sc");

            int iEvt = 0;
            while(reader.Next()){
                for(int det_i=0; det_i<13; det_i++){
                    hist[det_i]->Fill(adc[det_i]);
                    // hist[det_i]->GetYaxis()->SetRangeUser(0,1.2*hist[det_i]->GetMaximum())
                }
                iEvt++;
            }

            double qlive = 0.0;
            while(scaler_reader.Next()){
                qlive = scaler[3];
            }


            TVectorD pulses(1);
            pulses[0] = qlive;      // 10^-8 C/pulse
            pulses.Write(TString("pulses"));

            // TVectorD lastTimeStamp(1);
            // lastTimeStamp[0] = 100;      ///// THIS NEEDS TO BE FURTHER DEVELOPED
            // lastTimeStamp.Write(TString("lastTimeStampSeconds"));


            f1->Write();    // Writes all objects in memory to file

            f1->Close();
            delete f;
            delete f1;

            gROOT->Reset();

        }catch(...){}
    }

}
