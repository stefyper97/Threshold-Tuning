#include <iostream>
#include <fstream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TF1.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TObject.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
#define Nchip 10

using namespace std;

double row, column, th, sigma, chi2;
string prefix = "FitValues_hitmap_id0x";
string chipID[Nchip] = {"70", "71", "72", "73", "74", "78", "79", "7a", "7b", "7c"};
string suffix = ".txt";

void threshold(const char *directory){

  gSystem->cd(directory);

  TFile* hfile = new TFile("Threshold_noise.root", "recreate");

  // ******* Threshold and Noise Histo for every chip  ********

  vector<TH1D*> Thr {Nchip};
  vector<TH1D*> Noise {Nchip};
  vector<TF1*> fitThr {Nchip};
  vector<TF1*> fitNoise {Nchip};
  
  TCanvas *c = new TCanvas("Thr_distribution_HIC", "multipads", 80, 80, 1500, 1000);
  c->Divide(5, 2);

  TCanvas *cnoise = new TCanvas("Noise_distribution_HIC", "multinoise", 80, 80, 1500, 1000);
  cnoise->Divide(5, 2);
   
  for(int i = 0; i < Nchip; i++){
    string thn = "Thr" + chipID[i];
    const char* Thname = thn.c_str();
    string tit = "Threshold chip " + chipID[i];
    const char* Thtitle = tit.c_str();
    string no = "Noise" + chipID[i];
    const char* Noname = no.c_str();
    string notit = "Noise chip " + chipID[i];
    const char* Notitle = notit.c_str();
    
    Thr.at(i) = new TH1D(Thname, Thtitle, 50, 50, 350);
    Noise.at(i) = new TH1D(Noname, Notitle, 20, 0, 10);
    Thr.at(i)->SetDirectory(0);
    Noise.at(i)->SetDirectory(0);
    
    string fname = prefix + chipID[i] + suffix;
    const char* filename = fname.c_str();
    
    ifstream myfile (filename);
    if(myfile.good()){
      while (!myfile.eof()){
	myfile>>row>>column>>th>>sigma>>chi2;

	Thr.at(i)->Fill(th);
	Noise.at(i)->Fill(sigma);
      }
      
      myfile.close();
      
    }
    else{
      //scrivere qualcosa per evitare il crash nel caso in cui ci siano dei chip sbondati
      cout << "Unable to open file"<<endl;
      cout<<myfile.eof()<<endl;
      cout<<myfile.fail()<<endl;
      cout<<myfile.bad()<<endl;
    }

    c->cd(i + 1);
    Thr.at(i)->GetXaxis()->SetTitle("Threshold [e]");
    Thr.at(i)->GetYaxis()->SetTitle("Number of pixels");
    Thr.at(i)->Fit("gaus", "qrm+", "", 0, 300);
    fitThr.at(i) = (TF1*)Thr.at(i)->GetListOfFunctions()->FindObject("gaus");
    Thr.at(i)->Draw();
    Thr.at(i)->Write();
    TLegend* leg_thr = new TLegend(0.3, 0.7, 1.0, 0.95);
    leg_thr->AddEntry(Thr.at(i), Form("Threshold on chip 0x%s", chipID[i].c_str()));
    leg_thr->AddEntry(fitThr.at(i), Form("Fit on chip 0x%s", chipID[i].c_str()));
    leg_thr->AddEntry((TObject*)0, Form("p0 = %f", fitThr.at(i)->GetParameter(0)), " ");
    leg_thr->AddEntry((TObject*)0, Form("p1 = %f", fitThr.at(i)->GetParameter(1)), " ");
    leg_thr->AddEntry((TObject*)0, Form("p2 = %f", fitThr.at(i)->GetParameter(2)), " ");
    leg_thr->SetTextSize(.05);
    leg_thr->Draw("same");
    
    cnoise->cd(i + 1);
    Noise.at(i)->GetXaxis()->SetTitle("Threshold [e]");
    Noise.at(i)->GetYaxis()->SetTitle("Number of pixels");
    Noise.at(i)->Fit("gaus","q");
    fitNoise.at(i) = (TF1*)Noise.at(i)->GetListOfFunctions()->FindObject("gaus");
    Noise.at(i)->Draw();
    Noise.at(i)->Write();

    TLegend* leg_noise = new TLegend(0.3, 0.7, 1.0, 0.95);
    leg_noise->AddEntry(Noise.at(i), Form("Noise on chip 0x%s", chipID[i].c_str()));
    leg_noise->AddEntry(fitNoise.at(i), Form("Fit on chip 0x%s", chipID[i].c_str()));
    leg_noise->AddEntry((TObject*)0, Form("p0 = %f", fitNoise.at(i)->GetParameter(0)), " ");
    leg_noise->AddEntry((TObject*)0, Form("p1 = %f", fitNoise.at(i)->GetParameter(1)), " ");
    leg_noise->AddEntry((TObject*)0, Form("p2 = %f", fitNoise.at(i)->GetParameter(2)), " ");
    leg_noise->SetTextSize(.05);
    leg_noise->Draw("same");
  }

  c->Write();
  cnoise->Write();
  hfile->Write();
  hfile->Close();

  c->SaveAs("Threshold.png", "recreate");
  cnoise->SaveAs("Noise.png", "recreate");
  
  ofstream outfile;
  outfile.open("fit_data.txt");
  for(int k = 0; k < Nchip; k++){
    outfile<<"0x"<<chipID[k]<<" "<<fitThr.at(k)->GetParameter(1)<<" "<<fitThr.at(k)->GetParameter(2)<<" "<<fitNoise.at(k)->GetParameter(1)<<" "<<fitNoise.at(k)->GetParameter(2)<<"\n";
  }
  outfile.close();
}

void Threshold_temp(const char *dirFile){
  ifstream in_dir(dirFile);
  if(!in_dir){
    cout<<"Il file "<<dirFile<<" non esiste"<<endl;
  }

  TTree logtree("logbook", "Logbook file");
  int HIC, Vcasn, Ithr, Vcasn2;
  double Vbb, Temp;
  string Path;
  const char *pPath;
  string chip;
  double thr, erThr, noise, erNoise;
  double Thr[Nchip], ErThr[Nchip], Noise[Nchip], ErNoise[Nchip];
  double mean;
  double dev;
  
  logtree.Branch("V_casn", &Vcasn, "Vcasn/I");
  logtree.Branch("I_thr", &Ithr, "Ithr/I");
  logtree.Branch("V_back_bias", &Vbb, "Vbb/D");
  logtree.Branch("V_casn2", &Vcasn2, "Vcasn2/I");
  logtree.Branch("Temperature", &Temp, "Temp/D");
  logtree.Branch("HIC_number", &HIC, "HIC/D");
  logtree.Branch("Path_to_data", &Path, "Path/C");
  logtree.Branch("Threshold", &Thr[0], "Thr[10]/D");
  logtree.Branch("Error_Threshold", &ErThr[0], "ErThr[10]/D");
  logtree.Branch("Noise", &Noise[0], "Noise[10]/D");
  logtree.Branch("Error_Noise", &ErNoise[0], "ErNoise[10]/D");
  logtree.Branch("Mean_Threshold_on_HIC", &mean, "mean/D");
  logtree.Branch("Standard_Deviation_on_HIC", &dev, "dev/D");

  string prova;
  while(in_dir>>HIC>>Vcasn>>Ithr>>Vbb>>Vcasn2>>Temp>>Path){
    cout<<Path<<endl;
    pPath = Path.c_str();
    gSystem->cd(pPath);
    threshold(pPath);
    ifstream in("fit_data.txt");
    if(!in){
      cout<<"Il file "<<"fit_data.txt"<<" non esiste"<<endl;
    }
    int k = 0;
    mean = 0;
    dev = 0;
    while(in>>chip>>thr>>erThr>>noise>>erNoise){
      Thr[k] = thr;
      ErThr[k] = erThr;
      Noise[k] = noise;
      ErNoise[k] = erNoise;
      mean += thr;
      k++;
    }
    mean = mean / Nchip;
    for(int i = 0; i < Nchip; i++){
      dev += pow(Thr[i] - mean, 2.);
    }
    dev = sqrt(dev / Nchip);
    cout<<"Mean Threshold on HIC "<<mean<<endl;
    cout<<"Standard deviation on HIC "<<dev<<endl;
    in.close();
    /*    gSystem->cd("../");
    ifstream in_reg("register_dump.txt");
    if(!in_reg){
      cout<<"Il file "<<"register_dump.txt"<<" non esiste"<<endl;
    }
    int i = 0;
    string reg, val; 
    while(in_reg>>chip>>reg>>val){
      for(int i = 0; i < Nchip; i++){
	if(chip == chipID[i]){
	  if(reg == "0x604"){
	    Vcasn[i] = val;
	  }
	  if(reg == "0x607"){
	    Vcasn_2[i] = val;
	  }
	  if(reg == "0x60e"){
	    Ithr[i] = val;
	  }
	}
      }
      }*/
    gSystem->cd("../../../../");
    logtree.Fill(); 
  }
    
  in_dir.close();
  
  TFile* outfile = new TFile("Threshold_Temp_Parameters.root", "recreate");
  logtree.Write();
  outfile->Write();
  outfile->Close();
}
