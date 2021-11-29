#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TPad.h"


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

void threshold_tuning(const char* path_to_file="Threshold_Parameters.root",const double threshold=150.,const double VBB_ref=0 /*possible values: 0, -1, -3*/){
	
	//testing the selection of correct VBB value
	bool vbb_test=false;
	if(VBB_ref==0 || VBB_ref==-1 || VBB_ref==-3) vbb_test=true;

	if(!vbb_test){ 
		std::cout<<"Warning! The selected VBB value is not correct. The allowed values are 0, -1 and -3."<<std::endl;

		return;
	} //if the value of vbb is not correct, the software stops



	gROOT->SetBatch(kTRUE);

	const int NCHIP=10;

	//chipname list
	int chipname[NCHIP]={0x70,0x71,0x72,0x73,0x74,0x78,0x79,0x7a,0x7b,0x7c};

	//opening file and catching branches
	TFile* infile=new TFile(path_to_file);
	TTree* tree=(TTree*)infile->Get("logbook");

	int vcasn[NCHIP],ithr[NCHIP],vcasn2[NCHIP],HICnum; //values that should be converted to int
	double vbb;
	char* path;
	double thr[NCHIP],thr_err[NCHIP],noise[NCHIP],noise_err[NCHIP],meanHIC,stddevHIC;

	tree->SetBranchAddress("V_casn",&vcasn);
	tree->SetBranchAddress("V_casn_2",&vcasn2);
	tree->SetBranchAddress("V_back_bias",&vbb);
	tree->SetBranchAddress("I_thr",&ithr);
	tree->SetBranchAddress("HIC_number",&HICnum);
	tree->SetBranchAddress("Threshold",thr);
	tree->SetBranchAddress("Error_Threshold",thr_err);
	tree->SetBranchAddress("Noise",noise);
	tree->SetBranchAddress("Error_Noise",noise_err);

	//vectors to collect data from tree
	std::vector<int> vcasn_val[NCHIP];
	std::vector<int> vcasn2_val[NCHIP];
	std::vector<double> vbb_val;
	std::vector<int> ithr_val[NCHIP];
	std::vector<double> thr_val[NCHIP];
	std::vector<double> thr_err_val[NCHIP];
	std::vector<double> noise_val[NCHIP];
	std::vector<double> noise_err_val[NCHIP];
	//getting information from entry 0
	tree->GetEntry(0);

	std::cout<<"Running analysis on threshold on HIC "<<HICnum<<std::endl;


	//loop on the tree to get data and fill the vectors

	for(int iEn=0;iEn<tree->GetEntries();++iEn){
		tree->GetEntry(iEn);
		//std::cout<<thr[0]<<"\t"<<thr[9]<<std::endl;
		
		vbb_val.push_back(vbb);
		
		for(int iChip=0;iChip<NCHIP;++iChip){
			vcasn_val[iChip].push_back(vcasn[iChip]);
			vcasn2_val[iChip].push_back(vcasn2[iChip]);
			ithr_val[iChip].push_back(ithr[iChip]);
			thr_val[iChip].push_back(thr[iChip]);
			thr_err_val[iChip].push_back(thr_err[iChip]);
			noise_val[iChip].push_back(noise[iChip]);
			noise_err_val[iChip].push_back(noise_err[iChip]);
		}	

	}
	
	
	
	//open output configuration file

	std::stringstream outfile_name;
	outfile_name<<"config_HIC"<<HICnum<<"_vbb"<<VBB_ref<<"_threshold"<<threshold<<".conf";

	std::ofstream conf_output(outfile_name.str().c_str(), std::ofstream::out);

	conf_output<<"; Threshold configuration for HIC "<<HICnum<<std::endl
		<<"; Threshold value expected "<<threshold<<std::endl
		<<"; Values for VBB = "<<VBB_ref<<std::endl<<std::endl;	

	//inserting VCLIP values for different VBB values
	conf_output<<"[0xf0f]"<<std::endl;
	if(VBB_ref==0)
		conf_output<<"VCLIP = 0x0"<<std::endl<<std::endl;
	if (VBB_ref==-1)
		conf_output<<"VCLIP = 0x23"<<std::endl<<std::endl;
	if (VBB_ref==-3)
		conf_output<<"VCLIP = 0x3c"<<std::endl<<std::endl;


	//searching VCASN best value for chips: run on fixed value of ITHR (50 is the default for ALICE)

	const double ithr_ref=50.;

	double vcasn_bestvalue[NCHIP];

	for(int iChip=0;iChip<NCHIP;++iChip){
		vcasn_bestvalue[iChip]=-9999; //negative flag to identify missing chips
		}

		for(int iChip=0;iChip<NCHIP;++iChip){
		//find vcasn values available for the given Vbb and chip and create a list
		std::vector<double>vcasn_list;
		
		for(int iV=0;iV<vcasn_val[iChip].size();++iV){
			bool add_vcasn=true;
			//std::cout<<vcasn_val.at(iV)<<"\t";
			if(vbb_val.at(iV)==VBB_ref){
				for(int iVC=0;iVC<vcasn_list.size();++iVC){
					if(vcasn_val[iChip].at(iV)==vcasn_list.at(iVC)) add_vcasn=false;
				}
			} else
			add_vcasn=false;

			if(add_vcasn)vcasn_list.push_back(vcasn_val[iChip].at(iV));
		}

		std::cout<<"Thresholds for chip 0x"<<std::hex<<chipname[iChip]<<std::endl;
		//loop on each vcasn. Using ITHR=50 as reference. Comparing with threshold to find the closest value
		int closest_index_vcasn=-9999;
		double diff=999999;
		for(int iV=0;iV<vcasn_list.size();++iV){
			for(int iA=0;iA<vcasn_val[iChip].size();++iA){
				if(vcasn_val[iChip].at(iA)==vcasn_list.at(iV) &&ithr_val[iChip].at(iA)==ithr_ref && vbb_val.at(iA)==VBB_ref){
					if(abs(thr_val[iChip].at(iA)-threshold)<diff) {
						closest_index_vcasn=iA;
						diff=abs(thr_val[iChip].at(iA)-threshold);
					}
				}
			}
		}

		double vcasn_ref=vcasn_val[iChip].at(closest_index_vcasn);
		std::cout<<" The VCASN value that gives the threshold value closest to "<<threshold<<" is 0x"<<std::hex<<(int)vcasn_ref<<std::endl;
		vcasn_bestvalue[iChip]=vcasn_ref;
	}

	//prepare graph for fit
	TGraphErrors* graph_val[NCHIP];

	//loop on chips to find the best ITHR (VCASN fixed)
	for(int iChip=0;iChip<NCHIP;++iChip){
		if(vcasn_bestvalue[iChip]>0){
			conf_output<<"[0x"<<std::hex<<chipname[iChip]<<"]"<<std::endl<<
				"VCASN = 0x"<<std::hex<<(int)vcasn_bestvalue[iChip]<<std::endl<<
				"VCASN2 = 0x"<<std::hex<<(int)(vcasn_bestvalue[iChip]+12)<<std::endl;

			graph_val[iChip]=new TGraphErrors(12); //number of points will be fixed?
			int entry=0;
			for(int iA=0;iA<vcasn_val[iChip].size();++iA){
	
				if(vcasn_val[iChip].at(iA)==vcasn_bestvalue[iChip] && vbb_val.at(iA)==VBB_ref){
					graph_val[iChip]->SetPoint(entry,ithr_val[iChip].at(iA),thr_val[iChip].at(iA));
					graph_val[iChip]->SetPointError(entry,0,thr_err_val[iChip].at(iA));
					++entry;
				}
			}
			//labeling plots
			graph_val[iChip]->SetTitle(Form("0x%x VCASN %0.0f",chipname[iChip],vcasn_bestvalue[iChip]));
			graph_val[iChip]->GetXaxis()->SetTitle("ITHR");
			graph_val[iChip]->GetYaxis()->SetTitle("threshold [e-]");
	
			//graph fit with linear function
			graph_val[iChip]->Fit("pol1","S");
	
			double chisquare=graph_val[iChip]->GetFunction("pol1")->GetChisquare();

			//TO DO: improve check on fit convergence.
			if(chisquare){
				//if fit converges, the parameters are used to finf the best guess for ITHR
				double p0=graph_val[iChip]->GetFunction("pol1")->GetParameter(0);
				double p1=graph_val[iChip]->GetFunction("pol1")->GetParameter(1);
		
				double p0_err=graph_val[iChip]->GetFunction("pol1")->GetParError(0);
				double p1_err=graph_val[iChip]->GetFunction("pol1")->GetParError(1);
		
				//std::cout<<"p0 "<<p0<<" p0_err "<<p0_err<<" p1 "<<p1<<" p1_err "<<p1_err<<std::endl;
		
				//best value calculation
				double ithr_threshold=(threshold-p0)/p1;
		
				std::cout<<"ITHR best value to have a threshold of "<<threshold <<" e- should be 0x"<<std::hex<<(int)round(ithr_threshold)<<std::endl;
				//writing conf file
				conf_output<<"ITHR = 0x"<<std::hex<<(int)round(ithr_threshold)<<std::endl<<std::endl;
				std::cout<<"Expected threshold for chip 0x"<<std::hex<<chipname[iChip]<<" "<<std::dec<<p0+p1*ithr_threshold<<" (FIT)"<<std::endl;
			}
			else{
				//if fit does not converge, the existing data are used to find the best value of ITHR. A warning is inserted as a comment on conf file
				std::cout<<"Fit does not converge. Using available data to find the best ITHR value"<<std::endl;
				int closest_index_ithr=-9999;
				double diff=999999;
				for(int iA=0;iA<ithr_val[iChip].size();++iA){
					if(vcasn_val[iChip].at(iA)==vcasn_bestvalue[iChip] && vbb_val.at(iA)==VBB_ref){
						if(abs(thr_val[iChip].at(iA)-threshold)<diff) {
							closest_index_ithr=iA;
							diff=abs(thr_val[iChip].at(iA)-threshold);
						}
					}
				}
			std::cout<<"ITHR best value to have a threshold of "<<threshold <<" e- should be 0x"<<std::hex<<(int)ithr_val[iChip].at(closest_index_ithr)<<std::endl;
			conf_output<<"ITHR = 0x"<<std::hex<<(int)ithr_val[iChip].at(closest_index_ithr)<<std::endl<<
			"; Warning! Fit doesn't converge. Using available data to find the best value"<<std::endl<<std::endl;
			std::cout<<"Expected threshold for chip 0x"<<std::hex<<chipname[iChip]<<" "<<std::dec<<(int)thr_val[iChip].at(closest_index_ithr)<<" err "<<(int)thr_err_val[iChip].at(closest_index_ithr)<<" (MEASURED)"<<std::endl;
			
			}
		}
	}
	
	//printout of graphs to check fit results
	
	TCanvas* can=new TCanvas();
	can->Divide(5,2);

	for(int iChip=0;iChip<NCHIP;++iChip){
		can->cd(iChip+1);
		graph_val[iChip]->SetMarkerStyle(22);
		graph_val[iChip]->Draw("AP");

	}
	can->Print(Form("fit_HIC%d_vbb%0.0f_threshold%0.0f.pdf",HICnum,VBB_ref,threshold));
	
	//test: fit on VCASN distribution
	//fixed ithr=50
	TGraphErrors* vcasn_graph[NCHIP];
	for(int iChip=0;iChip<NCHIP;++iChip){
		vcasn_graph[iChip]=new TGraphErrors(10);
		int entry=0;
		for(int iA=0;iA<vcasn_val[iChip].size();++iA){
			if(ithr_val[iChip].at(iA)==ithr_ref && vbb_val.at(iA)==VBB_ref && thr_val[iChip].at(iA)>0){
				vcasn_graph[iChip]->SetPoint(entry,vcasn_val[iChip].at(iA),thr_val[iChip].at(iA));
				vcasn_graph[iChip]->SetPointError(entry,0,thr_err_val[iChip].at(iA));
				++entry;
			}
		}
		vcasn_graph[iChip]->Fit("expo");

		double p0_exp=vcasn_graph[iChip]->GetFunction("expo")->GetParameter(0);
		double p1_exp=vcasn_graph[iChip]->GetFunction("expo")->GetParameter(1);

		int vcasn_guess=round((log(threshold)-p0_exp)/p1_exp);

		std::cout<<"Fit guess of VCASN for chip 0x"<<std::hex<<chipname[iChip]<<" is "<<std::dec<<vcasn_guess<<std::endl;

		can->cd(iChip+1);
		vcasn_graph[iChip]->SetTitle(Form("0x%x;VCASN;threshold [e-]",chipname[iChip]));
		vcasn_graph[iChip]->SetMarkerStyle(22);
		vcasn_graph[iChip]->Draw("AP");
	}

	can->Print("test_vcasn_fit.pdf");


	//test: optimisation of procedure: interpolation on VCASN and ITHR and surface design to find the best value
/*
	TH2I* graph_surf[NCHIP];
	TGraph2D* tgraph_surf[NCHIP];
	for(int iChip=0;iChip<NCHIP;++iChip){
		graph_surf[iChip]=new TH2I(Form("graph_surf_0x%x",chipname[iChip]),Form("0x%x;VCASN;ITHR;threshold [e-]",chipname[iChip]),10,50,80,10,50,80);
		tgraph_surf[iChip]=new TGraph2D(12);
		int entry=0;
		for(int iA=0;iA<vcasn_val.size();++iA){
			int xbin=graph_surf[iChip]->GetXaxis()->FindBin(vcasn_val.at(iA));
			int ybin=graph_surf[iChip]->GetYaxis()->FindBin(ithr_val.at(iA));
			
			if(vbb_val.at(iA)==VBB_ref && graph_surf[iChip]->GetBinContent(xbin,ybin)==0){
				graph_surf[iChip]->Fill(vcasn_val.at(iA),ithr_val.at(iA),thr_val[iChip].at(iA));
				//std::cout<<vcasn_val.at(iA)<<"\t"<<ithr_val.at(iA)<<"\t"<<thr_val[iChip].at(iA)<<std::endl;
				
				tgraph_surf[iChip]->SetPoint(entry,vcasn_val.at(iA),ithr_val.at(iA),thr_val[iChip].at(iA));
				++entry;
				
	
			}
		}
		//interpolate
		//for(int iVCASN=50;iVCASN<80;++iVCASN){
			//graph_surf[iChip]->Interpolate(56,75);
		//}
	
	}

	can->Clear();
	can->Divide(5,2);
	can->Print("test_surf.pdf[");
	for(int iChip=0;iChip<NCHIP;++iChip){
		TPad* pad=(TPad*)can->cd(iChip+1);
		pad->SetPhi(-45);
		graph_surf[iChip]->SetStats(0);
		graph_surf[iChip]->Draw("surf1");

		
	}
	can->Print("test_surf.pdf");
	
	for(int iChip=0;iChip<NCHIP;++iChip){
		TPad* pad=(TPad*)can->cd(iChip+1);
		pad->SetPhi(-45);
		
		tgraph_surf[iChip]->SetTitle(Form("0x%x;VCASN;ITHR",chipname[iChip]));
		tgraph_surf[iChip]->Draw("PCOL");
	}
	can->Print("test_surf.pdf");
	
	can->Print("test_surf.pdf]");
*/
	return;
}