#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TROOT.h"
#include "TFitResultPtr.h"


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

void threshold_tuning(const char* path_to_file="Threshold_Parameters.root",const double threshold=150.,const double VBB_ref=0 /*possible values: 0, -3*/){
	gROOT->SetBatch(kTRUE);

	const int NCHIP=10;

	//chipname list
	int chipname[NCHIP]={0x70,0x71,0x72,0x73,0x74,0x78,0x79,0x7a,0x7b,0x7c};

	//opening file and catching branches
	TFile* infile=new TFile(path_to_file);
	TTree* tree=(TTree*)infile->Get("logbook");

	double vcasn,ithr,vbb,vcasn2,HICnum; //values that should be converted to int
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
	std::vector<double> vcasn_val;
	std::vector<double> vcasn2_val;
	std::vector<double> vbb_val;
	std::vector<double> ithr_val;
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
		vcasn_val.push_back(vcasn);
		vcasn2_val.push_back(vcasn2);
		vbb_val.push_back(vbb);
		ithr_val.push_back(ithr);
		for(int iChip=0;iChip<NCHIP;++iChip){
			thr_val[iChip].push_back(thr[iChip]);
			thr_err_val[iChip].push_back(thr_err[iChip]);
			noise_val[iChip].push_back(noise[iChip]);
			noise_err_val[iChip].push_back(noise_err[iChip]);
		}	

	}
	
	//find vcasn values available for the given Vbb and create a list
	std::vector<double>vcasn_list;
	
	for(int iV=0;iV<vcasn_val.size();++iV){
		bool add_vcasn=true;
		//std::cout<<vcasn_val.at(iV)<<"\t";
		if(vbb_val.at(iV)==VBB_ref){
			for(int iVC=0;iVC<vcasn_list.size();++iVC){
				if(vcasn_val.at(iV)==vcasn_list.at(iVC)) add_vcasn=false;
			}
		} else
		add_vcasn=false;

		if(add_vcasn)vcasn_list.push_back(vcasn_val.at(iV));
	
	}
	
	//open output configuration file

	std::stringstream outfile_name;
	outfile_name<<"config_HIC"<<HICnum<<"_vbb"<<VBB_ref<<"_threshold"<<threshold<<".conf";

	std::ofstream conf_output(outfile_name.str().c_str(), std::ofstream::out);

	conf_output<<"; Threshold configuration for HIC "<<HICnum<<std::endl
		<<"; Threshold value expected "<<threshold<<std::endl
		<<"; Values for VBB = "<<VBB_ref<<std::endl<<std::endl;	

	//searching VCASN best value for chips: run on fixed value of ITHR (50 is the default for ALICE)

	const double ithr_ref=50.;

	double vcasn_bestvalue[NCHIP];

	for(int iChip=0;iChip<NCHIP;++iChip){
		vcasn_bestvalue[iChip]=-9999; //negative flag to identify missing chips
	}

	for(int iChip=0;iChip<NCHIP;++iChip){
		std::cout<<"Thresholds for chip 0x"<<std::hex<<chipname[iChip]<<std::endl;
		//loop on each vcasn. Using ITHR=50 as reference. Comparing with threshold to find the closest value
		int closest_index_vcasn=-9999;
		double diff=999999;
		for(int iV=0;iV<vcasn_list.size();++iV){
			for(int iA=0;iA<vcasn_val.size();++iA){
				if(vcasn_val.at(iA)==vcasn_list.at(iV) &&ithr_val.at(iA)==ithr_ref && vbb_val.at(iA)==VBB_ref){
					if(abs(thr_val[iChip].at(iA)-threshold)<diff) {
						closest_index_vcasn=iA;
						diff=abs(thr_val[iChip].at(iA)-threshold);
					}
				}
			}
		}

		double vcasn_ref=vcasn_val.at(closest_index_vcasn);
		std::cout<<" The VCASN value that gives the threshold value closest to "<<threshold<<" is "<<vcasn_ref<<std::endl;
		vcasn_bestvalue[iChip]=vcasn_ref;
	}

	//prepare graph for fit
	TGraphErrors* graph_val[NCHIP];

	//loop on chips to find the best ITHR (VCASN fixed)
	for(int iChip=0;iChip<NCHIP;++iChip){
		if(vcasn_bestvalue[iChip]>0){
			conf_output<<"[0x"<<std::hex<<chipname[iChip]<<"]"<<std::endl<<
				"VCASN = 0x"<<std::hex<<vcasn_bestvalue[iChip]<<std::endl<<
				"VCASN2 = 0x"<<std::hex<<vcasn_bestvalue[iChip]+12<<std::endl;

			graph_val[iChip]=new TGraphErrors(12); //number of points will be fixed?
			int entry=0;
			for(int iA=0;iA<vcasn_val.size();++iA){
	
				if(vcasn_val.at(iA)==vcasn_bestvalue[iChip] && vbb_val.at(iA)==VBB_ref){
					graph_val[iChip]->SetPoint(entry,ithr_val.at(iA),thr_val[iChip].at(iA));
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
		
				std::cout<<"ITHR best value to have a threshold of "<<threshold <<" e- should be "<<round(ithr_threshold)<<std::endl;
				//writing conf file
				conf_output<<"ITHR = 0x"<<std::hex<<round(ithr_threshold)<<std::endl<<std::endl;
			}
			else{
				//if fit does not converge, the existing data are used to find the best value of ITHR. A warning is inserted as a comment on conf file
				std::cout<<"Fit does not converge. Using available data to find the best ITHR value"<<std::endl;
				int closest_index_ithr=-9999;
				double diff=999999;
				for(int iA=0;iA<ithr_val.size();++iA){
					if(vcasn_val.at(iA)==vcasn_bestvalue[iChip] && vbb_val.at(iA)==VBB_ref){
						if(abs(thr_val[iChip].at(iA)-threshold)<diff) {
							closest_index_ithr=iA;
							diff=abs(thr_val[iChip].at(iA)-threshold);
						}
					}
				}
			std::cout<<"ITHR best value to have a threshold of "<<threshold <<" e- should be "<<ithr_val.at(closest_index_ithr)<<std::endl;
			conf_output<<"ITHR = 0x"<<std::hex<<ithr_val.at(closest_index_ithr)<<std::endl<<
			"; Warning! Fit doesn't converge. Using available data to find the best value"<<std::endl<<std::endl;
			
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
	can->Print(Form("fit_HIC%0.0f_vbb%0.0f_threshold%0.0f.pdf",HICnum,VBB_ref,threshold));
	
	
	return;
}