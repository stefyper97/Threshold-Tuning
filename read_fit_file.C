#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void readfile(const char* filelist){
	//info in filelist: stave #, Vbb, PATH (VCASN and ITHR are stored on the register dump file)


	const int NCHIP=10;
	const int chipname[NCHIP]={0x70,0x71,0x72,0x73,0x74,0x78,0x79,0x7a,0x7b,0x7c};

	//register mapping
	const int ithr_reg=0x60e;
	const int vcasn_reg=0x604;
	const int vcasn2_reg=0x607;

	//opening tree

	int HIC, Vcasn[NCHIP], Ithr[NCHIP], Vcasn2[NCHIP];
 	double Vbb;
  	//string Path;
  	//const char *pPath;
  	//string chip;
  	double thr, erThr, noise, erNoise;
  	double Thr[NCHIP], ErThr[NCHIP], Noise[NCHIP], ErNoise[NCHIP];
  	double mean;
  	double dev;

	TTree logtree("logbook", "Logbook file");
	logtree.Branch("V_casn", &Vcasn[0], "Vcasn[10]/I");
  	logtree.Branch("I_thr", &Ithr[0], "Ithr[10]/I");
  	logtree.Branch("V_back_bias", &Vbb, "Vbb/D");
  	logtree.Branch("V_casn_2", &Vcasn2[0], "Vcasn2[10]/I");
  	logtree.Branch("HIC_number", &HIC, "HIC/I");
  	//logtree.Branch("Path_to_data", &Path, "Path/C");
  	logtree.Branch("Threshold", &Thr[0], "Thr[10]/D");
  	logtree.Branch("Error_Threshold", &ErThr[0], "ErThr[10]/D");
  	logtree.Branch("Noise", &Noise[0], "Noise[10]/D");
  	logtree.Branch("Error_Noise", &ErNoise[0], "ErNoise[10]/D");
  	logtree.Branch("Mean_Threshold_on_HIC", &mean, "mean/D");
  	logtree.Branch("Standard_Deviation_on_HIC", &dev, "dev/D");


	//opening filelist 
	std::ifstream input_file(filelist,std::ifstream::in);
	//int stave_id, Vbb;
	std::string path;
	while(!input_file.eof()){
		for(int iChip=0;iChip<NCHIP;++iChip){
			Vcasn[iChip]=-9999;
			Ithr[iChip]=-9999;
			Vcasn2[iChip]=-9999;
			Thr[iChip]=-9999;
			ErThr[iChip]=-9999;
			Noise[iChip]=-9999;
			ErNoise[iChip]=-9999;


		}
		input_file>>HIC>>Vbb>>path;

		//read registers from file. Opening file
		std::stringstream regfile_name;
		regfile_name<<path<<"/register_dump.txt";

		std::ifstream reg_file(regfile_name.str().c_str(),std::ifstream::in);
		int chip_id,reg_num,reg_val;
		//loop on register values
		while(!reg_file.eof()){
			reg_file>>std::hex>>chip_id>>std::hex>>reg_num>>std::hex>>reg_val;
			for(int iChip=0;iChip<NCHIP;++iChip){
				if ((chip_id&0x7f)==chipname[iChip] && reg_num==ithr_reg && reg_val!=0xffff){
					Ithr[iChip]=reg_val;
				}
				else if((chip_id&0x7f)==chipname[iChip] && reg_num==vcasn_reg && reg_val!=0xffff){
					Vcasn[iChip]=reg_val;
				}
				else if((chip_id&0x7f)==chipname[iChip] && reg_num==vcasn2_reg && reg_val!=0xffff){
					Vcasn2[iChip]=reg_val;
				}
			}
		}
		
		//read threshold values
		std::stringstream thrfile_name;
		thrfile_name<<path<<"/fit.txt";
		std::ifstream thr_file(thrfile_name.str().c_str(),std::ifstream::in);
		std::string buffer;
		getline(thr_file,buffer);

		int avg,stddev;
		while(!thr_file.eof()){
			thr_file>>std::hex>>chip_id>>std::dec>>avg>>std::dec>>stddev;
			//std::cout<<"0x"<<std::hex<<chip_id<<"\t"<<std::dec<<avg<<"\t"<<stddev<<std::endl;
			for(int iChip=0;iChip<NCHIP;++iChip){
				if((chip_id&0x7f)==chipname[iChip]){
					Thr[iChip]=avg;
					ErThr[iChip]=stddev;	
				}
			}
		}
		

		logtree.Fill(); 

	}
	//opening root file and saving data from tree
	TFile* outfile = new TFile(Form("Threshold_Parameters_HIC_0%d_Vbb_%0.0fV.root", HIC, Vbb), "recreate");
	logtree.Write();
  	outfile->Write();
  	outfile->Close();
	

	return;
}
