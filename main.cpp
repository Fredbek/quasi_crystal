#include <iostream>
#include <math.h>
#include <ostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <algorithm>

double h_temp = 0.25;
int qtd = 20;

std::string OutFilePath;
std::string lattices_path_init = "lattices/";
std::string data_arch_name = "appxt_", extension = "_connect.dat";
std::string Data_out_Path = "Results/Data/";
std::string sites[4] = {"239", "1393", "8119", "47321"};

//int TermSteps = 10; //Quantidade de passos para a termalização
const static int MedPerBlock = 1500, blocks = 20, Med = blocks*MedPerBlock;

double J = -1; //Determina a intensidade e tipo das interaçes interspin

//////////////// Definições referentes à geradores de números aleatório
std::random_device rd;
std::mt19937 eng(rd());
/////////////////////////////////////////

void print(auto texto){
	std::cout << texto << std::endl;
}

double Media(double* Array, int array_lenght){
	double soma = 0;
	
	for(int i = 0; i < array_lenght; i++){
		soma += Array[i]/array_lenght;
	}
	
	return soma;
}

double DesvPadMed(double* Array, int array_lenght){
	double soma = 0, media_dados = Media(Array, array_lenght);

	for(int i = 0; i < array_lenght; i++){
		soma += pow((Array[i] - media_dados), 2);
	}
	
	return sqrt(soma/(array_lenght));
}

class Lattice{

private:
	int** big_spin = NULL;
	bool* bf_rede = NULL;
	int* phase_lat = NULL;
	int** Viz = NULL;

	int* Termalization_Steps = NULL;

	double* Temperatures = NULL;
	int qtd_redes;
	double* Probabilities = NULL;

	std::string data_path = lattices_path_init + str_sites_qtd + "/" + data_arch_name + str_sites_qtd + extension;

	std::string str_sites_qtd;

	int int_sites_qtd = std::stoi(str_sites_qtd);

	//double data_crua[Med][3];// Pega os dados das medidas na ordem (en, en^2, m)
	double en_proc[blocks], mg_proc[blocks], cv_proc[blocks], sus_proc[blocks];
	
	//double raw_data[Med][2]; //Data crua com entradas (m^2, m^4)
	double mean_sqrd_mag[blocks]; //blocos com medidas das grandezas desejadas
	double mean_qdrt_mag[blocks];
	double binder_cumulant[blocks];
	
	std::vector<int> bf_cluster; //Buffer do cluster que será girado

	int WolffTermSteps = floor(int_sites_qtd/20); //Quantidade de passos para a termalização
	int MetTermSteps = 50*int_sites_qtd;
	int TermSteps;

	/// Definições do rng
	std::uniform_real_distribution<double> dist{ 0, 1.0 };
	std::uniform_int_distribution<int> dice{0, (int_sites_qtd - 1)};

	int got = 0;
	int tries = 0;

public:

	Lattice(double fst_int_in, double fst_int_fim, double fst_int_step, double scd_int_fim, double scd_int_step, double trd_int_fim, double trd_int_step ,std::string Comp)
		:str_sites_qtd(Comp)
	{	
		qtd_redes = 0;
		for(double i = fst_int_in; i < fst_int_fim; i += fst_int_step){
			qtd_redes += 1;
		}
		for(double i = fst_int_fim; i < scd_int_fim; i += scd_int_step){
			qtd_redes += 1;
		}
		for(double i = scd_int_fim; i < trd_int_fim; i += trd_int_step){
			qtd_redes += 1;
		}

		print(qtd_redes);

		Temperatures = new double[qtd_redes];
		Termalization_Steps = new int[qtd_redes];
		big_spin = new int*[qtd_redes];
		bf_rede = new bool[int_sites_qtd];

		phase_lat = new int[int_sites_qtd];

		Viz = new int*[int_sites_qtd];

		for(int i = 0; i < qtd_redes; i++){
			big_spin[i] = new int[int_sites_qtd];
		}

		for(int i = 0; i < int_sites_qtd; i++){
			Viz[i] = new int[8];
			bf_rede[i] = false;
			phase_lat[i] = 0;
		}

		int indx = 0;
		for(double i = fst_int_in; i < fst_int_fim; i += fst_int_step){
			Temperatures[indx] = i;
			Probabilities[indx] = 1 - exp(-2*abs(J)/Temperatures[indx]);

			Termalization_Steps[indx] = (Temperatures[indx] > 2 and Temperatures[indx] < 2.44) ? (MetTermSteps) : (MetTermSteps);
			indx += 1;
		}
		for(double i = fst_int_fim; i < scd_int_fim; i += scd_int_step){
			Temperatures[indx] = i;
			Probabilities[indx] = 1 - exp(-2*abs(J)/Temperatures[indx]);

			Termalization_Steps[indx] = (Temperatures[indx] > 2 and Temperatures[indx] < 2.44) ? (MetTermSteps) : (MetTermSteps);
			indx += 1;
		}
		for(double i = scd_int_fim; i <= trd_int_fim; i += trd_int_step){
			Temperatures[indx] = i;
			Probabilities[indx] = 1 - exp(-2*abs(J)/Temperatures[indx]);

			Termalization_Steps[indx] = (Temperatures[indx] > 2 and Temperatures[indx] < 2.44) ? (MetTermSteps) : (MetTermSteps);
			indx += 1;
		}

		this->InitLattices();
		this->getViz();
		this->getPhase(0, 1);
		this->Parallel_Termalization();
	}

	~Lattice(){
		delete[] Temperatures;
		delete[] Termalization_Steps;
		delete[] Probabilities;
		delete[] phase_lat;
		delete[] bf_rede;

		for(int i = 0; i < int_sites_qtd; i++){
			delete[] Viz[i];
		}

		for(int i = 0; i < qtd_redes; i++){
			delete[] big_spin[i];
		}

		delete[] Viz;
		delete[] big_spin;

	}

	double i_rand() { ///Função que gera um inteiro aleatório entre 0 e int_sites_qtd
		return dice(eng);
	}
	
	double r_rand() { ///Função que gera uma real alatório entre 0 e 1
		return dist(eng);
	}

	void reset_bf_rede(){
		for(int i = 0; i < int_sites_qtd; i++){
			bf_rede[i] = false;
		}
	}

	void InitLattices(){
		for(int i = 0; i < qtd_redes; i++){
			for(int j = 0; j < int_sites_qtd; j++){
				big_spin[i][j] = 1;
			}
		}
	}

	void getViz(){ //Lê a tabela de vizinhos
	
		std::string VizPath = lattices_path_init + str_sites_qtd + "/appxt_" + str_sites_qtd + "_connect.dat";

		std::ifstream Data;

		Data.open(VizPath);

		for(int i = 0; i < int_sites_qtd; i++){
			for(int j = 0; j < 8; j++){
				Data >> Viz[i][j];
				Viz[i][j] = Viz[i][j] - 1;
			}
		}

		Data.close();
	}

	void getPhase(int site, int phase){
		phase_lat[site] = phase;
		for(int i = 0; i < 8; i++){
			if(Viz[site][i] >= 0 and phase_lat[Viz[site][i]] == 0){
				this->getPhase(Viz[site][i], -phase);
			}
		}
	}

	double Ttl_magnetization(int* Spin_Config){ //Devolve a  magnetização total para uma configuração de spin
		
		double mag = 0;
		
		for(int i = 0; i < int_sites_qtd; i++){
			mag += Spin_Config[i]*phase_lat[i];
		}

		return abs(mag);
	}
	
	double Site_energy(int* Spin_Config, int site){

		double energy = 0;

		for(int i = 0; i < 8; i++){
			if(Viz[site][i] >= 0){
				energy += -J*Spin_Config[site]*Spin_Config[Viz[site][i]];
			}
		}

		return energy;
	}

	double Lattice_Energy(int* Spin_Config){
		
		double energy = 0;

		for(int i = 0; i < int_sites_qtd; i++){
			energy += this->Site_energy(Spin_Config, i);
		}

		return energy/2;
	}

	void MetropolisStep(int lat_indx){// Faz uma varredura de Metropolis na rede de índice lat_indx
		int* Spin_Config = big_spin[lat_indx];

		for(int i = 0; i < int_sites_qtd; i++){
			double DE = -2*Site_energy(Spin_Config, i);

			if(DE <= 0){
				Spin_Config[i] = -Spin_Config[i];
			}
			else{
				double z = exp(-DE/Temperatures[lat_indx]);
				if(z >= r_rand()){
					Spin_Config[i] = -Spin_Config[i];
				}
			}
		}

	}

	void Parallel_Termalization(){

		int* Passos = new int[qtd_redes];

		for(int i = 0; i < qtd_redes; i++){
			Passos[i] = 0;
		}

		for(int cont = 0; cont < MetTermSteps/5; cont++){
			for(int cont = 0; cont < 5; cont++){
				for(int indx = 0; indx < qtd_redes; indx++){
					
					if(Passos[indx] > Termalization_Steps[indx]){
						this->MetropolisStep(indx);
						Passos[indx] += 1;
					}
				}

				for(int i = qtd_redes - 1; i >= 1; i -= 1){
					int rede_ant = i - 1;
					this->Parallel_Tempering(rede_ant, i);
				}
			
			}
		}


		delete[] Passos;
		std::cout << "Parallel Tempering Acceptance Ratio: " << 100.0*(got/tries) << "%" << std::endl;
	}

	int get_qtd_redes(){
		return qtd_redes;
	}

	void get_mean_values(int i, double* values){
		this->getProcData(i);
		values[0] = this->energy();
		values[1] = this->senergy();
		values[2] = this->magnet();
		values[3] = this->smagnet();
		values[4] = this->cv();
		values[5] = this->scv();
		values[6] = this->susc();
		values[7] = this->ssusc();
		values[8] = this->BinderCumulant();
		values[9] = this->sBinderCumulant();
	}

	double get_Temp(int ind){
		return Temperatures[ind];
	}

	void getProcData(int indx){

		int* Spin_Config = big_spin[indx];

		double Temp = Temperatures[indx];

		double med_en = 0, med_en_sq = 0, med_mg = 0, med_mg_sq = 0, med_mg_qd = 0, binder_ratio = 0;

		if(Temperatures[indx] < 2.0 or Temperatures[indx] > 2.44){//Seta a quantidade de passos de termalização
			TermSteps = MetTermSteps;
		}
		else{
			TermSteps = MetTermSteps;
		}

		for(int vez = 0; vez < blocks; vez++){

			med_en = 0;
			med_en_sq = 0;
			med_mg = 0;
			med_mg_sq = 0;
			med_mg_qd = 0;

			en_proc[vez] = 0;
			mg_proc[vez] = 0;
			cv_proc[vez] = 0;
			sus_proc[vez] = 0;
			
			for(int medida = 0; medida < MedPerBlock; medida++){
				en_proc[vez] += this->Lattice_Energy(Spin_Config)/MedPerBlock;
				mg_proc[vez] += this->Ttl_magnetization(Spin_Config)/MedPerBlock;
				mean_sqrd_mag[vez] += pow(this->Ttl_magnetization(Spin_Config), 2)/MedPerBlock;
				mean_qdrt_mag[vez] += pow(this->Ttl_magnetization(Spin_Config), 4)/MedPerBlock;

				med_en += this->Lattice_Energy(Spin_Config)/MedPerBlock;
				med_mg += this->Ttl_magnetization(Spin_Config)/MedPerBlock;
				med_en_sq += pow(this->Lattice_Energy(Spin_Config), 2)/MedPerBlock;
				med_mg_sq += pow(this->Ttl_magnetization(Spin_Config), 2)/MedPerBlock;
				med_mg_qd += pow(this->Ttl_magnetization(Spin_Config), 4)/MedPerBlock;

				this->MetropolisStep(indx);
			}

			cv_proc[vez] = (med_en_sq - pow(med_en, 2))/pow(Temp, 2);
			sus_proc[vez] = (med_mg_sq - pow(med_mg, 2))/Temp;
			binder_ratio = (med_mg_qd)/pow(med_mg_sq, 2);
			binder_cumulant[vez] = (3.0/2.0)*(1.0 - binder_ratio/(3.0));
		}
	}

	void Parallel_Tempering(int ind_ant, int ind_nxt){
		int* ant = big_spin[ind_ant];
		int* nxt = big_spin[ind_nxt];
		double DS = (1/Temperatures[ind_nxt] - 1/Temperatures[ind_ant])*(this->Lattice_Energy(ant) - this->Lattice_Energy(nxt));

		tries += 1;

		if(DS <= 0){
			this->troca(ind_ant, ind_nxt);
			got += 1;
		}

		else{
			if(this->r_rand() < exp(-DS)){
				this->troca(ind_ant, ind_nxt);
				got += 1;
			}
		}
	}

	void troca(int ind_ant, int ind_nxt){
		int* ant_cp = new int[int_sites_qtd];
		int* nxt_cp = new int[int_sites_qtd];

		for(int i = 0; i < int_sites_qtd; i++){
			ant_cp[i] = big_spin[ind_ant][i];
			nxt_cp[i] = big_spin[ind_nxt][i];
		}

		for(int i = 0; i < int_sites_qtd; i++){
			big_spin[ind_ant][i] = nxt_cp[i];
			big_spin[ind_nxt][i] = ant_cp[i];
		}

		delete[] ant_cp;
		delete[] nxt_cp;
	}	

	double energy(){
		return Media(en_proc, blocks)/int_sites_qtd;
	}
	
	double senergy(){
		return DesvPadMed(en_proc, blocks)/int_sites_qtd;
	}
	
	double magnet(){
		return Media(mg_proc, blocks)/int_sites_qtd;
	}
	
	double smagnet(){
		return DesvPadMed(mg_proc, blocks)/int_sites_qtd;
	}
	
	double cv(){
		return Media(cv_proc, blocks)/int_sites_qtd;
	}
	
	double scv(){
		return DesvPadMed(cv_proc, blocks)/int_sites_qtd;
	}
	
	double susc(){
		return Media(sus_proc, blocks)/int_sites_qtd;
	}
	
	double ssusc(){
		return DesvPadMed(sus_proc, blocks)/int_sites_qtd;
	}
	
	double BinderRatio(){ //Calcula o BinderRatio da rede

		double meanSqMag = Media(mean_sqrd_mag, blocks)/int_sites_qtd;
		double meanQqMag = Media(mean_qdrt_mag, blocks)/int_sites_qtd;
		
		return meanQqMag/pow(meanSqMag, 2);
	}
	
	double BinderCumulant(){ //Calcula o BinderCumulant da rede
		return Media(binder_cumulant, blocks);
	}
	
	double sBinderCumulant(){ // O erro do Cumulante de Binder
		return DesvPadMed(binder_cumulant, blocks);
	}

};


void Parallel_exporter(std::string ComprimentoRede, std::ofstream& file){
	Lattice Redes(0.1, 2, 0.5, 3, 0.5, 5, 0.5, ComprimentoRede);
	int pontos = Redes.get_qtd_redes();

	double* values = new double[10];

	for(int i = 0; i < pontos; i++){
		Redes.get_mean_values(i, values);

		file << Redes.get_Temp(i) << "," << values[0] << "," << values[1] << ",";
		file << values[2] << "," << values[3] << ",";
		file << values[4] << "," << values[5] << ",";
		file << values[6] << "," << values[7] << "," << values[8] << "," << 0 << std::endl;
	}

	delete[] values;
}

int main(int argc, char** argv){

	std::string ComprimentoRede;

	if(argc >= 2){
		std::cout << "Comprimento da rede definido como: " << sites[std::stoi(argv[1])] << std::endl;
		ComprimentoRede = sites[std::stoi(argv[1])];
		
		OutFilePath = Data_out_Path + sites[std::stoi(argv[1])] + ".csv";
	}
	else {
		std::cout << "Nenhum comprimento de rede definido. Utilizando o comprimento padrão 239.";
		ComprimentoRede = "239";
		
		OutFilePath = Data_out_Path + "239" + "_par_temp.csv";
	}

	std::ofstream myfile;
	myfile.open(OutFilePath);

	Parallel_exporter(ComprimentoRede, myfile);

	myfile.close();

	std::cout << "Programa finalizado com sucesso!" << '\n';
	
	return 0;
}
