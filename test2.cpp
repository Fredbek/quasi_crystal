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

class Net{
	
private:
	
	std::string NetLength = "239";
	int NumLat = std::stoi(NetLength);
	
	std::string data_path = lattices_path_init + NetLength + "/" + data_arch_name + NetLength + extension;

	int*** big_spin = NULL;
	int** big_bf_spin = NULL;
	double* big_energy = NULL;
	double* Temperatures = NULL;
	double* probabilities = NULL;

	int* phase_net = NULL;

	void(Net::*Term_Algthm)(int);

	double* nets_tries = NULL;
	double* nets_got = NULL;

	int qtd_redes = qtd;
	double tries = 0; //Tentativas do parallel tempering
	double got = 0; //Vezes que rodou o parallel tempering

	int** Rede = NULL;
	bool* bf_Rede = NULL; // buffer que vê quais spins já foram adicionados no clsuter
	int** Viz = NULL;//Viz[NumLat][4];
	
	//double data_crua[Med][3];// Pega os dados das medidas na ordem (en, en^2, m)
	double en_proc[blocks], mg_proc[blocks], cv_proc[blocks], sus_proc[blocks];
	
	//double raw_data[Med][2]; //Data crua com entradas (m^2, m^4)
	double mean_sqrd_mag[blocks]; //blocos com medidas das grandezas desejadas
	double mean_qdrt_mag[blocks];
	double binder_cumulant[blocks];
	
	std::vector<int> bf_cluster; //Buffer do cluster que será girado
	
	double T; //Temperatura da rede
	
	double Padd = 1 - exp(-2*abs(J)/T); // Probabilidade de adicionar um spin na rede
	
	/// Definições do rng
	std::uniform_real_distribution<double> dist{ 0, 1.0 };
	std::uniform_int_distribution<int> dice{0, (NumLat - 1)};
	
	int WolffTermSteps = floor(std::stoi(NetLength)); //Quantidade de passos para a termalização
	int MetTermSteps = 50*std::stoi(NetLength);
	int TermSteps;
	
	bool bool_Tc;
	
public:

	Net(){}

	Net(std::string Comp)
		:NetLength(Comp)
	{
		Temperatures = new double[qtd_redes]; // Inicia dinamicamente as temperaturas
		big_spin = new int**[qtd_redes]; // Inicia dinâmicamente as redes
		bf_Rede = new bool[NumLat]; //Inicia dinâmicamente os buffers das redes
		big_energy = new double[qtd_redes]; //Matriz que guarda a energia pras temperaturas

		phase_net = new int[NumLat];

		nets_tries = new double[qtd_redes];
		nets_got = new double[qtd_redes];

		probabilities = new double[qtd_redes];

		double init_temp = h_temp;

		for(int i = 0; i < qtd_redes; i++){
			Temperatures[i] = init_temp;

			probabilities[i] = 1 - exp(-2*abs(J)/Temperatures[i]);

			big_spin[i] = new int*[NumLat];

			for(int j = 0; j < NumLat; j++){
				big_spin[i][j] = new int[2];
				big_spin[i][j][0] = 0;
			}

			init_temp += h_temp;
		}
		
		Viz = new int*[NumLat]; //Inicia dinamicamente a tabela de vizinhos
		for(int i = 0; i < NumLat; i++){
			Viz[i] = new int[8];
			phase_net[i] = 0;
		}

		if(T < 2 or T > 2.44){//Seta a quantidade de passos de termalização
			bool_Tc = false;
			TermSteps = MetTermSteps;
			Term_Algthm = &Net::MetropolisStep;
		}
		else{
			bool_Tc = true;
			TermSteps = WolffTermSteps;
			Term_Algthm = &Net::WollfStep;
		}

		this->InitNets();
		this->Reset_bf_Rede();
		this->getViz();
		this->getPhase(0, 1);
		this->Parallel_Termalization();
	}

	Net(double Temp, std::string Comp)
		:T(Temp), NetLength(Comp)
	{
		
		Rede = new int*[NumLat]; //Inicia dinamicamente a rede
		
		bf_Rede = new bool[NumLat]; //Inicia dinamicamente o buffer da rede
		
		Viz = new int*[NumLat]; //Inicia dinamicamente a tabela de vizinhos
		for(int i = 0; i < NumLat; i++){
			Rede[i] = new int[2];
			Viz[i] = new int[8];
			Rede[i][0] = 0;
		}
		
		if(T < 2 or T > 2.44){//Seta a quantidade de passos de termalização
			bool_Tc = false;
			TermSteps = MetTermSteps;
			Term_Algthm = &Net::MetropolisStep;
		}
		else{
			bool_Tc = true;
			TermSteps = WolffTermSteps;
			Term_Algthm = &Net::WollfStep;
		}
		
		this->InitNet();
		this->getViz();
		this->getPhase(0, 1);
		this->Termalization();
		this->getProcData(Rede);
	}
	
	~Net(){

		if(big_spin != NULL){
			for(int i = 0; i < qtd_redes; i++){
				for(int j = 0; j < NumLat; j++){
					delete[] big_spin[i][j];
				}

				delete[] big_spin[i];
			}

			delete[] big_spin;
		
			delete[] Temperatures;
			delete[] big_energy;
		}

		if(Rede != NULL){
			delete[] Rede;
		}
		if(bf_Rede != NULL){
			delete[] bf_Rede;
		}
		
		if(Viz != NULL){
			for(int i = 0; i < NumLat; i++){
				delete[] Viz[i];
			}
			delete[] Viz;	
		}

		delete[] phase_net;
	}
	
	void getPhase(int i, int initial){
		if(phase_net[i] == 0){
			phase_net[i] = initial;
			for(int ind = 0; ind < 8; ind++){
				if(Viz[i][ind] >= 0){
					this->getPhase(Viz[i][ind], -initial);
				}
			}
		}
	}
	
	void getPhases(int i, int initial){

		if(big_spin[0][i][0] == 0){
			big_spin[0][i][0] = initial;
			for(int ind = 0; ind < 8; i++){
				if(Viz[i][ind] >= 0 and big_spin[0][Viz[i][ind]][0] == 0){
					this->getPhases(Viz[i][ind], -initial);
				}
			}
		}
	}

	double i_rand() { ///Função que gera um inteiro aleatório entre 0 e NumLat
		return dice(eng);
	}
	
	double r_rand() { ///Função que gera uma real alatório entre 0 e 1
		return dist(eng);
	}
	
	void InitNet(){ //Inicia a rede colocando todos os spins para cima
	
		for(int i = 0; i < NumLat; i++){
			Rede[i][1] = 1;
			bf_Rede[i] = false;
		}
	}

	void InitNets(){ //Inicia a rede colocando todos os spins para cima
		for(int k = 0; k < qtd_redes; k++){
			for(int i = 0; i < NumLat; i++){
				big_spin[k][i][1] = 1;
			}
		}
	}

	void Reset_bf_Rede(){ //reinicia o buffer da rede
		for(int i = 0; i < NumLat; i++){
			bf_Rede[i] = false;
		}
	}

	void getViz(){ //Lê a tabela de vizinhos
	
		std::string VizPath = lattices_path_init + NetLength + "/appxt_" + NetLength + "_connect.dat";

		std::ifstream Data;

		Data.open(VizPath);

		for(int i = 0; i < NumLat; i++){
			for(int j = 0; j < 8; j++){
				Data >> Viz[i][j];
				Viz[i][j] = Viz[i][j] - 1;
			}
		}

		Data.close();
	}

	void addCluster(int ind){
		
		bf_cluster.push_back(ind);//Adiciona o spin no cluster
		bf_Rede[ind] = true;
		
		for(int i = 0; i < 8; i++){ //Percorre os vizinhos
			if(Viz[ind][i] >= 0){
				if(Padd >= this->r_rand() and bf_Rede[Viz[ind][i]] == false and phase_net[ind]*Rede[ind][1] == phase_net[Viz[ind][i]]*Rede[Viz[ind][i]][1]){ //Adiciona caso já não esteja com prob Padd 
					
					this->addCluster(Viz[ind][i]);
				}
			}
		}
	}
	
	void WollfStep(){
		int ind = this->i_rand(); //Sorteia um spin aleatório
		
		this->addCluster(ind); //Adiciona ele e os vizinhos no buffer de acordo com a probabilidade
		
		for(auto i : bf_cluster){ //Gira todos os spins no cluster
			Rede[i][1] = -Rede[i][1];
		}
		
		this->Reset_bf_Rede(); //Reseta o cluster e o buffer da rede
		bf_cluster.clear();
		
	}
	
	void MetropolisStep(){ //Definição do passo de metropolis #1 Varredura
		for(int i = 0; i < NumLat; i++){ //Varre todos os sítios da rede
			double DE = -2*this->LatticeEnergy(Rede, i); //Variação de energia ao girar o spin
			
			if(DE <= 0){
				Rede[i][1] = -Rede[i][1];
			}
			else{
				double z = exp(-DE/T);
				if(z >= r_rand()){
					Rede[i][1] = -Rede[i][1];
				}
			}
		}
	}
	
	void Termalization(){
		for(int i = 0; i < TermSteps; i++){

			double j = 10;

			if(bool_Tc == true){ //Se estiver na regiãao de transição usa o algoritmo de wolff
				if(j = 0){
					this->WollfStep();
					j = 11;
					}
				else{
					this->MetropolisStep();
				}
			}
			else{ //fora da região usa o de wolff
				this->MetropolisStep();
			}

			j -= 1;
		}
	}

	double SpinPerLat(int** Spin_Config){ //Calcula a magnetização por spin da rede atual
		double soma = 0;
		
		for(int i = 0; i < NumLat; i++){
			soma += Spin_Config[i][1]*phase_net[i];
		}
		
		return abs(soma);
	}
	
	double LatticeEnergy(int** Spin_Config, int i){ //Energia associada à um 
		double soma = 0;

		for(int j = 0; j < 8; j++){
			if(Viz[i][j] >= 0){
				soma += -J*Spin_Config[i][1]*Spin_Config[Viz[i][j]][1];
			}
		}

		return soma;
	}
	
	double EnergyPerLat(int** Spin_Config){
		double soma = 0;
		
		for(int i = 0; i < NumLat; i++){
			soma += this->LatticeEnergy(Spin_Config, i);
		}
		
		return soma/2;
	}
	
	void getProcData(int** Spin_Config){

		double med_en = 0, med_en_sq = 0, med_mg = 0, med_mg_sq = 0, med_mg_qd = 0, binder_ratio = 0;

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
				en_proc[vez] += this->EnergyPerLat(Spin_Config)/MedPerBlock;
				mg_proc[vez] += this->SpinPerLat(Spin_Config)/MedPerBlock;
				mean_sqrd_mag[vez] += pow(this->SpinPerLat(Spin_Config), 2)/MedPerBlock;
				mean_qdrt_mag[vez] += pow(this->SpinPerLat(Spin_Config), 4)/MedPerBlock;

				med_en += this->EnergyPerLat(Spin_Config)/MedPerBlock;
				med_mg += this->SpinPerLat(Spin_Config)/MedPerBlock;
				med_en_sq += pow(this->EnergyPerLat(Spin_Config), 2)/MedPerBlock;
				med_mg_sq += pow(this->SpinPerLat(Spin_Config), 2)/MedPerBlock;
				med_mg_qd += pow(this->SpinPerLat(Spin_Config), 4)/MedPerBlock;

				if(bool_Tc == true){
					this->WollfStep();
				}
				else{
					this->MetropolisStep();
				}
			}

			cv_proc[vez] = (med_en_sq - pow(med_en, 2))/pow(T, 2);
			sus_proc[vez] = (med_mg_sq - pow(med_mg, 2))/T;
			binder_ratio = (med_mg_qd)/pow(med_mg_sq, 2);
			binder_cumulant[vez] = (3.0/2.0)*(1.0 - binder_ratio/(3.0));
		}
	}
	
	double energy(){
		return Media(en_proc, blocks)/NumLat;
	}
	
	double senergy(){
		return DesvPadMed(en_proc, blocks)/NumLat;
	}
	
	double magnet(){
		return Media(mg_proc, blocks)/NumLat;
	}
	
	double smagnet(){
		return DesvPadMed(mg_proc, blocks)/NumLat;
	}
	
	double cv(){
		return Media(cv_proc, blocks)/NumLat;
	}
	
	double scv(){
		return DesvPadMed(cv_proc, blocks)/NumLat;
	}
	
	double susc(){
		return Media(sus_proc, blocks)/NumLat;
	}
	
	double ssusc(){
		return DesvPadMed(sus_proc, blocks)/NumLat;
	}
	
	double BinderRatio(){ //Calcula o BinderRatio da rede

		double meanSqMag = Media(mean_sqrd_mag, blocks)/NumLat;
		double meanQqMag = Media(mean_qdrt_mag, blocks)/NumLat;
		
		return meanQqMag/pow(meanSqMag, 2);
	}
	
	double BinderCumulant(){ //Calcula o BinderCumulant da rede
		return Media(binder_cumulant, blocks);
	}
	
	double sBinderCumulant(){ // O erro do Cumulante de Binder
		return DesvPadMed(binder_cumulant, blocks);
	}

	//////////////////////// Definições para o Parallel Tempering ////////////////////////////////

	void MetropolisStep(int ind){ //Definição do passo de metropolis #1 Varredura
		
		int** Spin_Config = big_spin[ind];
		double Temp = Temperatures[ind];
		
		for(int i = 0; i < NumLat; i++){ //Varre todos os sítios da rede
			double DE = -2*this->LatticeEnergy(Spin_Config, i); //Variação de energia ao girar o spin
			
			if(DE <= 0){
				Spin_Config[i][1] = -Spin_Config[i][1];
			}
			else{
				double z = exp(-DE/Temp);
				if(z >= r_rand()){
					Spin_Config[i][1] = -Spin_Config[i][1];
				}
			}
		}
	}

	void Parallel_Termalization(){

		for(int cont = 0; cont < TermSteps/5; cont++){
			this->TermStep();
		}

		std::cout << "Parallel Tempering Acceptance Ratio: " << 100.0*(got/tries) << "%" << std::endl;
	}

	void TermStep(){ //For now just with the metropolis algorithm
		for(int cont = 0; cont < 5; cont++){
			for(int indx = 0; indx < qtd_redes; indx++){

				if(Temperatures[indx] < 2.0 or Temperatures[indx] > 2.44){//Seta a quantidade de passos de termalização
					bool_Tc = false;
					TermSteps = MetTermSteps;
					Term_Algthm = &Net::MetropolisStep;
				}
				else{
					bool_Tc = true;
					TermSteps = WolffTermSteps;
					Term_Algthm = &Net::WollfStep;
				}

				(this->*Term_Algthm)(indx); //Roda o algoritmo atual

				void(Net::*test_wolff)(int) = &Net::WollfStep;

				if(Term_Algthm == test_wolff){ //Se o algoritmo atual for o de wolff
					this->MetropolisStep(indx); // Roda uma varredura de Metropolis
					print("1");
				}
			}
		}

		for(int i = qtd_redes - 1; i >= 1; i -= 1){
			int rede_ant = i - 1;
			this->Parallel_Tempering(rede_ant, i);
		}
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

		int** Spin_Config = big_spin[indx];

		double Temp = Temperatures[indx];

		double med_en = 0, med_en_sq = 0, med_mg = 0, med_mg_sq = 0, med_mg_qd = 0, binder_ratio = 0;

		if(Temperatures[indx] < 2.0 or Temperatures[indx] > 2.44){//Seta a quantidade de passos de termalização
			bool_Tc = false;
			TermSteps = MetTermSteps;
			Term_Algthm = &Net::MetropolisStep;
		}
		else{
			bool_Tc = true;
			TermSteps = WolffTermSteps;
			Term_Algthm = &Net::WollfStep;
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
				en_proc[vez] += this->EnergyPerLat(Spin_Config)/MedPerBlock;
				mg_proc[vez] += this->SpinPerLat(Spin_Config)/MedPerBlock;
				mean_sqrd_mag[vez] += pow(this->SpinPerLat(Spin_Config), 2)/MedPerBlock;
				mean_qdrt_mag[vez] += pow(this->SpinPerLat(Spin_Config), 4)/MedPerBlock;

				med_en += this->EnergyPerLat(Spin_Config)/MedPerBlock;
				med_mg += this->SpinPerLat(Spin_Config)/MedPerBlock;
				med_en_sq += pow(this->EnergyPerLat(Spin_Config), 2)/MedPerBlock;
				med_mg_sq += pow(this->SpinPerLat(Spin_Config), 2)/MedPerBlock;
				med_mg_qd += pow(this->SpinPerLat(Spin_Config), 4)/MedPerBlock;

				(this->*Term_Algthm)(indx);
				void(Net::*test_wolff)(int) = &Net::WollfStep;
				if(Term_Algthm == test_wolff){ //Se o algoritmo atual for o de wolff
					print("2");
				}
			}

			cv_proc[vez] = (med_en_sq - pow(med_en, 2))/pow(Temp, 2);
			sus_proc[vez] = (med_mg_sq - pow(med_mg, 2))/Temp;
			binder_ratio = (med_mg_qd)/pow(med_mg_sq, 2);
			binder_cumulant[vez] = (3.0/2.0)*(1.0 - binder_ratio/(3.0));
		}
	}

	void Parallel_Tempering(int ind_ant, int ind_nxt){
		int** ant = big_spin[ind_ant];
		int** nxt = big_spin[ind_nxt];
		double DS = (1/Temperatures[ind_nxt] - 1/Temperatures[ind_ant])*(this->EnergyPerLat(ant) - this->EnergyPerLat(nxt));

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
		int* ant_cp = new int[NumLat];
		int* nxt_cp = new int[NumLat];

		for(int i = 0; i < NumLat; i++){
			ant_cp[i] = big_spin[ind_ant][i][1];
			nxt_cp[i] = big_spin[ind_nxt][i][1];
		}

		for(int i = 0; i < NumLat; i++){
			big_spin[ind_ant][i][1] = nxt_cp[i];
			big_spin[ind_nxt][i][1] = ant_cp[i];
		}

		delete[] ant_cp;
		delete[] nxt_cp;
	}

	void addCluster(int ind_net, int ind){
		
		int** Spin_Config = big_spin[ind_net];

		double prob = 1 - exp(-2*abs(J)/Temperatures[ind_net]);//probabilities[ind_net];

		bf_cluster.push_back(ind);//Adiciona o spin no cluster
		bf_Rede[ind] = true;
		
		for(int i = 0; i < 8; i++){ //Percorre os vizinhos
			if(Viz[ind][i] >= 0){
				if(prob >= this->r_rand() and bf_Rede[Viz[ind][i]] == false and phase_net[ind]*Spin_Config[ind][1] == phase_net[Viz[ind][i]]*Spin_Config[Viz[ind][i]][1]){ //Adiciona caso já não esteja com prob Padd 
					this->addCluster(ind_net, Viz[ind][i]);
				}
			}
		}
	}
	
	void WollfStep(int ind_net){
		int ind = this->i_rand(); //Sorteia um spin aleatório
		int** Spin_Config = big_spin[ind_net];
		this->addCluster(ind_net, ind); //Adiciona ele e os vizinhos no buffer de acordo com a probabilidade

		for(auto i : bf_cluster){ //Gira todos os spins no cluster
			Spin_Config[i][1] = -Spin_Config[i][1];
		}
		
		this->Reset_bf_Rede(); //Reseta o cluster e o buffer da rede
		bf_cluster.clear();
		
	}
};

void exporter(std::string ComprimentoRede, std::ofstream& file, double T_in, double T_fin, double h){
	for(double i = T_in; i <= T_fin; i += h){
		Net Camp(i, ComprimentoRede);
		
		file << i << "," << Camp.energy() << "," << Camp.senergy() << ",";
		file << Camp.magnet() << "," << Camp.smagnet() << ",";
		file << Camp.cv() << "," << Camp.scv() << ",";
		file << Camp.susc() << "," << Camp.ssusc() << "," << Camp.BinderCumulant() << "," << 0 << std::endl;
		
		std::cout << (i/5)*100 << "%" << std::endl;
	}
}

void Parallel_exporter(std::string ComprimentoRede, std::ofstream& file){
	Net Redes(ComprimentoRede);
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
