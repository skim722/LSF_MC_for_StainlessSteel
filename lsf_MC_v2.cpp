// LSF  Model
#include <string>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

//GLOBAL VARIABLES
string material;
int N=1;                          // number of spins
double T;                       // temperature
double H=0.0;                       //magnetic field
const int E_V_T0=-636.82394252;   // unit [eV]
int n_Fe, n_Cr, n_Ni; //42 
int cr_i, fe_i, ni_i;
double c_Fe=double(n_Fe)/double(N); double c_Cr=double(n_Cr)/double(N); double c_Ni=double(n_Ni)/double(N);
double *m_Cr=new double[n_Cr];  // Magnetic moments
double *m_Fe =new double[n_Fe]; //dynamic array
double *m_Ni=new double[n_Ni];
double *j_Cr=new double[n_Cr];  //J parameters
double *j_Fe=new double[n_Fe];
double *j_Ni=new double[n_Ni];
double** Cr_coeff= new double*[9];
double** Fe_coeff= new double*[42];
double** Ni_coeff= new double*[13];
int steps=0;
//FUNCTIONS
double Random_m(double min, double max);
double J_Fe(double m_Fe, int i);
double J_Cr(double m_Cr, int i );
double J_Ni(double m_Ni, int i);
void initialize ();
double energyPerComp( int element);
double energy();
double magPerComp( int element);
double magnetization();
bool MetropolisStep();
double oneMonteCarloStepPerSpin ( );
void ReadFile(string filename, double **Cr_coeff, double **Fe_coeff, double **Ni_coeff);
void PrintCoeff();
void DeleteCoeff();

void ReadFile(string filename, double **Cr_coeff, double **Fe_coeff, double **Ni_coeff){
     for(int i = 0 ; i < 9 ; i++){
        Cr_coeff[i] = new double[6];
     }
     for(int j = 0; j < 42 ; j++){
        Fe_coeff[j] = new double[6];
     }
     for(int k = 0; k < 13; k++){
        Ni_coeff[k] = new double[6];
     }
     ifstream file(filename.c_str());
     string line, name, temp;
     int number = 0;
     if(file.is_open()){
       while(!file.eof()){
           getline(file,line);
           for (int i = 0 ; i < 8 ; i++){
                   size_t index=line.find(" ");
                   temp = line.substr(0,index); // 0 to index
                   line=line.substr(index+1); //index+1 to end
                   if(i==0) name = temp.c_str();
                   else if(i==1) number = atoi(temp.c_str());
                   else{
                       if(name=="Cr") Cr_coeff[number-1][i-2]=atof(temp.c_str());
                       else if(name=="Fe") Fe_coeff[number-1][i-2]=atof(temp.c_str());
                       else if(name=="Ni") Ni_coeff[number-1][i-2]=atof(temp.c_str());
                  }
             }
        }
    }
}

void PrintCoeff(){
    cout<<" COEFFICENTS in OUR MEMORY \n";
    for(int i=0; i<9; i++){
         cout<<"Cr"<<i+1<<" ";
         for(int j=0;j<6;j++) cout<<Cr_coeff[i][j]<<"\t";
         cout<<"\n";
    }
    for(int j=0; j<42;j++){
        cout<<"Fe"<<j+1<<" ";
        for(int jj=0;jj<6;jj++) cout<<Fe_coeff[j][jj]<<"\t";
        cout<<"\n";
    }
    for(int i=0; i<13;i++){
        cout<<"Ni"<<i+1<<" ";
        for(int j=0;j<6;j++) cout<<Ni_coeff[i][j]<<"\t";
        cout<<"\n";
    }
}

void DeleteCoeff(){
     for(int i = 0 ; i < 9 ; i++){
        delete [] Cr_coeff[i];
     }
     delete Cr_coeff;
     for(int j = 0; j < 42 ; j++){
        delete [] Fe_coeff[j];
     }
     delete Fe_coeff;
     for(int k = 0; k < 13; k++){
        delete [] Ni_coeff[k];
     }
     delete Ni_coeff;
}


double Random_m(double min, double max){
	double r = (double)rand() / (double)RAND_MAX;
	return min + r * (max - min);
}

double J_Fe_pp(double m_Fe){   
       double fe_p1,fe_p2,fe_p3,fe_p4,fe_p5,fe_p6;
       fe_p1 =    0.002009;
       fe_p2 =   -0.002029;
       fe_p3 =     0.00919;
       fe_p4 =    -0.02375;
       fe_p5 =   0.0008004;
       fe_p6 =  -3.251e-05; 
return pow((fe_p1*pow(m_Fe,5) + fe_p2*pow(m_Fe,4) + fe_p3*pow(m_Fe,3) + fe_p4*pow(m_Fe,2) + fe_p5*pow(m_Fe,1)+fe_p6),1);
}


double J_Ni_pp(double m_Ni){
    double p1,p2,p3,p4;
    p1 =      0.0254;
    p2 =   -0.004665;
    p3 =     0.00988;
    p4 =   2.916e-05;
    return pow((p1*pow(m_Ni,3) + p2*pow(m_Ni,2) + p3*pow(m_Ni,1)+ p4),1);


}

double J_Cr_pp(double m_Cr){
    double cr_p1,cr_p2,cr_p3,cr_p4;
    cr_p1 =    0.003225;
    cr_p2 =    0.009891;
    cr_p3 =     0.00212;
    cr_p4 =    0.000189;
   return  pow((cr_p1*pow(m_Cr,3) + cr_p2*pow(m_Cr,2) + cr_p3*pow(m_Cr,1)+cr_p4),1);
}

double J_Fe(double m_Fe, int i ){   
    return (Fe_coeff[i][0]*pow(m_Fe,5) + Fe_coeff[i][1]*pow(m_Fe,4) + Fe_coeff[i][2]*pow(m_Fe,3) + Fe_coeff[i][3]*pow(m_Fe,2) + Fe_coeff[i][4]*pow(m_Fe,1)+Fe_coeff[i][5]);
}

/*
double J_Fe(double m_Fe,int i){
    return fe1[i]*pow(m_Fe,2);
} //TEST FOR EQUIPARTITION THEOREM
*/

double J_Ni(double m_Ni, int i){
    return (Ni_coeff[i][0]*pow(m_Ni,5) + Ni_coeff[i][1]*pow(m_Ni,4) + Ni_coeff[i][2]*pow(m_Ni,3) + Ni_coeff[i][3]*pow(m_Ni,2) + Ni_coeff[i][4]*pow(m_Ni,1)+Ni_coeff[i][5]); 
}

double J_Cr(double m_Cr, int i){
    return (Cr_coeff[i][0]*pow(m_Cr,5) + Cr_coeff[i][1]*pow(m_Cr,4) + Cr_coeff[i][2]*pow(m_Cr,3) + Cr_coeff[i][3]*pow(m_Cr,2) + Cr_coeff[i][4]*pow(m_Cr,1)+Cr_coeff[i][5]);
}

//Unit of Energy [eV]
//Unit of T [eV]=[11.604*10^3 Kelvin]

void initialize ( ) {
    for (int i=0; i<n_Cr;i++){
        m_Cr[i]=Random_m(0,2.903); //2.903 
        j_Cr[i]=J_Cr(abs(m_Cr[i]), cr_i);
        }                             
    for (int j=0; j<n_Fe;j++){
        m_Fe[j]=Random_m(0,3.218); //3.218
        j_Fe[j]=J_Fe(abs(m_Fe[j]),fe_i);
        }
    for (int k=0; k<n_Ni;k++){
        m_Ni[k]=Random_m(0,1.303);  //1.303 
        j_Ni[k]=J_Ni(abs(m_Ni[k]), ni_i);
        }   
    steps = 0;

}

double energyPerComp( int element){
     double e=0.0;
     if (element==1) //Cr
     { 
         for (int i=0; i<n_Cr;i++)
              e+=j_Cr[i]; 
     }
     else if (element==2) //Fe
     {    
         for (int j=0; j<n_Fe;j++)
              e+=j_Fe[j];  
     }
     else if (element==3) //Ni
     {    
         for (int k=0; k<n_Ni;k++)
              e+=j_Ni[k];  
     }
    return e;
}

double energy(){
    double etot=0;
    for (int ii=0;ii<3;ii++)
       etot+=energyPerComp(ii+1);
    return etot;
}

double magPerComp( int element){
     double m=0.0; double N_ele;
     
     if (element==1) //Cr
     {  
         N_ele=n_Cr; 
         for (int i=0; i<n_Cr;i++)
              m+=m_Cr[i];
     }
     else if (element==2) //Fe
     {   
         N_ele=n_Fe;
         for (int j=0; j<n_Fe;j++)
              m+=m_Fe[j];
     }
     else if (element==3) //Ni
     {   
         N_ele=n_Ni;
         for (int k=0; k<n_Ni;k++)
              m+=m_Ni[k];
     }
  return m/double(N_ele);
}

double RMS_magPerComp( int element){
      double m=0.0; double N_ele;
 
      if (element==1) //Cr
      {
          N_ele=n_Cr;
          m=0.0;
          for (int i=0; i<n_Cr;i++)
               m+=m_Cr[i]*m_Cr[i];
      }
      else if (element==2) //Fe
      {
          N_ele=n_Fe;
          m=0.0;
          for (int j=0; j<n_Fe;j++){
               m+=(abs(m_Fe[j]))*(abs(m_Fe[j])); 
              }
      }
      else if (element==3) //Ni
      {
          N_ele=n_Ni;
          for (int k=0; k<n_Ni;k++)
              m+= m_Ni[k]*m_Ni[k];
      }
    return m/N_ele;
 }


double magnetization(){
   return magPerComp(1)*c_Cr+magPerComp(2)*c_Fe+magPerComp(3)*c_Ni;
}

bool MetropolisStep () {
     double E_bf, E_af;
     double m_bf,j_bf,max_m,min_m;
     E_bf=energy();

     if (material=="Cr") //Cr
        {   
            m_bf=m_Cr[0]; j_bf=j_Cr[0]; //Save old m&j
            m_Cr[0]=Random_m(-2*2.903,2*2.903);      //Random_m(min_m,max_m);
            j_Cr[0]=J_Cr(abs(m_Cr[0]), cr_i);
        }
     else if (material=="Fe") //Fe
        {
            m_bf=m_Fe[0]; j_bf=j_Fe[0];
            m_Fe[0]=Random_m(-2*3.218,2*3.218);     //Random_m(min_m,max_m);
            j_Fe[0]=J_Fe(abs(m_Fe[0]), fe_i);
        }
     else if (material=="Ni") //Ni
        {
            m_bf=m_Ni[0]; j_bf=j_Ni[0];
            m_Ni[0]=Random_m(-2*1.303,2*1.303); //1.303 
            j_Ni[0]=J_Ni(abs(m_Ni[0]),ni_i);  
        }  
   
      E_af=energy(); //cout<<"Ebf  "<<E_bf<<"\t Eaf  "<<E_af<<endl;
      //ratio of Boltzmann factors
      long double ratio = -(E_af-E_bf)/(T+1e-9);
//     cout<<"\tE_bf\t"<<E_bf+583.3<<"\tE_af\t"<<E_af+583.3<<"\t-(E_af-E_bf)/(T+1e-9)\t"<<ratio<<"\tRATIO\t"<<expl(ratio);
      if (std_rand()<expl(ratio))
          {
          // cout<<"\tACCEPT!"<<endl;
           return true;
          }
      else
         {
          //  cout<<"\tREJECT!"<<endl;
           //Return to previous m
             if (material=="Cr"){ m_Cr[0]=m_bf; j_Cr[0]=j_bf; }
             else if (material=="Fe"){m_Fe[0]=m_bf;j_Fe[0]=j_bf;}
             else if(material=="Ni"){m_Ni[0]=m_bf; j_Ni[0]=j_bf;} 
           return false; //reject
         }     
}

double oneMonteCarloStepPerSpin ( ) {
     double acceptanceRatio;
     int accepts = 0;
     for (int i = 0; i < N; i++)
          if (MetropolisStep())
               ++accepts;
     acceptanceRatio = accepts/double(N);
     ++steps;
     return acceptanceRatio;
}

int main (int argc, char *argv[]) {
 ReadFile("coeff.txt", Cr_coeff, Fe_coeff, Ni_coeff);
 PrintCoeff();
 cout.precision(4);                        //set precision
 cout.setf(ios::fixed);
 cout << " Simplified LSF- Metropolis simulation\n"
      << " ---------------------------------------------------\n";
 cout << " Enter number of Monte Carlo steps: ";
 int MCSteps; cin >> MCSteps;
 double T_high;
 cout<< "Enter the highest tempearture : "; cin>>T_high; //Get T as Kelvin
 ofstream file("summation_of_J.data");
 ofstream afile("E_V_mT.data"); 
 file << " N (# of atoms): " << N <<'\t' << " H(Magnetic Field): " << H <<'\t'<<" MCSteps: "<< MCSteps << endl;
 size_t size=64;
 std::vector <double> J_T0(size), J_i(size),J_ave(size), J_T0_ave(size);
 for(int i=0; i<size; ++i) { J_T0[i]=0.0; J_i[i]=0.0; J_ave[i]=0.0; J_T0_ave[i]=0.0;}
 cr_i=0; ni_i=0; fe_i=0; 
 for (int index=0; index<64; index++){
    if(index<9){ material="Cr"; n_Cr=1; n_Fe=0; n_Ni=0; cr_i=index; fe_i=-1; ni_i=-1;}
    else if( index <51){ material="Fe"; n_Cr=0; n_Fe=1; n_Ni=0; cr_i = -1; fe_i=index-9; ni_i=-1;}
    else{material="Ni"; n_Cr=0; n_Fe=0; n_Ni=1; cr_i=-1; fe_i=-1; ni_i=index-51;}
    double sum_J_i=0.0;
    for(int t=0; t<(T_high/10); t++){
         T=double(t*10+0.1)*(8.6173324*0.00001); //Unit of temperature i K->eV 
         initialize();
         int thermSteps = int(0.2 * MCSteps); 
         for (int s = 0; s < thermSteps; s++) oneMonteCarloStepPerSpin();
         double mAv = 0, m2Av = 0, eAv = 0, e2Av = 0,RMS_mAv_Cr=0, RMS_mAv_Fe=0, RMS_mAv_Ni=0, mAv_Cr = 0, m2Av_Cr = 0, eAv_Cr = 0, e2Av_Cr = 0,mAv_Fe = 0, m2Av_Fe = 0, eAv_Fe = 0, e2Av_Fe = 0,mAv_Ni = 0, m2Av_Ni = 0, eAv_Ni = 0, e2Av_Ni = 0;
         J_ave[index]=0.0;     
         for (int s = 0; s < MCSteps; s++) {
              for(int inner_loop=0; inner_loop<100; inner_loop++) oneMonteCarloStepPerSpin();
              double m = magnetization(); double m_Cr=magPerComp(1); double m_Fe=magPerComp(2); double m_Ni=magPerComp(3);
              double RMS_m_Cr=RMS_magPerComp(1); double RMS_m_Fe=RMS_magPerComp(2); double RMS_m_Ni=RMS_magPerComp(3);
              double e = energy(); double e_Cr=energyPerComp(1); double e_Fe=energyPerComp(2); double e_Ni=energyPerComp(3);
              if(index<9) J_ave[index]+=J_Cr(abs(m_Cr),cr_i);
              else if(index<51) {J_ave[index]+=J_Fe(abs(m_Fe), fe_i); } 
              else J_ave[index]+=J_Ni(abs(m_Ni),ni_i);
              mAv += m; m2Av += m * m;  eAv += e; e2Av += e * e;
              RMS_mAv_Cr += RMS_m_Cr; RMS_mAv_Fe+=RMS_m_Fe; RMS_mAv_Ni+=RMS_m_Ni;
              mAv_Cr += abs(m_Cr); m2Av_Cr += m_Cr * m_Cr;  eAv_Cr += e_Cr; e2Av_Cr += e_Cr * e_Cr;
              mAv_Fe += abs(m_Fe); m2Av_Fe += m_Fe * m_Fe;  eAv_Fe += e_Fe; e2Av_Fe += e_Fe * e_Fe;
              mAv_Ni += abs(m_Ni); m2Av_Ni += m_Ni * m_Ni;  eAv_Ni += e_Ni; e2Av_Ni += e_Ni * e_Ni;
         }
         RMS_mAv_Cr = sqrt(RMS_mAv_Cr/MCSteps); RMS_mAv_Fe=sqrt( RMS_mAv_Fe/MCSteps); RMS_mAv_Ni= sqrt(RMS_mAv_Ni/MCSteps);
         mAv /= MCSteps; m2Av /= MCSteps; eAv /= MCSteps; e2Av /= MCSteps;
         mAv_Cr /= MCSteps; m2Av_Cr /= MCSteps; eAv_Cr /= MCSteps; e2Av_Cr /= MCSteps;
         mAv_Fe /= MCSteps; m2Av_Fe /= MCSteps; eAv_Fe /= MCSteps; e2Av_Fe /= MCSteps;
         mAv_Ni /= MCSteps; m2Av_Ni /= MCSteps; eAv_Ni /= MCSteps; e2Av_Ni /= MCSteps;
         J_ave[index]/=MCSteps; 
         if (t==0){
            if(index<9) { J_T0[index]=J_Cr(RMS_mAv_Cr, cr_i); J_T0_ave[index]=J_ave[index]; }
            else if(index<51) {J_T0[index]=J_Fe(RMS_mAv_Fe, fe_i); J_T0_ave[index]= J_ave[index];}
            else {J_T0[index]=J_Ni(RMS_mAv_Ni, ni_i); J_T0_ave[index]=J_ave[index];}      
         }
         if(index<9){
            file <<index+1<<" "<<material<<cr_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Cr<<"\t"<<RMS_mAv_Cr << "\t"<<J_Cr(RMS_mAv_Cr,cr_i)<<"\t"<<J_T0[index]<<"\t"<<J_Cr(RMS_mAv_Cr, cr_i)-J_T0[index]<<"\t"<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
            cout<<index+1<<" "<<material<<cr_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Cr<<"\t"<<RMS_mAv_Cr << "\t"<<J_Cr(RMS_mAv_Cr,cr_i)<<"\t"<<J_T0[index]<<"\t"<<J_Cr(RMS_mAv_Cr, cr_i)-J_T0[index]<<"\t"<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
            sum_J_i+=(J_Cr(RMS_mAv_Cr, cr_i)-J_T0[index]);
         }
         else if(index<51){
            file <<index+1<<" "<<material<<fe_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Fe<<"\t"<<RMS_mAv_Fe << "\t"<<J_Fe(RMS_mAv_Fe,fe_i)<<"\t"<<J_T0[index]<<"\t"<<J_Fe(RMS_mAv_Fe, fe_i)-J_T0[index]<<"\t"<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
            cout <<index+1<<" "<<material<<fe_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Fe<<"\t"<<RMS_mAv_Fe << "\t"<<J_Fe(RMS_mAv_Fe,fe_i)<<"\t"<<J_T0[index]<<"\t"<<J_Fe(RMS_mAv_Fe, fe_i)-J_T0[index]<<"\t"<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
            sum_J_i+=J_Fe(RMS_mAv_Fe, fe_i)-J_T0[index];
         }
         else{
            file <<index+1<<" "<<material<<ni_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Ni<<"\t"<<RMS_mAv_Ni << "\t"<<J_Ni(RMS_mAv_Ni,ni_i)<<"\t"<<J_T0[index]<<"\t"<<J_Ni(RMS_mAv_Ni, ni_i)-J_T0[index]<<"\t"<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
            cout<<index+1<<" "<<material<<ni_i+1<<" "<<T*11.60496*1000 <<"\t"<<mAv_Ni<<"\t"<<RMS_mAv_Ni << "\t"<<J_Ni(RMS_mAv_Ni,ni_i)<<"\t" <<J_T0[index]<<"\t"<<J_Ni(RMS_mAv_Ni, ni_i)-J_T0[index]<<J_ave[index]<<"\t"<<J_T0_ave[index]<<"\t"<<J_ave[index]-J_T0_ave[index]<<"\t"<<(J_ave[index]-J_T0_ave[index])<<"\t"<<T<<endl;
         }
        
      }
      file <<"\t"<<"T=\t"<<T*11.60496*1000<<"\tSum_i(J(m_i, T))\t"<<sum_J_i<<"\tE(V,{m(T)}=\t"<<E_V_T0 + sum_J_i<<"\n    "<<endl;
      cout<<"\t"<<"T=\t"<<T*11.60496*1000<<"\tSum_i(J(m_i, T))\t"<<sum_J_i<<"\tE(V,{m(T)}=\t"<<E_V_T0 + sum_J_i<<"\n    "<<endl;
      afile<<"\t"<<"T=\t"<<T*11.60496*1000<<"\tSum_i(J(m_i, T))\t"<<sum_J_i<<"\tE(V,{m(T)}=\t"<<E_V_T0 + sum_J_i<<"\n    "<<endl; }
      file.close();
      afile.close();
}
