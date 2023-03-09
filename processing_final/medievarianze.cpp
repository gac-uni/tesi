#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double ** medvar(int npassi, int nwalk, double tp1, string workdir, string label){
    double * medwalk;
    medwalk = new double[nwalk];
    double * varwalk;
    varwalk = new double[nwalk];
    double ** out = new double*[2];
    out[0]=medwalk;
    out[1]=varwalk;

    double tempvalori[npassi]={0};
    for(int i=0;i<40;i++){
        ifstream temp (workdir+"/mt_val_dim/"+label+"_"+to_string(i)+".dat");
        if(temp.is_open()){
            string line;
            int count=0;
            while(getline(temp,line)){
                tempvalori[count]=stod(line);
                medwalk[i]+=tempvalori[count];
                count++;
            }
            medwalk[i]/=count;
            temp.close();
        } 
        for(int j=0;j<npassi;j++){
            varwalk[i]+=(medwalk[i]-tempvalori[j])*(medwalk[i]-tempvalori[j]);
        }
        varwalk[i]*=(tp1/(npassi));
        varwalk[i]/=npassi-1;
    }
    return out;
}

double * medscema(int nelem, double * vmed, double * vvar){
    double * out = new double[2];
    for(int i=0;i<nelem;i++){
        out[0] += vmed[i];
        out[1] += vvar[i];
    }
    out[0] /= nelem;
    out[1] /= nelem;
    return out;
}

double * medpesata(int nelem, double * vmed, double * vvar){
    double * out = new double[2];
    for(int i=0;i<nelem;i++){
        out[0] += vmed[i]/vvar[i];
        out[1] += 1/vvar[i];
    }
    out[0] /= out[1];
    out[1] = 1/out[1];
    return out;
}

int main(int argc, char ** argv){
    string workdir;
    workdir = argv[1];
    double tp1ken;
    double tp1pot;
    double tp1stim1;
    double tp1stim2;
    int nwalk;
    int npassi;
    sscanf(argv[2], "%lf", &tp1ken);
    sscanf(argv[3], "%lf", &tp1pot);
    sscanf(argv[4], "%lf", &tp1stim1);
    sscanf(argv[5], "%lf", &tp1stim2);
    sscanf(argv[6], "%d", &nwalk);
    sscanf(argv[7], "%d", &npassi);
    // double tp1ken=36.649819;
    // double tp1pot=94.523932;
    // double tp1stim1=98.890308;
    // double tp1stim2=262.960069;
    double * medwalkken;
    double * varwalkken;
    double * medwalkpot;
    double * varwalkpot;
    double * medwalkstim1;
    double * varwalkstim1;
    double * medwalkstim2;
    double * varwalkstim2;
    double ** tempken=medvar(npassi, nwalk, tp1ken, workdir, "ken");
    medwalkken = tempken[0];
    varwalkken = tempken[1];
    double ** temppot=medvar(npassi, nwalk, tp1pot, workdir, "pot");
    medwalkpot = temppot[0];
    varwalkpot = temppot[1];
    double ** tempstim1=medvar(npassi, nwalk, tp1stim1, workdir, "stim1");
    medwalkstim1 = tempstim1[0];
    varwalkstim1 = tempstim1[1];
    double ** tempstim2=medvar(npassi, nwalk, tp1stim2, workdir, "stim2");
    medwalkstim2 = tempstim2[0];
    varwalkstim2 = tempstim2[1];

    cout << "---ken---" << endl;
    double * finken = medpesata(nwalk, medwalkken, varwalkken);
    cout << finken[0] << " " << sqrt(finken[1]) << endl;
    double * finken2 = medscema(nwalk, medwalkken, varwalkken);
    cout << finken2[0] << " " << sqrt(finken2[1]) << endl;
    cout << "---pot---" << endl;
    double * finpot = medpesata(nwalk, medwalkpot, varwalkpot);
    cout << finpot[0] << " " << sqrt(finpot[1]) << endl;
    double * finpot2 = medscema(nwalk, medwalkpot, varwalkpot);
    cout << finpot2[0] << " " << sqrt(finpot2[1]) << endl;
    cout << "---stim1---" << endl;
    double * finstim1 = medpesata(nwalk, medwalkstim1, varwalkstim1);
    cout << finstim1[0] << " " << sqrt(finstim1[1]) << endl;
    double * finstim12 = medscema(nwalk, medwalkstim1, varwalkstim1);
    cout << finstim12[0] << " " << sqrt(finstim12[1]) << endl;
    cout << "---stim2---" << endl;
    double * finstim2 = medpesata(nwalk, medwalkstim2, varwalkstim2);
    cout << finstim2[0] << " " << sqrt(finstim2[1]) << endl;
    double * finstim22 = medscema(nwalk, medwalkstim2, varwalkstim2);
    cout << finstim22[0] << " " << sqrt(finstim22[1]) << endl;
    
    return 0;
}
