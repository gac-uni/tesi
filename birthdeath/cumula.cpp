#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

double * mediacumulata(double * vettore, int nelementi, int elementoiniziale){
    double * cumulata;
    cumulata=new double[nelementi-elementoiniziale];
    cumulata[0]=vettore[elementoiniziale];
    for(int i=1;i<nelementi-elementoiniziale;i++){
        cumulata[i]+=cumulata[i-1]+vettore[elementoiniziale+i];
        cumulata[i-1]/=i;
    }
    cumulata[nelementi-elementoiniziale-1]/=nelementi-elementoiniziale;
    return cumulata;
}

int main(int argc, char ** argv){
    int piniziale = 0;
    string nome = argv[1];
    if(argc == 3){
        sscanf(argv[2], "%d", &piniziale);
    }
    cout << piniziale << endl;
    ifstream fin(nome+".dat");
    string temp;
    vector<double> v;
    if(fin.is_open()){
        while(getline(fin, temp)){
            v.push_back(stod(temp));
        }
    }
    fin.close();
    double * vcum = mediacumulata(v.data(), v.size(), piniziale);
    ofstream fout(nome+"cumman.dat");
    if(fout.is_open()){
        for(int i=0;i<v.size()-piniziale;i++){
            fout << vcum[i] << endl;
        }
    }
    fout.close();


    return 0;
}