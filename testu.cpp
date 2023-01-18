#include <iostream>
#include <cmath>
#include "gpd.h"


using namespace std;

double dzeroU(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    */
    return pow(parametri[0]/r, parametri[1]);
}
double dprimaU(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    */
    return -pow(parametri[0], parametri[1])*parametri[1]*pow(r, -parametri[1]-1);
}
double dsecondaU(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    */
    return -pow(parametri[0], parametri[1])*parametri[1]*(parametri[1]+1)*pow(r, -parametri[1]-2);
}
double UT(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    2) lscatola
    */
    if(r<=parametri[2]/2.){
        return dzeroU(r, parametri)+dzeroU(parametri[2]-r, parametri)-2*dzeroU(parametri[2]/2., parametri);
    }else{
        return 0;
    }
}
double dprimaUT(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    2) lscatola
    */
    if(r<=parametri[2]/2.){
        return dprimaU(r, parametri)-dprimaU(parametri[2]-r, parametri);
    }else{
        return 0;
    }
}
double dsecondaUT(double r, double * parametri){
    /*
    parametri: 
    0) a1
    1) a2
    2) lscatola
    */
    if(r<=parametri[2]/2.){
        return dsecondaU(r, parametri)+dsecondaU(parametri[2]-r, parametri);
    }else{
        return 0;
    }
}
GnuplotDriver gp;
int main(){
    const int npunti = 1000;
    double parametriU[3] = {2.5, 5, 20};
    double eallau[npunti]={0};
    double dist[npunti]={0};
    double estremi[2]={0,10};
    double passo=(estremi[1]-estremi[0])/npunti;
    for(int i=0;i<npunti;i++){
        dist[i]=estremi[0]+passo*i;
        eallau[i]=exp(-UT(dist[i], parametriU));
    }
    gp.conf("ls", "line");
    gp.conf("p");
    gp.plot(dist, eallau, npunti);
    return 0;
}
