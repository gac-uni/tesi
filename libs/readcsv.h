#include <cstdio>
#include <iostream>
#include <cmath>

class contenitorepazzo{
    public:
        std::string * datostr = NULL;
        double * datolf = NULL;
        bool isdatolf;
        int nrighe;
        int ncolonne;
        contenitorepazzo(std::string * puntatore, int nrow, int ncol);
        contenitorepazzo(double * puntatore, int nrow, int ncol);
        ~contenitorepazzo();
        bool ugualepezzo(int coordA, int dirA, contenitorepazzo * B, int coordB, int dirB);
        int getpezzo(int coord, int dir, double * dati);
        int getpezzo(int coord, int dir, std::string * dati);
        int transpose();
};
class readcsv{ 
    private:
        contenitorepazzo * headercontainer;
        contenitorepazzo * indexcontainer;
        contenitorepazzo * daticontainer;
        bool strheader=false;
        bool strindex=false;
        bool strdati=false;
        int nheader=0;
        int nindex=0;
        char sep=',';
        char newline='\n';
        //funzioni che servono per il constructor
        //trovo le dimensioni del file
        int getdim(int * nrighe, int * ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo l'header se è double
        double * readheaderlf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo l'header se è std::string
        std::string * readheaderstr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo l'indice se è double
        double * readindexlf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo l'indice se è std::string
        std::string * readindexstr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo i dati se sono double
        double * readdatilf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //leggo i dati se sono std::string
        std::string * readdatistr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept);
        //metto tutto insieme per il constructor così nell'overloading copio e incollo meno righe
        int preconstructor(FILE * filept);
    public:
        int nrighe=0;
        int ncolonne=0;                
        //constructor vari per fare cose
        readcsv(FILE * filept);
        readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom);
        readcsv(FILE * filept, int nindexcustom, int nheadercustom);
        readcsv(FILE * filept, int nindexcustom, int nheadercustom, char sepcustom, char newlinecustom);
        readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom, int nindexcustom, int nheadercustom);
        readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom, int nindexcustom, int nheadercustom, char sepcustom, char newlinecustom);
        //destructor
        ~readcsv();
        //funzioni per la visualizzazione a terminale
        int dispheader();
        int dispindex();
        int dispdati();
        //funzione per trasporre i dati
        int transpose();
        //funzioni per restituire in una variabile header, index e dati
        int header(std::string ** headerrichiesto);
        int header(double ** headerrichiesto);
        int index(std::string ** indexrichiesto);
        int index(double ** indexrichiesto);
        int dati(std::string ** datirichiesto);
        int dati(double ** datirichiesto);
        //lettura di dati usando index e header
        //prendo una riga a partire da un index
        template<class Ind, class Data>
        Data * getrow(Ind * indexcercato){
            Ind * temp = new Ind[nindex];
            for(int i=0; i<nindex; i++){
                temp[i]=indexcercato[i];
            }
            contenitorepazzo indextemporaneo(temp, 1, nindex);
            Data * out;
            out = new Data[ncolonne];
            bool trovato=false;
            for(int i=0; i<nrighe; i++){
                if(indextemporaneo.ugualepezzo(0, 0, indexcontainer, i, 0)){
                    (*daticontainer).getpezzo(i, 0, out);
                    trovato=true;
                    break;
                }
            }
            if(trovato){
                return out;
            }else{
                delete [] out;
                std::cerr<<"readcsv.getrow: nessuna corrispondenza trovata nei dati\n";
                return NULL;
            }
        }
        //prendo una colonna a partire da un header
        template<class Head, class Data>
        Data * getcol(Head * headcercato){
            Head * temp = new Head[nheader];
            for(int i=0; i<nheader; i++){
                temp[i]=headcercato[i];
            }
            contenitorepazzo headtemporaneo(temp, 1, nheader);
            Data * out;
            out = new Data[nrighe];
            bool trovato=false;
            for(int i=0; i<ncolonne; i++){
                if(headtemporaneo.ugualepezzo(0, 0, headercontainer, i, 1)){
                    (*daticontainer).getpezzo(i, 1, out);
                    trovato=true;
                    break;
                }
            }
            if(trovato){
                return out;
            }
            else{
                std::cerr<<"readcsv.getcol: nessuna corrispondenza trovata nei dati\n";
                delete [] out;
                return NULL;
            }
        }
        //prendo un dato a partire da index e header
        template<class Ind, class Head, class Data>
        Data getvalue(Ind * indexcercato, Head * headercercato){
            Data out;
            Data * rigacorretta = new Data[ncolonne];
            Ind * tempIndex = new Ind[nindex];
            for(int i=0; i<nindex; i++){
                tempIndex[i]=indexcercato[i];
            }
            contenitorepazzo indextemporaneo(tempIndex, 1, nindex);
            Head * tempHead = new Head[nheader];
            for(int i=0; i<nheader; i++){
                tempHead[i]=headercercato[i];
            }
            contenitorepazzo headtemporaneo(tempHead, 1, nheader);
            bool trovato=false;
            for(int i=0; i<nrighe; i++){
                if(indextemporaneo.ugualepezzo(0, 0, indexcontainer, i, 0)){
                    (*daticontainer).getpezzo(i, 0, rigacorretta);
                    trovato=true;
                    break;
                }
            }
            if(trovato){
                trovato=false;
                for(int i=0; i<ncolonne; i++){
                    if(headtemporaneo.ugualepezzo(0, 0, headercontainer, i, 1)){
                        out=rigacorretta[i];
                        trovato=true;
                        break;
                    }
                }
                delete [] rigacorretta;
                if(trovato){
                    return out;
                }else{
                    std::cerr<<"readcsv.getvalue: nessuna corrispondenza per l'index trovata nei dati\n";
                    return nanl("");
                }
            }else{
                std::cerr<<"readcsv.getvalue: nessuna corrispondenza per l'header trovata nei dati\n";
                delete [] rigacorretta;
                return nanl("");
            }
        }
};