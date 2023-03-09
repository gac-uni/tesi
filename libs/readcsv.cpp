#include "readcsv.h"

contenitorepazzo::contenitorepazzo(std::string * puntatore, int nrow, int ncol){
    datostr=puntatore;
    isdatolf=false;
    nrighe=nrow;
    ncolonne=ncol;
}
contenitorepazzo::contenitorepazzo(double * puntatore, int nrow, int ncol){
    datolf=puntatore;
    isdatolf=true;
    nrighe=nrow;
    ncolonne=ncol;
}
contenitorepazzo::~contenitorepazzo(){
    isdatolf ? delete [] datolf : delete [] datostr;
}
bool contenitorepazzo::ugualepezzo(int coordA, int dirA, contenitorepazzo * B, int coordB, int dirB){                
    if(isdatolf != (*B).isdatolf){
        std::cerr<<"contenitorepazzo.ugualepezzo: sto confrontando oggetti con tipi diversi\n";
        return false;
    }
    if((dirA == 0) && (dirB == 0)){
        if(ncolonne != (*B).ncolonne){
            return false;
        }
        bool uguale = true;
        if(isdatolf){
            for(int i=0; i<ncolonne; i++){
                if(datolf[ncolonne*coordA+i]!=(*B).datolf[ncolonne*coordB+i]){
                    uguale=false;
                    break;
                }
            }
        }else{
            for(int i=0; i<ncolonne; i++){
                if(datostr[ncolonne*coordA+i]!=(*B).datostr[ncolonne*coordB+i]){
                    uguale=false;
                    break;
                }
            }
        }
        if(uguale){
            return true;
        }else{
            return false;
        }
    }
    if((dirA == 1) && (dirB == 1)){
        if(nrighe != (*B).nrighe){
            return false;
        }
        bool uguale = true;
        if(isdatolf){
            for(int i=0; i<nrighe; i++){
                if(datolf[ncolonne*i+coordA]!=(*B).datolf[(*B).ncolonne*i+coordB]){
                    uguale=false;
                    break;
                }
            }
        }else{
            for(int i=0; i<ncolonne; i++){
                if(datostr[ncolonne*i+coordA]!=(*B).datostr[(*B).ncolonne*i+coordB]){
                    uguale=false;
                    break;
                }
            }
        }
        if(uguale){
            return true;
        }else{
            return false;
        }
    }
    if((dirA == 0) && (dirB == 1)){
        if(ncolonne != (*B).nrighe){
            return false;
        }
        bool uguale = true;
        if(isdatolf){
            for(int i=0; i<ncolonne; i++){
                if(datolf[ncolonne*coordA+i]!=(*B).datolf[(*B).ncolonne*i+coordB]){
                    uguale=false;
                    break;
                }
            }
        }else{
            for(int i=0; i<ncolonne; i++){
                if(datostr[ncolonne*coordA+i]!=(*B).datostr[(*B).ncolonne*i+coordB]){
                    uguale=false;
                    break;
                }
            }
        }
        if(uguale){
            return true;
        }else{
            return false;
        }                    
    }
    if((dirA == 1) && (dirB == 0)){
        if(nrighe != (*B).ncolonne){
            return false;
        }
        bool uguale = true;
        if(isdatolf){
            for(int i=0; i<nrighe; i++){
                if(datolf[ncolonne*i+coordA]!=(*B).datolf[(*B).ncolonne*coordB+i]){
                    uguale=false;
                    break;
                }
            }
        }else{
            for(int i=0; i<ncolonne; i++){
                if(datostr[ncolonne*i+coordA]!=(*B).datostr[(*B).ncolonne*coordB+i]){
                    uguale=false;
                    break;
                }
            }
        }
        if(uguale){
            return true;
        }else{
            return false;
        }
    }
    return false;
}
int contenitorepazzo::getpezzo(int coord, int dir, double * dati){
    if(isdatolf==false){
        std::cerr<<"contenitorepazzo.getpezzo: il tipo richiesto è double ma il tipo disponibile è std::string\n";
        return -1;
    }
    if(dir==0){
        for(int i=0; i<ncolonne; i++){
            dati[i]=datolf[ncolonne*coord+i];
        }
    }else{
        for(int i=0; i<nrighe; i++){
            dati[i]=datolf[ncolonne*i+coord];
        }
    }
    return 0;
}
int contenitorepazzo::getpezzo(int coord, int dir, std::string * dati){
    if(isdatolf==true){
        std::cerr<<"contenitorepazzo.getpezzo: il tipo richiesto è std::string ma il tipo disponibile è double\n";
        return -1;
    }
    if(dir==0){
        for(int i=0; i<ncolonne; i++){
            dati[i]=datostr[ncolonne*coord+i];
        }
    }else{
        for(int i=0; i<nrighe; i++){
            dati[i]=datostr[ncolonne*i+coord];
        }
    }
    return 0;
}
int contenitorepazzo::transpose(){
    if(isdatolf){
        double * tempdato = new double[nrighe*ncolonne];
        for(int i=0; i<nrighe; i++){
            for(int j=0; j<ncolonne; j++){
                tempdato[nrighe*j+i] = datolf[ncolonne*i+j];
            }
        }
        delete [] datolf;
        datolf=tempdato;
    }else{
        std::string * tempdato = new std::string[nrighe*ncolonne];
        for(int i=0; i<nrighe; i++){
            for(int j=0; j<ncolonne; j++){
                tempdato[nrighe*j+i] = datolf[ncolonne*i+j];
            }
        }
        delete [] datostr;
        datostr=tempdato;
    }
    int temp=nrighe;
    nrighe=ncolonne;
    ncolonne=temp;
    return 0;
}
int readcsv::getdim(int * nrighe, int * ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    int ncolonnetemp=0;
    char ch;
    int newlbefore=0;
    while(EOF!=(ch=getc(filept))){
        if(ch==sep || ch==newline){
            (*nrighe) == 0 ? (*ncolonne)++ : ncolonnetemp++;
        }
        if(ch==newline){
            if((*nrighe)!=0 && ncolonnetemp!=(*ncolonne)){
                std::cerr<<"ci sono righe di lunghezza diversa nei dati forniti."<<std::endl;
                return -1;
            }
            (*nrighe)++;
            ncolonnetemp=0;
        }
        if(ch==newline){
            newlbefore++;
        }else{
            newlbefore=0;
        }
    }
    (*nrighe)-=nheader;
    (*nrighe)++;
    (*nrighe)-=newlbefore;
    (*ncolonne)-=nindex;
    fseek(filept, 0, SEEK_SET);
    return 0;
}
double * readcsv::readheaderlf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    double * header;
    header=new double[nheader*ncolonne];
    for(int i=0; i<nheader;i++){
        //tiro via dall'inizio le prime caselle vuote
        ch=getc(filept);
        int vuoterimanenti=nindex;
        if(ch==sep){
            vuoterimanenti-=1;
        } 
        while(vuoterimanenti>0){
            ch=getc(filept);
            if(ch==sep){
                vuoterimanenti-=1;
            } 
        }
        //prendo i valori dell'header
        for(int j=0;j<ncolonne;j++){
            fscanf(filept, "%lf,", &header[i*ncolonne+j]);
        }
    }
    fseek(filept, 0, SEEK_SET);
    return header;
}
std::string * readcsv::readheaderstr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    char nome[1000];
    std::string * header;
    header=new std::string[nheader*ncolonne];
    for(int i=0; i<nheader;i++){
        //tiro via dall'inizio le prime caselle vuote
        ch=getc(filept);
        int vuoterimanenti=nindex;
        if(ch==sep){
            vuoterimanenti-=1;
        } 
        while(vuoterimanenti>0){
            ch=getc(filept);
            if(ch==sep){
                vuoterimanenti-=1;
            } 
        }
        //prendo i valori dell'header
        int ncolonna=0;
        while(ncolonna<ncolonne){
            int lunghezzanome=0;
            ch=getc(filept);
            while(ch!=sep && ch!=newline){
                nome[lunghezzanome]=ch;
                lunghezzanome++;
                ch=getc(filept);
            }
            header[i*ncolonne+ncolonna].assign(nome, lunghezzanome);
            ncolonna++;
        }

    }
    fseek(filept, 0, SEEK_SET);
    return header;
}
double * readcsv::readindexlf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    double * index;
    index=new double[nrighe*nindex];
    //tiro via le righe iniziali
    for(int i=0; i<nheader; i++){
        ch=getc(filept);
        while(ch!=newline){ 
            ch=getc(filept);
        }
    }
    //prendo i valori dell'index
    for(int i=0; i<nrighe; i++){
        //prendo i valori dell'indice all'inizio della riga
        for(int j=0;j<nindex;j++){
            fscanf(filept, "%lf,", &index[i*nindex+j]);
        }
        //finisco di leggere la riga
        ch=getc(filept);
        while(ch!=newline && ch!=EOF){
            ch=getc(filept);
        }
    }
    fseek(filept, 0, SEEK_SET);
    return index;
}
std::string * readcsv::readindexstr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    char nome[1000];
    std::string * index;
    index=new std::string[nrighe*ncolonne];
    //tiro via le righe iniziali
    for(int i=0; i<nheader; i++){
        ch=getc(filept);
        while(ch!=newline){ 
            ch=getc(filept);
        }
    }
    //prendo i valori dell'index
    for(int i=0; i<nrighe; i++){
        //prendo i valori dell'indice all'inizio della riga
        int nindice=0;
        while(nindice<nindex){
            int lunghezzanome=0;
            ch=getc(filept);
            while(ch!=sep && ch!=newline){
                nome[lunghezzanome]=ch;
                lunghezzanome++;
                ch=getc(filept);
            }
            index[i*nindex+nindice].assign(nome, lunghezzanome);
            nindice++;
        }
        //finisco di leggere la riga
        ch=getc(filept);
        while(ch!=newline && ch!=EOF){
            ch=getc(filept);
        }
    }
    fseek(filept, 0, SEEK_SET);
    return index;
}
double * readcsv::readdatilf(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    double * dati;
    dati=new double[nrighe*ncolonne];
    //tiro via le righe iniziali
    for(int i=0; i<nheader; i++){
        ch=getc(filept);
        while(ch!=newline){ 
            ch=getc(filept);
        }
    }
    //leggo i dati
    for(int i=0; i<nrighe; i++){
        //tiro via dall'inizio l'indice
        ch=getc(filept);
        int vuoterimanenti=nindex;
        if(ch==sep){
            vuoterimanenti-=1;
        } 
        while(vuoterimanenti>0){
            ch=getc(filept);
            if(ch==sep){
                vuoterimanenti-=1;
            } 
        }
        //prendo i dati
        for(int j=0;j<ncolonne;j++){
            fscanf(filept, "%lf,", &dati[i*ncolonne+j]);
        }
    }
    fseek(filept, 0, SEEK_SET);
    return dati;
}
std::string * readcsv::readdatistr(int nrighe, int ncolonne, int nheader, int nindex, char sep, char newline, FILE * filept){
    char ch;
    char nome[1000];
    std::string * dati;
    dati=new std::string[nrighe*ncolonne];
    //tiro via le righe iniziali
    for(int i=0; i<nheader; i++){
        ch=getc(filept);
        while(ch!=newline){ 
            ch=getc(filept);
        }
    }
    //leggo i dati
    for(int i=0; i<nrighe; i++){
        //tiro via dall'inizio l'indice
        ch=getc(filept);
        int vuoterimanenti=nindex;
        if(ch==sep){
            vuoterimanenti-=1;
        } 
        while(vuoterimanenti>0){
            ch=getc(filept);
            if(ch==sep){
                vuoterimanenti-=1;
            } 
        }
        //prendo i dati
        int ncolonna=0;
        while(ncolonna<ncolonne){
            int lunghezzanome=0;
            ch=getc(filept);
            while(ch!=sep && ch!=newline && ch!=EOF){
                nome[lunghezzanome]=ch;
                lunghezzanome++;
                ch=getc(filept);
            }
            dati[i*ncolonne+ncolonna].assign(nome, lunghezzanome);
            ncolonna++;
        }
    }
    fseek(filept, 0, SEEK_SET);
    return dati;
}
int readcsv::preconstructor(FILE * filept){
    getdim(&nrighe, &ncolonne, nheader, nindex, sep, newline, filept);
    if(strheader){
        std::string * temp = readheaderstr(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        headercontainer=new contenitorepazzo(temp, nheader, ncolonne);
    }else{
        double * temp = readheaderlf(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        headercontainer=new contenitorepazzo(temp, nheader, ncolonne);
    }
    if(strindex){
        std::string * temp = readindexstr(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        indexcontainer=new contenitorepazzo(temp, nrighe, nindex);
    }else{
        double * temp= readindexlf(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        indexcontainer=new contenitorepazzo(temp, nrighe, nindex);
    }
    if(strdati){
        std::string * temp = readdatistr(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        daticontainer=new contenitorepazzo(temp, nrighe, ncolonne);
    }else{
        double * temp = readdatilf(nrighe, ncolonne, nheader, nindex, sep, newline, filept);
        daticontainer=new contenitorepazzo(temp, nrighe, ncolonne);
    }
    return 0;
}
readcsv::readcsv(FILE * filept){
    preconstructor(filept);     
}
readcsv::readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom){
    strheader=strheadercustom;
    strindex=strindexcustom;
    strdati=strdaticustom;
    preconstructor(filept);           
}
readcsv::readcsv(FILE * filept, int nindexcustom, int nheadercustom){
    nindex=nindexcustom;
    nheader=nheadercustom;
    preconstructor(filept);            
}
readcsv::readcsv(FILE * filept, int nindexcustom, int nheadercustom, char sepcustom, char newlinecustom){
    nindex=nindexcustom;
    nheader=nheadercustom;
    sep=sepcustom;
    newline=newlinecustom;
    preconstructor(filept);            
}
readcsv::readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom, int nindexcustom, int nheadercustom){
    strheader=strheadercustom;
    strindex=strindexcustom;
    strdati=strdaticustom;
    nindex=nindexcustom;
    nheader=nheadercustom;
    preconstructor(filept);            
}
readcsv::readcsv(FILE * filept, bool strheadercustom, bool strindexcustom, bool strdaticustom, int nindexcustom, int nheadercustom, char sepcustom, char newlinecustom){
    strheader=strheadercustom;
    strindex=strindexcustom;
    strdati=strdaticustom;
    nindex=nindexcustom;
    nheader=nheadercustom;
    sep=sepcustom;
    newline=newlinecustom;
    preconstructor(filept);            
}
readcsv::~readcsv(){
    delete headercontainer;
    delete indexcontainer;
    delete daticontainer;
}
int readcsv::dispheader(){
    if(strheader){
        std::cout<<"---header---"<<std::endl;
        for(int i=0; i<nheader; i++){
            for(int j=0; j<ncolonne; j++){
                std::cout<<(*headercontainer).datostr[ncolonne*i+j]<<"  ";
            }
            std::cout<<std::endl;
        }
    }else{
        std::cout<<"---header---"<<std::endl;
        for(int i=0; i<nheader; i++){
            for(int j=0; j<ncolonne; j++){
                std::cout<<(*headercontainer).datolf[ncolonne*i+j]<<"  ";
            }
            std::cout<<std::endl;
        }
    }
    return 0;
}
int readcsv::dispindex(){
    if(strindex){
        std::cout<<"---index---"<<std::endl;
        for(int i=0; i<nrighe; i++){
            for(int j=0; j<nindex; j++){
                std::cout<<(*indexcontainer).datostr[nindex*i+j]<<"  ";
            }
            std::cout<<std::endl;
        }
    }else{
        std::cout<<"---index---"<<std::endl;
        for(int i=0; i<nrighe; i++){
            for(int j=0; j<nindex; j++){
                std::cout<<(*indexcontainer).datolf[nindex*i+j]<<"  ";
            }
            std::cout<<std::endl;
        }
    }
    return 0;
}
int readcsv::dispdati(){
    if(strdati){
        std::cout<<"---dati---"<<std::endl;
        for(int i=0;i<nrighe; i++){
            for(int j=0;j<ncolonne;j++){
                std::cout<<(*daticontainer).datostr[ncolonne*i+j]<< " ";
            }
            std::cout<<std::endl;
        }
    }else{
        std::cout<<"---dati---"<<std::endl;
        for(int i=0;i<nrighe; i++){
            for(int j=0;j<ncolonne;j++){
                std::cout<<(*daticontainer).datolf[ncolonne*i+j]<< " ";
            }
            std::cout<<std::endl;
        }
    }
    return 0;
}
int readcsv::transpose(){
    (*headercontainer).transpose();
    (*indexcontainer).transpose();
    (*daticontainer).transpose();
    contenitorepazzo *tempcontenitore;
    tempcontenitore=headercontainer;
    headercontainer=indexcontainer;
    indexcontainer=tempcontenitore;
    //sistemo le variabili di dimensione
    int temp=nindex;
    nindex=nheader;
    nheader=temp;
    temp=ncolonne;
    ncolonne=nrighe;
    nrighe=temp;
    return 0;
}
int readcsv::header(std::string ** headerrichiesto){
    *headerrichiesto=(*headercontainer).datostr;
    return 0;
}
int readcsv::header(double ** headerrichiesto){
    *headerrichiesto=(*headercontainer).datolf;
    return 0;
}
int readcsv::index(std::string ** indexrichiesto){
    *indexrichiesto=(*indexcontainer).datostr;
    return 0;
}
int readcsv::index(double ** indexrichiesto){
    *indexrichiesto=(*indexcontainer).datolf;
    return 0;
}
int readcsv::dati(std::string ** datirichiesto){
    *datirichiesto=(*daticontainer).datostr;
    return 0;
}
int readcsv::dati(double ** datirichiesto){
    *datirichiesto=(*daticontainer).datolf;
    return 0;
}