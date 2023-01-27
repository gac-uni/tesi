/**
 * @file gpd.cpp
 * @author Giacomo Marcon (marcongiacomo@gmail.com)
 * @brief body of all the member function of the GnuplotDriver class
 * @version 0.0
 * @date 2023-01-20
 */
#include "gpd.h"
using std::cerr;
using std::string;
using std::endl;
using std::stringstream;
using std::ios;
 
int GnuplotDriver::sub(string * input, string toSub, string substitute){
	size_t toSubPos = input->find(toSub);
	while(toSubPos != string::npos){
		input->replace(toSubPos, substitute.length(), substitute);
		toSubPos=input->find(toSub);
	}
	return 0;
}

int GnuplotDriver::makepath(string path){
    string sep = "/";
    int depth=0;
    size_t possep=path.find(sep);
    size_t lastslash=-1;
    // finding last /
    while(possep!=string::npos){
        depth++;
        lastslash=possep;
        possep=path.find(sep, possep+1);
    }
    if(depth!=0){
        path=path.substr(0, lastslash); // deleting everything after last /
        string comando = "mkdir -p "+path; // creating required path 
        if(system(comando.c_str())==0){
            return 1;
        }else{
            return -1;
        }
        
    }
    return 0;
}

GnuplotDriver::GnuplotDriver() {
	ylabel = "Y [unità di Y]";
	xlabel = "X [unità di X]";
	zlabel = "Z [unità di Z]";
	title = "Titolo";
	raise = false;
	fitting=false;
	persist = false;
	matrix=false;
	logX=false;
	logY=false;
	save=false;
	overlayFunc=false;
	grid=true;
	toExec=false;
	leg=false;
	auto_color=true;
	printData=0; // 0 default (print solo se salvo il grafico in pdf) , 1 sempre print, -1 mai print
	trace_color = "blue";
	format="pdf";
	background_color = "white";
	limitx="";
	limity="";
	limitz="";
	limitx_fit="";
	limity_fit="";
	fileName_def="gpd_fig";
	fileName=fileName_def;
	printDataOutputName_default="gpd_data";
	printDataOutputName=printDataOutputName_default;
	format="pdf";
	lineStyle="linespoints";
	lineTitle = "";
	posLegend = "top right";
	gp = NULL;
	num_line=0;
	col_code = new string[4];
	col_code[0]="3F";
	col_code[1]="7F";
	col_code[2]="BF";
	col_code[3]="FF";
} 

GnuplotDriver::~GnuplotDriver() {
	int legal_stream=config_stream();
	if(gp){
		if(legal_stream){
			stringstream string_final;
			string_final << config_stream_string.str();
			if(!matrix && !fitting){
					if(colNum==1 || colNum==2){ //not doing a 3d plot
						string_final << "plot ";
					}else{ //doing a 3d plot
						string_final << "splot ";
					}
					string_final << buf.str() << endl;
					string_final << data.str();
			}else{
				string_final << buf.str() << data.str();
			}
			fprintf(gp, "%s",string_final.str().c_str());
			// printf("%s",string_final.str().c_str());
			if((printData == 1) || ((printData == 0) && save == true)){ 
				// saving data in file
				string printDataONext;
				if(printData == 1){ // file name was explicitly declared
					printDataONext = printDataOutputName+".dat";
				}else{ // file name is the same as the figure name
					printDataONext = fileName+".dat";
				}
				if(makepath(printDataONext)!=-1){
					outFile=fopen(printDataONext.c_str(), "w");
				}else{
					// on error generating path, fallback to standard name
					cerr<<"required path not legal, fallback to default name: " << printDataOutputName_default  << endl;
					string temp=printDataOutputName_default+".dat";
					outFile=fopen(temp.c_str(), "w");
				}
				string outData = data.str();
				sub(&outData, "e\n", "\n\n");
				fprintf(outFile, "%s", outData.c_str());
				fclose(outFile);
			}
		}else{
			cerr << "not possible to create a gnuplot request with selected options"<< endl;
		}
		fclose(gp);
	}
}  

// configurations that request an input
int GnuplotDriver::t(string arg){
	title=arg;
	return 0;
}
int GnuplotDriver::x(string arg){
	xlabel=arg;
	return 0;
}
int GnuplotDriver::y(string arg){
	ylabel=arg;
	return 0;
}
int GnuplotDriver::z(string arg){
	zlabel=arg;
	return 0;
}
int GnuplotDriver::c(string arg){
	trace_color=arg;
	auto_color=false;
	return 0;
}
int GnuplotDriver::fpath(string arg){
	save=true;
	fileName=arg;
	return 0;
}
int GnuplotDriver::fext(string arg){
	save=true;
	format=arg;
	return 0;
}
int GnuplotDriver::limx(double min, double max){
	limitx = "["+std::to_string(min)+":"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limx(string vuota, double max){
	limitx = "[:"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limx(double min, string vuota){
	limitx = "["+std::to_string(min)+":]";
	return 0;
}
int GnuplotDriver::limy(double min, double max){
	limity = "["+std::to_string(min)+":"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limy(string vuota, double max){
	limity = "[:"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limy(double min, string vuota){
	limity = "["+std::to_string(min)+":]";
	return 0;
}
int GnuplotDriver::limz(double min, double max){
	limitz = "["+std::to_string(min)+":"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limz(string vuota, double max){
	limitz = "[:"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limz(double min, string vuota){
	limitz = "["+std::to_string(min)+":]";
	return 0;
}
int GnuplotDriver::limxfit(double min, double max){
	limitx_fit = "["+std::to_string(min)+":"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limxfit(string vuota, double max){
	limitx_fit = "[:"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limxfit(double min, string vuota){
	limitx_fit = "["+std::to_string(min)+":]";
	return 0;
}
int GnuplotDriver::limyfit(double min, double max){
	limity_fit = "["+std::to_string(min)+":"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limyfit(string vuota, double max){
	limity_fit = "[:"+std::to_string(max)+"]";
	return 0;
}
int GnuplotDriver::limyfit(double min, string vuota){
	limity_fit = "["+std::to_string(min)+":]";
	return 0;
}
int GnuplotDriver::func(string arg){
	overlayFunc = true;
	funcToOverlay=arg;
	return 0;
}
int GnuplotDriver::ls(string arg){
	lineStyle=arg;
	return 0;
}
int GnuplotDriver::lt(string arg){
	leg=true;
	lineTitle=arg;
	return 0;
}
int GnuplotDriver::keyloc(string arg){
	posLegend=arg;
	return 0;
}
int GnuplotDriver::print(string arg){
	printData=1;
	printDataOutputName=arg;
	return 0;
}

// configurations that require no input
int GnuplotDriver::h(){
	cerr << "Usage:	[c (trace color, default \"" << trace_color << "\")]" << endl; 		
	cerr << "		[fExt (Extension of the figure file, default \"" <<format<< "\")]" << endl;
	cerr << "		[fPath (path for the figure file, default \""<<fileName_def<<"\")]"<<endl;
	cerr << "		[func (overlay function for the figure)]" << endl;
	cerr << "		[h (opens help prompt)]" << endl;
	cerr << "		[keyLoc (location of the key in the figure, default \""<<posLegend<<"\")]" <<endl;
	cerr << "		[limX (limits for x axis with notation [xmin:xmax])]" << endl;
	cerr << "		[limY (limits for y axis with notation [ymin:ymax])]" << endl;
	cerr << "		[limZ (limits for z axis with notation [zmin:zmax])]" << endl;
	cerr << "		[limXFit (limits for the points considered in the fitting, same notation as limX)]" << endl;
	cerr << "		[limYFit (limits for the points considered in the fitting, same notation as limY)]" << endl;
	cerr << "		[logX (sets logarithmic scale for x axis)]" << endl;
	cerr << "		[logY (sets logarithmic scale for y axis)]" << endl;
	cerr << "		[ls (style of the line, one of: \"points\", \"line\", \"linepoints\", default \""<<lineStyle<<"\")]" << endl;
	cerr << "		[lt (title of the line, appearing in the key)]" << endl;
	cerr << "		[m (interprets data as a heatmap)]" << endl;
	cerr << "		[noGrid (turns off grid)]" << endl;
	cerr << "		[noPrint (forces to not save data to file even if the figure is saved to file)]" << endl;
	cerr << "		[print (forces to save data to file even if the figure is not saved to file)]" << endl;
	cerr << "       [p (persist)]" <<endl;
	cerr << "       [t (chart title, default \"" << title << "\"]" << endl;
	cerr << "       [x (x axis title, default \"" << xlabel << "\")]" << endl;
	cerr << "       [y (y axis title, default \"" << ylabel << "\")]" << endl;
	cerr << "       [z (z axis title, default \"" << zlabel << "\")]" << endl;
	return 0;
}
int GnuplotDriver::p(){
	persist=true;
	return 0;
}
int GnuplotDriver::m(){
	matrix=true;
	return 0;
}
int GnuplotDriver::logx(){
	logX=true;
	return 0;
}
int GnuplotDriver::logy(){
	logY=true;
	return 0;
}
int GnuplotDriver::nogrid(){
	grid=false;
	return 0;
}
int GnuplotDriver::noprint(){
	printData=-1;
	return 0;
}

bool GnuplotDriver::config_stream(){
	if(!toExec && !fitting){
		return false;
	}
	//opening gnuplot stream with requested options
	if(!gp) {
		string com = (string) "gnuplot -background " + background_color + "";
		gp = popen(com.c_str(),"w");
	}
	// verify that the requested options are compatible
	//logarithmic scale makes no sense in heatmaps
	if((matrix==true)&&((logX==true)||(logY==true))){
		cerr<<"-m and -z/-w cannot coexist"<<endl;
		return false;
	}
	//pdfcairo will not accept -persist option
	if(persist==true && save==true){
		cerr<<"-f e -p cannot coexist"<<endl;
		return false;
	} 
	//-o e -m cannot coexist 
	if(matrix==true && overlayFunc==true){
		cerr<<"-m and -p cannot coexist, you cannot overlay a function to a heatmap"<< endl;
		return false;
	}
	//checking that requested linestyle is supported
	if(lineStyle != "line" && lineStyle != "points" && lineStyle != "linespoints"){
		cerr <<"the specified line style is not supported" << endl;
		return false;
	}
	// 2d plots do not accept limits on the z axis
	if(colNum<3 && (limitz != "" || matrix)){
		cerr << "2d plots do not accept limits on the z axis" << endl;
		return false;
	}
	if(fitting == false && (limitx_fit!="" || limity_fit!="")){
		cerr << "limits on the fitting data are not accepted if no fit is requested" << endl;
	}

	config_stream_string << "set title '" << title << "'" << endl;
	config_stream_string << "set xlabel '" << xlabel << "'" << endl;
	config_stream_string << "set ylabel '" << ylabel << "'" << endl;
	config_stream_string << "set zlabel '" << zlabel << "'" << endl;
	config_stream_string << "unset cblabel" << endl;
	config_stream_string << "unset key" << endl;
	//setting gnuplots options as requested
	if(leg){
		config_stream_string << "set key " << posLegend << endl;
	}
	if(grid){
		config_stream_string << "set grid " << endl;
	}	
	if(logX){
		config_stream_string<<"set logscale x"<<endl;
	}
	if(logY){
		config_stream_string<<"set logscale y"<<endl;
	}
	if(fitting){
		if(limitx_fit!=""){
			config_stream_string << "set xrange " << limitx_fit << endl;
		}
		if(limity_fit!=""){
			config_stream_string << "set yrange " << limity_fit << endl;
		}
	}else{
		if(limitx!=""){
			config_stream_string << "set xrange " << limitx << endl;
		}
		if(limity!=""){
			config_stream_string << "set yrange " << limity << endl;
		}
		if(limitz!=""){
			config_stream_string << "set zrange " << limitz << endl;
		}
	}
	if(persist){
		config_stream_string << "set terminal qt persist" << endl;
	}
	if(save){
		if(format=="pdf"){
			config_stream_string<<"set terminal pdfcairo"<<endl;
			if(makepath(fileName)==-1){
				cerr<<"requested path is not legal, fallback to default name: " << fileName_def << endl;
				fileName=fileName_def;
			}
			config_stream_string<<"set out \""<<fileName<<".pdf\""<<endl;
		}else if(format=="jpg"){
			config_stream_string<<"set terminal jpeg"<<endl;
			if(makepath(fileName)==-1){
				cerr<<"requested path is not legal, fallback to default name: " << fileName_def <<endl;
				fileName=fileName_def;
			}
			config_stream_string<<"set out \""<<fileName<<".jpg\""<<endl;
		}else{
			cerr<<"requested format is not supported" << endl;
			return false;
		}
	}
	//declaring the function to overlay if requested
	if(overlayFunc){
		config_stream_string << funcToOverlay << endl;
	}
	return true;
}

// creates the line that plots the data
int GnuplotDriver::plotcommand(){
	toExec=true;
	try{
		if(matrix==false){ // not plotting a heatmap
			buf << " \"-\"";
			switch(colNum){
				case 1:
					buf << " u 0:1 w ";
				break;
				case 2:
					buf << " u 1:2 w ";
				break;
				case 3:
					buf << " u 1:2:3 w ";
				break;
			}
			buf << lineStyle;
			if(!leg || lineTitle == ""){
				buf << " notitle";
			}else{
				buf << " t \"" << lineTitle << "\"";
				lineTitle = ""; // if not requested the title of the following line must have the default title
			}
			if(auto_color){
				int b=num_line%4;
				int g=(num_line/4)%4;
				int r=(num_line/16)%4;
				std::string coloretemp= "#"+col_code[r]+col_code[g]+col_code[b];
				buf << " lc rgb \"" << coloretemp << "\"";
				num_line < 63 ? num_line++ : num_line = 0;	
			}else{
				buf << " lc rgb \"" << trace_color << "\"";
				auto_color=true;
			}
			if(overlayFunc){
				int end=funcToOverlay.find("=");
				buf << ", " << funcToOverlay.substr(0, end);
			}
			buf << ", ";
			data << "e" << endl;
		}else{ // plotting a heatmap
			buf << "plot " <<" \"-\" matrix with image" << endl;
			data << "e" << endl;
		}
	}
	catch (ios::failure const &problem) { 
		cerr << "gnuplot_driver: " << problem.what() << endl;
	}
	return 0;
}

int GnuplotDriver::plot(double * dati, int nRows, int nCols){
	bool legalData=true;	
	colNum=nCols;
	rowNum=nRows;	
	if(matrix==false && colNum>3){
		cerr<<"cannot plot a graph with dimension greater than 3"<<endl;
		legalData=false;
	}
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			for(int j=0;j<colNum;j++){
				data<<dati[i*colNum+j] << " ";
			}
			data<<endl;
		}	
		plotcommand();	
		return 0;
	}else{
		cerr << "cannot create a stream with the requested options"<< endl;
		return -1;
	}
}
int GnuplotDriver::plot(double * dati1, double * dati2, int nRows){
	bool legalData=true;
	rowNum=nRows;
	colNum=2;
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			data<< dati1[i] << " " << dati2[i] << endl;
		}	
		plotcommand();	
		return 0;
	}else{
		cerr << "cannot create a stream with the requested options"<< endl;
		return -1;
	}
}
int GnuplotDriver::plot(double * dati1, double * dati2, double * dati3, int nRows){
	bool legalData=true;
	rowNum=nRows;
	colNum=3;	
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			data<< dati1[i] << " " << dati2[i] << " " << dati3[i] << endl;
		}	
		plotcommand();	
		return 0;
	}else{
		cerr << "cannot create a stream with the requested options"<< endl;
		return -1;
	}
}

int GnuplotDriver::fitcommand(string inputFitFunc, string params){
	int nVar=0;
	int posEquals=inputFitFunc.find("=");
	string funcName=inputFitFunc.substr(0, posEquals);
	int openPar=funcName.find("(");
	int closedPar=funcName.find(")");
	if(closedPar-openPar==2){
		nVar=1;
	}else if(closedPar-openPar==4){
		nVar=2;
	}
	string fitFunc=inputFitFunc.substr(openPar, inputFitFunc.length());
	string varFunc=inputFitFunc.substr(openPar, closedPar);
	try{
		buf << "fitFunc" << fitFunc << endl;
		buf << "set fit logfile \"" << fileName << "_fitlog.txt\" " << endl;
		buf << "set fit quiet" << endl;
		buf << "fit " << " fitFunc" << varFunc << " \"-\" ";
		if(nVar=1 && colNum==1){
			buf << "u 0:1 "; 
		}else if(nVar=1 && colNum==2){
			buf << "u 1:2 ";
		}else if(nVar=2 && colNum==2){
			buf << "u 0:1:2 ";
		}else if(nVar=2 && colNum==3){
			buf << "u 1:2:3 ";
		}
		buf << "via " << params << endl;
		buf << data.str();
		buf << "e" << endl;

		if(nVar==1){
			buf << "unset xrange" << endl << "unset yrange" << endl;
			if(limitx!=""){
				buf << "set xrange " << limitx << endl;
			}
			if(limity!=""){
				buf << "set yrange " << limity << endl;
			}
			buf << "set key " << posLegend << endl;
			buf << "plot " << " \"-\" ";
			if(colNum==1){
				buf << "u 0:1 ";
			}else if(colNum==2){
				buf << "u 1:2 ";
			}
			buf << "w " << lineStyle << " t \"data\", fitFunc" << varFunc << " t \"fit\" "<<endl;
			data << "e" << endl;
		}
	}
	catch (ios::failure const &problem) {
		cerr << "gnuplot_driver: " << problem.what() << endl;
	}
	return 0;
}
int GnuplotDriver::fit(double * dati, int nRows, int nCols, string fitFunc, string params){
	fitting=true;
	bool legalData=true;
	rowNum=nRows;
	colNum=nCols;
	if(colNum>3){
		cerr<<"cannot create plots with dimension greater than 3"<<endl;
		legalData=false;
	}
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			for(int j=0;j<colNum;j++){
				data<<dati[i*colNum+j] << " ";
			}
			data<<endl;
		}	
		fitcommand(fitFunc, params);
		return 0;
	}else{
		cerr << "cannot create a valid stream with requested options"<< endl;
		return -1;
	}
}
int GnuplotDriver::fit(double * dati1, double * dati2, int nRows, string fitFunc, string params){
	fitting=true;
	bool legalData=true;
	rowNum=nRows;
	colNum=2;
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			data<< dati1[i] << " " << dati2[i] << endl;
		}		
		fitcommand(fitFunc, params);
		return 0;
	}else{
		cerr << "cannot create a valid stream with requested options"<< endl;
		return -1;
	}
}
int GnuplotDriver::fit(double * dati1, double * dati2, double * dati3, int numRighe, string fitFunc, string params){
	fitting=true;
	bool legalData=true;
	rowNum=numRighe;
	colNum=3;
	if(legalData){	
		for(int i=0;i<rowNum;i++){
			data<< dati1[i] << " " << dati2[i] << " " << dati3[i] << endl;
		}	
		fitcommand(fitFunc, params);
		return 0;
	}else{
		cerr << "cannot create a valid stream with requested options"<< endl;
		return -1;
	}
}