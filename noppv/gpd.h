/**
 * @file gpd.h
 * @author Giacomo Marcon (marcongiacomo@gmail.com)
 * @brief Header for the gpd.cpp class
 * @version 0.0
 * @date 2023-01-20
 */
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <cstdlib>

/**
 * @brief Main class, used to plot data, fit data and change preferences.
 * 
 * @latexonly 
 * \input{../tex_sources/classe}
 * @endlatexonly
 * 
 * 
 * 
 */
class GnuplotDriver{
    std::string xlabel,ylabel,zlabel, title, trace_color,background_color, fileName, format, limitx, limity, limitz, limitx_fit, limity_fit, funcToOverlay, lineStyle, lineTitle, posLegend, printDataOutputName, fileName_def, printDataOutputName_default;
    std::string * col_code;
    bool inverse,raise,persist,logX,logY,matrix, save, overlayFunc,fitting, grid, toExec, leg, auto_color;
    std::stringstream config_stream_string, buf, data;
    int colNum, rowNum, num_line, printData;
    FILE *gp, *outFile;
    public:

    /**
     * @brief Construct a new Gnuplot Driver:: Gnuplot Driver object, setting default values.
     * 
     */
    GnuplotDriver(); 

    /**
     * @brief Destroy the Gnuplot Driver:: Gnuplot Driver object, asseble (merges command stream with data stream) and send command stream to gnuplot and saves data to file if required.
     * 
     */
    ~GnuplotDriver();

    /**
     * @brief Sets title of the figure.
     * 
     * @param title title to be given to the figure.
     */
    int t(std::string title);

    /**
     * @brief Sets label of the x axis.
     * 
     * @param x_title title to be given to the x axis.
     */
    int x(std::string x_title);

    /**
     * @brief Sets label of the y axis.
     * 
     * @param y_title title to be given to the y axis.     
     */
    int y(std::string y_title);

    /**
     * @brief Sets label of the z axis.
     * 
     * @param z_title title to be given to the z axis.     
     */
    int z(std::string z_title);

    /**
     * @brief Sets color of the line.
     * 
     * The color must be specified in exadecimal format RRGGBB. If Multiple lines are plotted into the same figure
     * then the color applies only to the first line coming after the GnuplotDriver.c() call.
     * 
     * @param color color to be given to the line.     
     */
    int c(std::string color);

    /**
     * @brief Path to witch the figure will be saved.
     * 
     * If the path includes not existing directories those will be created without an explicit notice.
     * 
     * @param filepath Path to witch the figure will be saved.     
     */
    int fpath(std::string filepath);

    /**
     * @brief Extension of the file in witch the figure will be saved.
     * 
     * This can be either "pdf" (default one) or "jpg".
     * 
     * @param filetype "jpg" or "pdf" (default).     
     */
    int fext(std::string filetype);

    /**
     * @brief Sets limits for the x axis.
     * 
     * @param min Bottom limit.
     * @param max Upper limit.     
     */
    int limx(double min, double max);

    /**
     * @brief Sets Upper limit for the x axis.
     * 
     * If it is necessary to set only the upper limit for the x axis (thus leaving the bottom limit the 
     * default one) then this function can be called with "" as the @latexonly \textit{empty} @endlatexonly
     * argument. 
     * 
     * @param empty ""
     * @param max Upper limit.     
     */
    int limx(std::string empty, double max);

    /**
     * @brief Sets Bottom for the x axis.
     * 
     * If it is necessary to set only the Bottom limit for the x axis (thus leaving the Upper limit the 
     * default one) then this function can be called with "" as the @latexonly \textit{empty} @endlatexonly
     * argument. 
     * 
     * @param empty ""
     * @param min Bottom limit.     
     */
    int limx(double min, std::string empty);

    /**
     * @brief Same as limx functions but for the y axis. 
     */
    int limy(double min, double max);

    /**
     * @brief Same as limx functions but for the y axis. 
     */
    int limy(std::string empty, double max);

    /**
     * @brief Same as limx functions but for the y axis. 
     */
    int limy(double min, std::string empty);

    /**
     * @brief Same as limx functions but for the z axis. 
     */
    int limz(double min, double max);

    /**
     * @brief Same as limx functions but for the z axis. 
     */
    int limz(std::string empty, double max);

    /**
     * @brief Same as limx functions but for the z axis. 
     */
    int limz(double min, std::string empty);

    /**
     * @brief Sets limits to the considered data in a fit.
     * 
     * When fitting a set of data must be given to fit on. This option lets the user limit the range
     * of data considered in the fit.
     * 
     * @param min Bottom limit.
     * @param max Upper limit.     
     */
    int limxfit(double min, double max);

    /**
     * @brief Sets upper limit to the considered data in a fit.
     * 
     * If is necessary to consider all the data given to fit up to a certain value this option can be used
     * with "" as the @latexonly \textit{empty} @endlatexonly argument.
     * 
     * @param empty ""
     * @param max Upper limit.
     */
    int limxfit(std::string empty, double max);

    /**
     * @brief Sets bottom limit to the considered data in a fit.
     * 
     * If is necessary to consider all the data given to fit from a certain value this option can be used
     * with "" as the @latexonly \textit{empty} @endlatexonly argument.
     * 
     * @param empty ""
     * @param min Bottom limit.
     */
    int limxfit(double min, std::string empty);

    /**
     * @brief Same as limxfit function but for the y axis.  
     */
    int limyfit(double min, double max);

    /**
     * @brief Same as limxfit function but for the y axis.  
     */
    int limyfit(std::string empty, double max);

    /**
     * @brief Same as limxfit function but for the y axis.  
     */
    int limyfit(double min, std::string empty);

    /**
     * @brief Specify a function to be overlayed to the figure.
     * 
     * The function must be inputted with the gnuplot syntax (e.g.: "f(x)=3*x+4" or "f(x,y)=12*x**2+3*y**2")
     * spaces in the declaration of the function will probably badly break the program :(.
     * 
     * @param function function that will be overlayed to the figure.     
     */
    int func(std::string function);

    /**
     * @brief Specify style of the line.
     * 
     * @param linestyle either "points", "line" or "linespoints" (default).
     */
    int ls(std::string linestyle);

    /**
     * @brief Specify the title of the line
     * 
     * If multiple lines are to be plotted then this title applies only to the first line plotted after the
     * GnuplotDriver.ls() call.
     * 
     * @param linetitle Title of the line.
     */
    int lt(std::string linetitle);
    int keyloc(std::string keylocation);
    int print(std::string dataout_path);
    int h();
    int p();
    int m();
    int logx();
    int logy();
    int nogrid();
    int noprint();
    int plot(double * dati, int numRighe, int numColonne);
    int plot(double * dati1, double * dati2, int numRighe);
    int plot(double * dati1, double * dati2, double * dati3, int numRighe);
    int fit(double * dati, int numRighe, int numColonne, std::string funzione, std::string parametri);
    int fit(double * dati1, double * dati2, int numRighe, std::string funzione, std::string parametri);
    int fit(double * dati1, double * dati2, double * dati3, int numRighe, std::string funzione, std::string parametri);
    private:

    /**
     * @brief substitues a substring with another in a bigger string
     * 
     * @param input string in witch the sustitution is effectuated
     * @param toSub substring to be substituted
     * @param substitute substring to put in place of toSub
     * @return int, 0 if the fuction exits correctly
     */
    int sub(std::string * input, std::string toSub, std::string substitute);

    /**
     * @brief creates the required directory to store output files
     * 
     * @param path the required path, works with and without the name of the file included
     * @return int: 0 if path is just a file name, 1 if the required directories already existed or were successfully created and -1 otherwise
     */
    int makepath(std::string path);
    bool config_stream();
    int plotcommand();
    int fitcommand(std::string funzione, std::string parametri);

};
