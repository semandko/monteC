// SiliconOxideMatrix.h
#ifndef SILICONOXIDEMATRIX_H
#define SILICONOXIDEMATRIX_H

#include <vector>
#include <fstream>


#define _CONST_A 1.38e-23
#define _CONST_B 1.602e-19
#define _K_B (_CONST_A / _CONST_B)
#define _X_0 0.1
#define _T 1000
#define _kT_EV(_T) (1.0 / 11604 * (_T + 273.15)) // 0.11 eV for 1000 C
#define _NUM_ROWS 30
#define _NUM_COLS 30

struct OxygenPosition
{
    bool valid;
    int row;
    int column;
};


class SiliconOxideMatrix {
	
public:
	
	SiliconOxideMatrix();
	~SiliconOxideMatrix();
    
    void fillMatrix();
    double Penalty(int i, int j);
    void initializeArrayO();
    void updateArrayOLoop(std::vector<std::pair<int, int>>& ArrayO);
    void updateDataMatrix(std::vector<std::pair<int, int>>& ArrayO);

    //
    void calculateNumberOfOxygenAtomsByInputPercentage(const int percentage);
    void putOxygenAtomsOnRdmPlaces();
    void debug_getNumber();
    //

    int calculationNumOxigen(int i, int j);
    void printMatrixToFile(const std::string& fileName);
    void printMatrixToImage(const std::string& fileName);
    bool metropolisCondition(double a, double b);
    double randomGenerator(unsigned int first_interval, unsigned int last_interval);
    void evolution(void);
    bool isOxygenInCell(int i, int j);
    bool findingOxigeninCeil(int &index_i, int & index_j);
	 
private:
    int maxi;
    int maxj;
    std::vector<double> Delta{0.0, 0.5, 0.51, 0.22, 0.0};
    std::vector<std::vector<double>> data;
    std::vector<std::pair<int, int>> ArrayO;

    //
    std::vector<std::pair<int, int>> allPossibleOxygenPlaces;
    std::vector<OxygenPosition> allOxygenPossitions;
    int generalNumberOfAtoms;
    int numberOfOxygenAtoms;
    //

    double penaltyValue;
    double x0 = _X_0;
    double kB = _K_B;
    double kT_eV = _kT_EV(_T);
    int T = _T;
    int maxIterations = static_cast<int>((maxi * maxj * x0) / 4);
        	  
};

#endif

