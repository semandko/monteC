// SiliconOxideMatrix.h
#ifndef SILICONOXIDEMATRIX_H
#define SILICONOXIDEMATRIX_H

#include <vector>
#include <fstream>
#include <mutex>


#define _CONST_A 1.38e-23
#define _CONST_B 1.602e-19
#define _K_B (_CONST_A / _CONST_B)
#define _X_0 0.1
#define _T 1000
#define _kT_EV(_T) (1.0 / 11604 * (_T + 273.15)) // 0.11 eV for 1000 C
#define _NUM_ROWS 10
#define _NUM_COLS 10

struct OxygenPos
{
    int row;
    int column;
    bool valid;
};

struct MatrixPossition
{
    double value;
    bool valid;
};


class SiliconOxideMatrix 
{
	
public:
	
	SiliconOxideMatrix();
	~SiliconOxideMatrix();
    
    void fillMatrix();
    double calculatePenaltyForKremniy(int i, int j);

    //
    void calculateNumberOfOxygenAtomsByInputPercentage(const int percentage);
    int getRdmIntNumber(int start, int notIncludedEnd);
    double getRdmDoubleNumber(int start, int notIncludedEnd);
    void putOxygenAtomsOnRdmPlaces();
    void evolutionFindOxygen();
    bool checkingCellOccupationAndJumping(int row, int column, bool method);
    std::vector<std::pair<int, int>> getJumpPossitions(int row, int column, bool method);
    void changeValidationToV(int row, int column, bool method);
    void debug_getNumber();
    int debug_getNumber_test();
    //

    int calculationNumberOfOxigens(int i, int j);
    void printMatrixToFile(const std::string& fileName);
    // void printMatrixToImage(const std::string& fileName);
    bool metropolisCondition(double a, double b);
	 
private:
    int N;
    std::vector<double> Delta{0.0, 0.5, 0.51, 0.22, 0.0};
    std::vector<std::vector<MatrixPossition>> data;

    //
    std::vector<OxygenPos> allPossibleOxygenPlaces;
    std::vector<OxygenPos> allOxygenPossitions;

    std::vector<std::pair<int, int>> horizontalKremniyPlacesOfset;
    std::vector<std::pair<int, int>> verticalKremniyPlacesOfset;
    std::vector<std::pair<int, int>> horizontalOxygenPlacesOfset;
    std::vector<std::pair<int, int>> verticalOxygenPlacesOfset;

    int generalNumberOfAtoms;
    int numberOfOxygenAtoms;
    long long evo_iterations;

    std::mutex iterations_m;
    std::mutex valid_m;
    //

    double penaltyValue;
    double x0 = _X_0;
    double kB = _K_B;
    double kT_eV = _kT_EV(_T);
    int T = _T;
};

#endif

