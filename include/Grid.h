
#ifndef GRID_H
#define GRID_H


#include"type.h"

class Grid
{
public:
	//data member
	real *data;
	uint nx , ny;

	//function member
	Grid(void){;}
	~Grid(void){;}

	Grid(const uint& para_totalGrid);
	Grid(const uint& para_nx,const uint& para_ny);

	//void applyBoundaryU(Grid& para_u, const uint& para_xStart, const uint& para_xEnd,const uint& para_yStart, const uint& para_yEnd, const real& para_hx);
	void applyBoundaryU(Grid& para_u, const uint& para_xStart, const uint& para_xEnd,const uint& para_yStart, const uint& para_yEnd,const real& para_hx,const int& numProcY,const int& mycoordX,const int& mycoordY);
	void applyBoundaryRHS(Grid& para_rhs, const uint& para_xStart, const uint& para_elementX,const uint& para_yStart, const uint& para_elementY,const real& para_hx,const real& para_hy,const int& mycoordX);
	void print(const Grid& para_print , const uint& para_nx , const uint& para_ny);

	real& operator()(const uint& x,const uint& y);
	real operator()(const uint& x,const uint& y) const;
	
	//void writeFile(const std::string& fileName, const Grid& para_storeToFile,const uint& para_nx,const uint& para_ny,const real& para_hx,const real& para_hy);
};

#endif