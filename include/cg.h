
#ifndef CG_H
#define CG_H

#include"type.h"
#include"Grid.h"

class cg
{
public:
	cg(void){;}
	~cg(void){;}

	void computeFirstStep(const Grid& u,const Grid& rhs,Grid& res,Grid& d,real& para_delta0,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC);
        void multiplyVector(const Grid& d,Grid& z,real& para_alpha,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC);

	void updateU(Grid& u, const Grid& d, const real& para_alpha, const uint& para_nx, const uint& para_ny);
	void updateRes(Grid& res, const Grid& z, const real& para_alpha, real& para_delta1, const uint& para_nx, const uint& para_ny);
	void updateD(Grid& d,const Grid& res,const real& para_beta,const uint& para_nx,const uint& para_ny);

	void findElementsinXY(int& elementX, int& elementY, const uint& nx, const uint& ny,const int& numProcX,const int& numProcY,const int& cartCordX,const int& cartCordY);
	void findGlobalCoord(int& xStart,int& xEnd,int& yStart,int& yEnd,const int& elementX,const int& elementY,const int& mycoordX,const int& mycoordY,const int& numProcX,const int& numProcY,const uint& nx,const uint& ny);

	void sendReceiveMsg(Grid& u,int& elementX,int& elementY,MPI_Comm& GRID_COMM_CART,MPI_Datatype& colType,MPI_Datatype& rowType);
	
};

#endif