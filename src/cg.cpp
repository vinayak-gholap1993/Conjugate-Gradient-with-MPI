#include "cg.h"

//r   f - A  u
// delta0 =  r * r
// d = r
void cg::computeFirstStep(const Grid& u,const Grid& rhs,Grid& res,Grid& d,real& para_delta0,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			res(i,j) = rhs(i,j) + para_invHx * (u(i,j-1) + u(i,j+1)) + para_invHy * (u(i-1,j) + u(i+1,j)) - para_weightC*u(i,j);
			d(i,j) = res(i,j);
			para_delta0 += res(i,j) * res(i,j);
		}
	}
}

// z = A * d
// alpha = d * A * d or d * z
void cg::multiplyVector(const Grid& d,Grid& z,real& para_alpha,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			z(i,j) =  (-1) * para_invHx * (d(i,j-1) + d(i,j+1)) - para_invHy * (d(i-1,j) + d(i+1,j)) + para_weightC*d(i,j);
			para_alpha += d(i,j) * z(i,j);
		}
	}
}

void cg::updateU(Grid& u,const Grid& d,const real& para_alpha,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			u(i,j) =  u(i,j) + para_alpha * d(i,j);
		}
	}
}

void cg::updateRes(Grid& res,const Grid& z,const real& para_alpha,real& para_delta1,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			res(i,j) =  res(i,j) - para_alpha * z(i,j);
			para_delta1 += res(i,j) * res(i,j);
		}
	}
}

void cg::updateD(Grid& d,const Grid& res,const real& para_beta,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			d(i,j) =  res(i,j) + para_beta * d(i,j);
		}
	}
}

//output : number in grid x / processor , same as y
void cg::findElementsinXY(int& elementX, int& elementY, const uint& nx, const uint& ny,const int& numProcX,const int& numProcY,const int& cartCordX,const int& cartCordY)
{
    elementX = nx / numProcX;
    if(cartCordX == numProcX-1)
    {
      elementX += nx % numProcX;
     // std::cout<<"coord in X "<<cartCordX<<"\t coord in Y "<<cartCordY<<"\t elementX "<<elementX<<"\t elementY  "<<elementY<<std::endl;
    }
    
    elementY = ny / numProcY;
    if(cartCordY == numProcY-1)
    {
      elementY += ny % numProcY;
      //std::cout<<"coord in X "<<cartCordX<<"\t coord in Y "<<cartCordY<<"\t elementX "<<elementX<<"\t elementY  "<<elementY<<std::endl;
    }
}

void cg::findGlobalCoord(int& xStart,int& xEnd,int& yStart,int& yEnd,const int& elementX,const int& elementY,const int& mycoordX,const int& mycoordY,const int& numProcX,const int& numProcY,const uint& nx,const uint& ny)
{
    if(mycoordX == numProcX-1)
      {
	xStart = mycoordX * (elementX - (nx % numProcX));
	xEnd = (xStart + elementX)-1;
      }
    else
      {
	xStart = mycoordX * elementX;
	xEnd = xStart + elementX;
      }
  //-------------------------------------------------------------------------------------
    if(mycoordY == numProcY-1)
      {
	yStart = mycoordY * (elementY - (ny % numProcY));
	yEnd = (yStart + elementY)-1;
      }
    
    else
      {
	yStart = mycoordY * elementY;
	yEnd = yStart + elementY; 
      }
//we don't want to iterate at boundary points,so,
    if(mycoordX == 0)
    {
      xStart += 1;
    }
    /*
    if(mycoordX == numProcX-1)
    {
      xEnd -= 1;
    }*/
    
    if(mycoordY == 0)
    {
      yStart += 1;
    }
    /*
    if(mycoordY == numProcY-1)
    {
      yEnd -= 1;
    }*/
      
}
