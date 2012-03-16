// 
//  ----------------------------------------------------------------------------
//  "THE BEER-WARE LICENSE" (invented by Poul-Henning Kamp):
//  Anders Nielsen <an@aqua.dtu.dk> wrote this file. As long as you retain this 
//  notice you can do whatever you want with this stuff. If we meet some day, 
//  and you think this stuff is worth it, you can buy me a beer in return. 
//  ----------------------------------------------------------------------------
//

// Anders Nielsen <anders.nielsen@hawaii.edu> Sep 2005 

#ifndef __nLogNormal_h__
#define __nLogNormal_h__

#include <fvar.hpp>
#include <math.h>

_CONST double log2pi = log(2.0*M_PI);


df1b2matrix bksolve(_CONST df1b2matrix& L, _CONST df1b2matrix& b)
{
  RETURN_ARRAYS_INCREMENT();
  int i, k, r, R, c, C;
  r=b.rowmin();
  R=b.rowmax();
  c=b.colmin();
  C=b.colmax();

  df1b2vector sumVec(c,C);
  df1b2matrix x(r,R,c,C);

  for(i=r; i<=R; ++i){
    sumVec=b(i);
    for(k=i-1; k>=r; k--){
      sumVec-=L(i,k)*x(k);
    }
    x(i)=sumVec/L(i,i);
  }
  RETURN_ARRAYS_DECREMENT();
  return x;
}

dvar_matrix bksolve(_CONST dvar_matrix& L, _CONST dvar_matrix& b)
{
  RETURN_ARRAYS_INCREMENT();
  int i, k, r, R, c, C;
  r=b.rowmin();
  R=b.rowmax();
  c=b.colmin();
  C=b.colmax();

  dvar_vector sumVec(c,C);
  dvar_matrix x(r,R,c,C);

  for(i=r; i<=R; ++i){
    sumVec=b(i);
    for(k=i-1; k>=r; k--){
      sumVec-=L(i,k)*x(k);
    }
    x(i)=sumVec/L(i,i);
  }
  RETURN_ARRAYS_DECREMENT();
  return x;
}

df1b2vector nLogNormal(_CONST df1b2vector& x, _CONST df1b2matrix& mu, _CONST df1b2matrix& S)
{
  RETURN_ARRAYS_INCREMENT();	
  int r, R, c, C, N;
  r=mu.rowmin();
  R=mu.rowmax();
  c=mu.colmin();
  C=mu.colmax();
  N=R-r+1;
  df1b2vector ret(c,C);
  df1b2matrix diff(r,R,c,C);
  for(int i=c; i<=C; ++i){
    for(int j=r; j<=R; ++j){
      diff(j,i)=x(j)-mu(j,i);
    }
  }
  df1b2variable logDet=0.0;
  df1b2matrix chol=choleski_decomp(S);
  for(int i=r; i<=R; ++i){logDet+=log(chol(i,i));}
  logDet*=2.0;
  df1b2matrix tmp=bksolve(chol,diff);
  ret=0.5*(log2pi*N+logDet+colsum(square(tmp)));
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

dvar_vector nLogNormal(_CONST dvar_vector& x, _CONST dvar_matrix& mu, _CONST dvar_matrix& S)
{
  RETURN_ARRAYS_INCREMENT();	
  int r, R, c, C, N;
  r=mu.rowmin();
  R=mu.rowmax();
  c=mu.colmin();
  C=mu.colmax();
  N=R-r+1;
  dvar_vector ret(c,C);
  dvar_matrix diff(r,R,c,C);
  for(int i=c; i<=C; ++i)diff.colfill(i,x-column(mu,i));
  dvariable logDet;
  dvar_matrix chol=choleski_decomp(S);
  logDet=2.0*sum(log(diagonal(chol)));
  dvar_matrix tmp=bksolve(chol,diff);
  ret=0.5*(log2pi*N+logDet+colsum(square(tmp)));
  RETURN_ARRAYS_DECREMENT();
  return ret;
}


df1b2vector nLogNormal(_CONST df1b2matrix& x, _CONST df1b2vector& mu, _CONST df1b2matrix& S)
{
  RETURN_ARRAYS_INCREMENT();
  df1b2vector ret=nLogNormal(mu, x, S);
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

dvar_vector nLogNormal(_CONST dvar_matrix& x, _CONST dvar_vector& mu, _CONST dvar_matrix& S)
{
  RETURN_ARRAYS_INCREMENT();
  dvar_vector ret=nLogNormal(mu, x, S);
  RETURN_ARRAYS_DECREMENT();
  return ret;
}


df1b2variable nLogNormal(_CONST df1b2vector& x, _CONST df1b2vector& mu, _CONST df1b2matrix& S)
{
  RETURN_ARRAYS_INCREMENT();	
  int r=mu.indexmin(), R=mu.indexmax();
  df1b2matrix MU(r,R,1,1);
  for(int i=r; i<=R; ++i){MU(i,1)=mu(i);}
  df1b2vector tmp=nLogNormal(x, MU, S);
  df1b2variable ret=tmp(tmp.indexmin());	
  RETURN_ARRAYS_DECREMENT();
  return ret;
}
dvariable nLogNormal(_CONST dvar_vector& x, _CONST dvar_vector& mu, _CONST dvar_matrix& S)
{
  RETURN_ARRAYS_INCREMENT();	
  dvar_matrix MU(mu.indexmin(),mu.indexmax(),1,1);
  MU.colfill(1,mu);
  dvar_vector tmp=nLogNormal(x, MU, S);
  dvariable ret=tmp(tmp.indexmin());	
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

//TEST
//  dvar_vector nLogNormalChol(_CONST dvar_vector& x, _CONST dvar_matrix& mu, _CONST dvar_matrix& chol)
//  {
//    RETURN_ARRAYS_INCREMENT();	
//    int r, R, c, C, N;
//    r=mu.rowmin();
//    R=mu.rowmax();
//    c=mu.colmin();
//    C=mu.colmax();
//    N=R-r+1;
//    dvar_vector ret(c,C);
//    dvar_matrix diff(r,R,c,C);
//    for(int i=c; i<=C; ++i)diff.colfill(i,x-column(mu,i));
//    dvariable logDet;
//    logDet=2.0*sum(log(diagonal(chol)));
//    dvar_matrix tmp=bksolve(chol,diff);
//    ret=0.5*(log2pi*N+logDet+colsum(square(tmp)));	
//    RETURN_ARRAYS_DECREMENT();
//    return ret;
//  }
//  
//  dvar_vector nLogNormalChol(_CONST dvar_matrix& x, _CONST dvar_vector& mu, _CONST dvar_matrix& chol)
//  {
//    RETURN_ARRAYS_INCREMENT();	
//    dvar_vector ret=nLogNormalChol(mu, x, chol);	
//    RETURN_ARRAYS_DECREMENT();
//    return ret;
//  }
//  
//  dvariable nLogNormalChol(_CONST dvar_vector& x, _CONST dvar_vector& mu, _CONST dvar_matrix& chol)
//  {
//    RETURN_ARRAYS_INCREMENT();	
//    dvar_matrix MU(mu.indexmin(),mu.indexmax(),1,1);
//    MU.colfill(1,mu);
//    dvar_vector tmp=nLogNormalChol(x, MU, chol);
//    dvariable ret=tmp(tmp.indexmin());	
//    RETURN_ARRAYS_DECREMENT();
//    return ret;
//  }
//  
//  void dmdm_fwsolve(void);
//  
//  dvar_matrix fwsolve(_CONST dvar_matrix& L, _CONST dvar_matrix& b)
//  {  
//    RETURN_ARRAYS_INCREMENT();	
//    int i, k, r, R, c, C;
//    r=b.rowmin();
//    R=b.rowmax();
//    c=b.colmin();
//    C=b.colmax();
//  
//    dvector sumVec(c,C);
//    dmatrix xVal(r,R,c,C);
//    dmatrix LVal=value(L);
//    dmatrix bVal=value(b);
//  
//    for(i=R; i>=r; --i){
//      sumVec=bVal(i); 
//      for(k=i+1; k<=R; k++){
//        sumVec-=LVal(k,i)*xVal(k);
//      }
//      xVal(i)=sumVec/LVal(i,i);
//    }
//  
//    dvar_matrix x=nograd_assign(xVal);
//    save_identifier_string("TEST1");
//    L.save_dvar_matrix_value();
//    L.save_dvar_matrix_position();
//    b.save_dvar_matrix_value();
//    b.save_dvar_matrix_position();
//    x.save_dvar_matrix_value();
//    x.save_dvar_matrix_position();
//    save_identifier_string("TEST6");
//    gradient_structure::GRAD_STACK1->
//             set_gradient_stack(dmdm_fwsolve);	
//    RETURN_ARRAYS_DECREMENT();
//    return x;
//  }
//  
//  void dmdm_fwsolve(void)
//  {
//    verify_identifier_string("TEST6");
//    dvar_matrix_position xpos=restore_dvar_matrix_position();
//    dmatrix dfx=restore_dvar_matrix_derivatives(xpos);
//    dmatrix x=restore_dvar_matrix_value(xpos);
//    dvar_matrix_position bpos=restore_dvar_matrix_position();
//    dmatrix b=restore_dvar_matrix_value(bpos);
//    dvar_matrix_position Lpos=restore_dvar_matrix_position();
//    dmatrix L=restore_dvar_matrix_value(Lpos);
//    verify_identifier_string("TEST1");
//  
//    int i, k, r, R, c, C;
//    r=b.rowmin();  
//    R=b.rowmax();  
//    c=b.colmin();
//    C=b.colmax();
//  
//    dmatrix dfL(Lpos);
//    dmatrix dfb(bpos);
//    dvector sumVec(c,C);
//    dvector dfsumVec(c,C);
//    dfL.initialize();
//    dfb.initialize();
//    dfsumVec.initialize();
//  
//  
//    for(i=r; i<=R; ++i){
//      //bring back sumVec 
//      sumVec=x(i)*L(i,i);
//  
//      //xVal(i)=sumVec/LVal(i,i);
//      dfsumVec+=dfx(i)/L(i,i);
//      dfL(i,i)+= -(dfx(i)*sumVec)/(L(i,i)*L(i,i)); 
//      dfx(i)=0.0;
//  
//      for(k=R; k>=i+1; k--){
//        //sumVec-=LVal(k,i)*xVal(k);
//        dfL(k,i)-=dfsumVec*x(k);
//        dfx(k)-=dfsumVec*L(k,i);
//      }
//      //sumVec=bVal(i); 
//      dfb(i)+=dfsumVec; 
//      dfsumVec=0.0;
//    }
//  
//    dfL.save_dmatrix_derivatives(Lpos);
//    dfb.save_dmatrix_derivatives(bpos);
//  }
//  
//  
//  
#endif



 
