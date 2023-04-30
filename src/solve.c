#include "fem.h"

void femMeshLocal(const femMesh *theMesh, const int iElem, int *map, double *x, double *y){
    int j,nLocal = theMesh->nLocalNode;
    for (j=0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j]   = theMesh->nodes->X[map[j]];
        y[j]   = theMesh->nodes->Y[map[j]]; }
}

double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,i,j,map[4];

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femMeshLocal(theMesh,iElem,map,x,y);
        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i]*dphidxsi[i];
                dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i];
                dydeta += y[i]*dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j]
                                                     + dphidy[i] * dphidy[j]) * jac * weight;
                    /*A[2*map[i]][2*map[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight ;
                    A[2*map[i]][2*map[j]+1] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[2*map[i]+1][2*map[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[2*map[i]+1][2*map[j]+1] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;*/

                    A[2*map[i]][2*map[j]] += (dphidx[i] * a * dphidx[j] * x[i] + dphidy[i] * c * dphidy[j] * x[i] + b * dphidx[i] * phi[i] + phi[i]*(b * dphidx[j] + (a * phi[j])/x[i] )) * jac * weight ;
                    A[2*map[i]][2*map[j]+1] += (dphidx[i] * b * dphidy[j] * x[i] + dphidy[i] * c * dphidx[j] * x[i] + b * dphidy[i] * phi[j]) * jac * weight;
                    A[2*map[i]+1][2*map[j]] += (dphidy[i] * b * dphidx[j] * x[i] + dphidx[i] * c * dphidy[j] * x[i] + b * phi[i] * dphidy[j]) * jac * weight;
                    A[2*map[i]+1][2*map[j]+1] += (dphidy[i] * a * dphidy[j] * x[i] + dphidx[i] * c * dphidx[j] * x[i]) * jac * weight;

                }
            }
            for (i = 0; i < theSpace->n; i++) {
                B[2*map[i]] += 0;
                B[2*map[i]+1] -= phi[i] * jac * weight * g * rho;
            }
        }
    }


  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
