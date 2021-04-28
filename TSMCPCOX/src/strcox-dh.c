#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>



SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
    return elmt;
}

// Cross product of y with jth column of X
double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of squares of jth column of X
double sqsum(double *X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += pow(X[nn+i], 2);
  return(val);
}

double sum(double *x, int n) {
  double val=0;
  for (int i=0;i<n;i++) val += x[i];
  return(val);
}

double MCP(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}

double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}

double lasso(double z, double l1, double l2, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else return(s*(fabs(z)-l1)/(v*(1+l2)));
}


// Memory handling, output formatting (Cox)
SEXP maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_) {
  
  // Declarations
  int n = nrows(X_);
  int p = length(v_);
  SEXP z;
  PROTECT(z = allocVector(REALSXP, 1));
  REAL(z)[0] = 0;
  double zz;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *v = INTEGER(v_);
  
  for (int j=0; j<p; j++) {
    zz = crossprod(X, y, n, v[j]-1) / m[v[j]-1];
    if (fabs(zz) > REAL(z)[0]) REAL(z)[0] = fabs(zz);
  }
  
  // Return list
  UNPROTECT(1);
  return(z);
}

SEXP cleanupCox(double *a, double *r, double *h, int *e, double *eta, double *haz, double *rsk, SEXP beta, SEXP Loss, SEXP iter, SEXP Eta) {
  Free(a);
  Free(r);
  Free(h);
  Free(e);
  Free(eta);
  Free(haz);
  Free(rsk);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, Loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, Eta);
  UNPROTECT(1);
  return(res);
}
// Coordinate descent for Cox models
SEXP cdfit_cox_dh(SEXP X_, SEXP d_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {

  // Lengths/dimensions
  int n = length(d_);
  int p = length(X_)/n;
  int L = length(lambda);

  // Pointers
  double *X = REAL(X_);
  double *d = REAL(d_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];
  int warn = INTEGER(warn_)[0];

  // Outcome
  SEXP res, beta, Loss, iter, Eta;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(Loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;

  // Intermediate quantities
  double *a = Calloc(p, double);    // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double *haz = Calloc(n, double);
  double *rsk = Calloc(n, double);
  double *r = Calloc(n, double);
  double *h = Calloc(n, double);
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double *eta = Calloc(n, double);
  for (int i=0; i<n; i++) eta[i] = 0;
  double xwr, xwx, u, v, l1, l2, shift, si, s, nullDev;
  int lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  rsk[n-1] = 1;
  for (int i=n-2; i>=0; i--) rsk[i] = rsk[i+1] + 1;
  nullDev = 0;
  for (int i=0; i<n; i++) nullDev -= d[i]*log(rsk[i]);
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Loss)[0] = nullDev;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a[j] != 0) nv++;
      }
      if ((nv > dfmax) | (tot_iter == max_iter)) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }
    }

    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        INTEGER(iter)[l]++;
        tot_iter++;
        REAL(Loss)[l] = 0;
        double maxChange = 0;

        // Calculate haz, risk
        for (int i=0; i<n; i++) haz[i] = exp(eta[i]);
        rsk[n-1] = haz[n-1];
        for (int i=n-2; i>=0; i--) {
          rsk[i] = rsk[i+1] + haz[i];
        }
        for (int i=0; i<n; i++) {
          REAL(Loss)[l] += d[i]*eta[i] - d[i]*log(rsk[i]);
        }

        // Approximate L
        h[0] = d[0]/rsk[0];
        for (int i=1; i<n; i++) {
          h[i] = h[i-1] + d[i]/rsk[i];
        }
        for (int i=0; i<n; i++) {
          h[i] = h[i]*haz[i];
          s = d[i] - h[i];
          if (h[i]==0) r[i]=0;
          else r[i] = s/h[i];
        }

        // Check for saturation
        if (REAL(Loss)[l]/nullDev < .01) {
          if (warn) warning("Model saturated; exiting...");
          for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }

        // Covariates
        for (int j=0; j<p; j++) {
          if (e[j]) {

            // Calculate u, v
            xwr = wcrossprod(X, r, h, n, j);
            xwx = wsqsum(X, h, n, j);
            u = xwr/n + (xwx/n)*a[j];
            v = xwx/n;

            // Update b_j
            l1 = lam[l] * m[j] * alpha;
            l2 = lam[l] * m[j] * (1-alpha);
            if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(u, l1, l2, gamma, v);
            if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(u, l1, l2, gamma, v);
            if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(u, l1, l2, v);

            // Update r
            shift = b[l*p+j] - a[j];
            if (shift !=0) {
              for (int i=0;i<n;i++) {
                si = shift*X[j*n+i];
                r[i] -= si;
                eta[i] += si;
              }
              if (fabs(shift)*sqrt(v) > maxChange) maxChange = fabs(shift)*sqrt(v);
            }
          }
        }

        // Check for convergence
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }

      // Scan for violations
      int violations = 0;
      for (int j=0; j<p; j++) {
        if (e[j]==0) {
          xwr = wcrossprod(X, r, h, n, j)/n;
          l1 = lam[l] * m[j] * alpha;
          if (fabs(xwr) > l1) {
            e[j] = 1;
            violations++;
          }
        }
      }
      if (violations==0) {
        for (int i=0; i<n; i++) {
          REAL(Eta)[l*n+i] = eta[i];
        }
        break;
      }
    }
  }
  res = cleanupCox(a, r, h, e, eta, haz, rsk, beta, Loss, iter, Eta);
  UNPROTECT(4);
  return(res);
}

  
SEXP standardize(SEXP X_) {
  // Declarations
  int n = nrows(X_);
  int p = ncols(X_);
  SEXP XX_, c_, s_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  PROTECT(c_ = allocVector(REALSXP, p));
  PROTECT(s_ = allocVector(REALSXP, p));
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *c = REAL(c_);
  double *s = REAL(s_);
  
  for (int j=0; j<p; j++) {
    // Center
    c[j] = 0;
    for (int i=0; i<n; i++) {
      c[j] += X[j*n+i];
    }
    c[j] = c[j] / n;
    for (int i=0; i<n; i++) XX[j*n+i] = X[j*n+i] - c[j];
    
    // Scale
    s[j] = 0;
    for (int i=0; i<n; i++) {
      s[j] += pow(XX[j*n+i], 2);
    }
    s[j] = sqrt(s[j]/n);
    for (int i=0; i<n; i++) XX[j*n+i] = XX[j*n+i]/s[j];
  }
  
  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_);
  SET_VECTOR_ELT(res, 1, c_);
  SET_VECTOR_ELT(res, 2, s_);
  UNPROTECT(4);
  return(res);
}
SEXP strat_cox_dh(SEXP X_, SEXP d_,SEXP qn_, SEXP mn_, SEXP fn_,SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {
  
  // Lengths/dimensions
  int n = length(d_);
  int p = length(X_)/n;
  int L = length(lambda);
  
  // Pointers
  double *X = REAL(X_);
  double *d = REAL(d_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int qn = INTEGER(qn_)[0];
  int mn = INTEGER(mn_)[0];
  int fn = INTEGER(fn_)[0];
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];
  int warn = INTEGER(warn_)[0];
  
  // Outcome
  SEXP res, beta, Loss, iter, Eta;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(Loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;
  
  // Intermediate quantities
  double *a = Calloc(p, double);    // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double *haz = Calloc(n, double);
  double *rsk = Calloc(n, double);
  double *r = Calloc(n, double);
  double *h = Calloc(n, double);
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double *eta = Calloc(n, double);
  for (int i=0; i<n; i++) eta[i] = 0;
  double xwr, xwx, u, v, l1, l2, shift, si, s, nullDev;
  int lstart;
  
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  for(int i=0; i<qn; i++) rsk[fn+mn*i-1]=1;
  for(int i=fn-2; i>=0; i--) rsk[i]=rsk[i+1]+1;
  
  for(int j=0; j< (qn-1); j++){
    for (int i=mn-2; i>=0; i--){
      rsk[fn+mn*j+i] = rsk[fn+mn*j+i+1] + 1;
    }
  }   
  nullDev = 0;
  for (int i=0; i<n; i++) nullDev -= d[i]*log(rsk[i]);
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Loss)[0] = nullDev;
  }
  
  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];
      
      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a[j] != 0) nv++;
      }
      if ((nv > dfmax) | (tot_iter == max_iter)) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }
    }
    
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        INTEGER(iter)[l]++;
        tot_iter++;
        REAL(Loss)[l] = 0;
        double maxChange = 0;
        
        // Calculate haz, risk
        for (int i=0; i<n; i++) haz[i] = exp(eta[i]);
        for(int i=0; i<qn; i++) rsk[fn+mn*i-1]=haz[fn+mn*i-1];
        for(int i=fn-2; i>=0; i--) rsk[i]=rsk[i+1]+haz[i];
        
        for(int j=0; j< (qn-1); j++){
          for (int i= (mn-2); i>=0; i--){
            rsk[fn+mn*j+i] = rsk[fn+mn*j+i+1] + haz[fn+mn*j+i];
          }
        }   
        for (int i=0; i<n; i++) {
          REAL(Loss)[l] += d[i]*eta[i] - d[i]*log(rsk[i]);
        }
        // Approximate L
        
       
       h[0] = d[0]/rsk[0];
        for(int i=0; i<(qn-1); i++) h[fn+mn*i]=d[fn+mn*i]/rsk[fn+mn*i];
        for(int i=1; i<fn; i++) h[i]=h[i-1]+ d[i]/rsk[i];
        for(int j=0; j< (qn-1); j++){
          for (int i= 1; i<mn; i++){
            h[fn+mn*j+i] = h[fn+mn*j+i-1] + d[fn+mn*j+i]/rsk[fn+mn*j+i];
          }
        }   
        
        for (int i=0; i<n; i++) {
          h[i] = h[i]*haz[i];
          s = d[i] - h[i];
          if (h[i]==0) r[i]=0;
          else r[i] = s/h[i];
        }
        
       /* h[0] = d[0]/rsk[0];
        for (int i=1; i<n; i++) {
          h[i] = h[i-1] + d[i]/rsk[i];
        }
        for (int i=0; i<n; i++) {
          h[i] = h[i]*haz[i];
          s = d[i] - h[i];
          if (h[i]==0) r[i]=0;
          else r[i] = s/h[i];
        }*/
        // Check for saturation
        if (REAL(Loss)[l]/nullDev < .01) {
          if (warn) warning("Model saturated; exiting...");
          for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }
        
        // Covariates
        for (int j=0; j<p; j++) {
          if (e[j]) {
            
            // Calculate u, v
            xwr = wcrossprod(X, r, h, n, j);
            xwx = wsqsum(X, h, n, j);
            u = xwr/n + (xwx/n)*a[j];
            v = xwx/n;
            
            // Update b_j
            l1 = lam[l] * m[j] * alpha;
            l2 = lam[l] * m[j] * (1-alpha);
            if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(u, l1, l2, gamma, v);
            if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(u, l1, l2, gamma, v);
            if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(u, l1, l2, v);
            
            // Update r
            shift = b[l*p+j] - a[j];
            if (shift !=0) {
              for (int i=0;i<n;i++) {
                si = shift*X[j*n+i];
                r[i] -= si;
                eta[i] += si;
              }
              if (fabs(shift)*sqrt(v) > maxChange) maxChange = fabs(shift)*sqrt(v);
            }
          }
        }
        
        // Check for convergence
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }
      
      // Scan for violations
      int violations = 0;
      for (int j=0; j<p; j++) {
        if (e[j]==0) {
          xwr = wcrossprod(X, r, h, n, j)/n;
          l1 = lam[l] * m[j] * alpha;
          if (fabs(xwr) > l1) {
            e[j] = 1;
            violations++;
          }
        }
      }
      if (violations==0) {
        for (int i=0; i<n; i++) {
          REAL(Eta)[l*n+i] = eta[i];
        }
        break;
      }
    }
  }
  res = cleanupCox(a, r, h, e, eta, haz, rsk, beta, Loss, iter, Eta);
  UNPROTECT(4);
  return(res);
}

