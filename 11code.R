##Rcpp related code
library(Rcpp)
library(inline)

src <-
      'NumericVector X(x);
       IntegerVector Dims(dims);
       int n=Dims(0); int m=Dims(1); int k=Dims(2);
       NumericMatrix xx(n-k,m*(k+1));
       for(int i=0;i<n-k;i++) {
         for(int j=0;j<m*(k+1);j++) {
           xx(i,j)=X(m*(k+1)-j-1+i*m);
         }
       }
       return xx;'

rcppreflow <- cxxfunction(signature(x="numeric",dims="integer"),                        src,plugin="Rcpp")
