! --------------------
! this function implements the CD for lasso
! algorithm for the following optimization problem
! para.in includes: epsilon, maxIter, lambda, 


      subroutine cd_lasso1(X, y, p, n, beta, paraIn, family, bInd, lam)
       
      logical :: bInd
      integer :: i, l, p, n, m, family
      double precision, dimension(n,p) :: X
      double precision, dimension(n) :: y, eta, ym, w, r, temp, Xl
      double precision, dimension(p) :: blast, beta
      double precision, dimension(3) :: paraIn
      double precision :: lam, lambda, v, z, sumEta, thresh
      lambda = lam
      epsilon = paraIn(1)
      maxIter = INT(paraIn(2))
      !print *, lambda
      epsilon = 1e-2
      
      do  i = 1, maxIter
      
           blast = beta
           
     
           DO   l = 1, p
              Xl = X(:,l)
              
                  temp= - MATMUL(X, beta)
               
                 eta= EXP(temp)
               
                
                 !print *, eta
                 ! STOP
                   !sumEta = SUM(eta)
                   !IF(sumEta == 0) THEN
                   !   bInd= .TRUE.
                   !   EXIT
                   !END IF

                  IF(family == 1) THEN
                      ym=1/eta
                      w=ym
                  END IF
                  IF(family == 2) THEN
                      ym=1/(1+eta)
                      w=ym * (1-ym)
                     
                  END IF
                  r= (y-ym)/w
                  !print *, DOT_PRODUCT(Xl, w)- SUM(Xl * w)
                  !STOP
                  !temp = DOT_PRODUCT(Xl, w)
                  !v=DOT_PRODUCT (temp, Xl)
                  !v=v/n
                  !z=DOT_PRODUCT(temp,r)/n
                  !z=z + v * beta(l)
                  v=SUM(Xl * w* Xl)
                  v=v/n
                  z=SUM(Xl *w*r)/n
                  z=z + v * beta(l)
                  thresh = MAX(dble(0), ABS(z)-lambda)
                  
                  !STOP
                  IF( l /= p) THEN
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE 
                      beta(l)=z/v
                  END IF
    
             ! END IF
        
          END DO      
          
          IF (maxIter /= 1) THEN
             IF  (bInd .EQV. .TRUE.) THEN
                  EXIT
                 ELSE IF(SUM(ABS(blast-beta))<epsilon) THEN
                  !print *, "fortran loop =         ", i, beta(1:5)
                  EXIT 
                  
              END IF
          END IF 
              
      end do

      end
