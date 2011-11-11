! --------------------
! this function implements the CD for mcp
! algorithm for the following optimization problem
! para.in includes: epsilon, maxIter, lambda, gamma


      subroutine cd_mcp1(X, y, p, n, beta, paraIn, family, bInd, lam)
       
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
      gamma = paraIn(3) 
      epsilon =1e-2
      ! print *, paraIn
      do  i = 1, maxIter

           blast = beta
     
     
           DO   l = 1, p
              Xl = X(:, l)
              
               temp = -MATMUL(X, beta)
               eta = EXP(temp)
               !   sumEta = SUM(eta)
               !    IF(sumEta == 0) THEN
               !       bInd= .TRUE.
               !       EXIT
               !   END IF

                  IF(family == 1) THEN
                      ym=1/eta
                      w=ym
                  END IF
                  IF(family == 2) THEN
                      ym=1/(1+eta)
                      w=ym * (1-ym)
                      !print *, w
                  END IF
                  
                  r=(y-ym)/w
                  v=SUM(Xl*w*Xl)
                  v=v/n
                  z=SUM(Xl*w*r)/n
                  z=z + v * beta(l)
                  thresh = MAX(dble(0), ABS(z)-lambda)
                  IF( l /= p) THEN
                     IF(ABS(z) < lambda*gamma) THEN
                        beta(l)=SIGN(dble(1),z)*thresh/(v*(1-1/gamma))
                     ELSE 
                        beta(l)=z/v
                     END IF
                  ELSE 
                      beta(l)=z/v
                  END IF
    
              !END IF
        
          END DO      
                      
       
          
          IF (maxIter /= 1) THEN
             IF  (bInd .EQV. .TRUE.) THEN
                  EXIT
                 ELSE IF(MAXVAL(ABS(blast-beta))<epsilon) THEN
                  EXIT
              END IF
          END IF 
              
      end do

      end
