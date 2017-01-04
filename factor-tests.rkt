#lang typed/racket
(require "random-var.rkt"
         "factors.rkt"
         math/matrix
         math/flonum
         typed/rackunit
         typed/rackunit/text-ui
         math/distributions)


(define x (make-GaussianRandomVar 'x))
(define y (make-GaussianRandomVar 'y))
(define z (make-GaussianRandomVar 'z))
(define non (make-GaussianRandomVar 'non))
(define multivariate-standard-gaussian
  (CanonicalFactor
   (set x y z)
   (matrix [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]])
   (vector 0.0 0.0 0.0)
   (- (fllog (flsqrt (* 2.0 pi))))))

(define standard-gaussian
  (partial-application-factor multivariate-standard-gaussian
                              (set (cons 'y 0.0) (cons 'x 0.0))))

(define another-standard-gaussian
  (partial-application-factor multivariate-standard-gaussian
                              (set (cons 'y 0.0)
                                   (cons 'non 0.0)
                                   (cons 'x 0.0))))
(define Canonical-Factor-Evidence-Reduction
  (test-suite "Evidence Reduction on Canonical Factors"
              
              (test-= "Evidence reduction for multivariate standard gaussian."
                      (standard-gaussian (set (cons 'z 0.0)))
                      (flnormal-pdf 0.0 1.0 0.0 #f)
                      0.00000001)

              (test-= "Evidence reduction for evidence that includes vars outside of scope"
                      (another-standard-gaussian (set (cons 'z 0.0)))
                      (flnormal-pdf 0.0 1.0 0.0 #f)
                      0.00000001)
         
              (test-exn "Evidence reduction actually removes vars from scope."
                        exn:fail?
                        (thunk (standard-gaussian (set (cons 'x 1.0)))))))
(run-tests Canonical-Factor-Evidence-Reduction)
         
