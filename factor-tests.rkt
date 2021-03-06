#lang typed/racket
(require "random-var.rkt"
         "factors.rkt"
         math/matrix
         math/flonum
         typed/rackunit
         typed/rackunit/text-ui
         math/distributions)


(define +TOLERANCE+ 0.0005)

(define x (make-GaussianRandomVar 'x))
(define y (make-GaussianRandomVar 'y))
(define z (make-GaussianRandomVar 'z))
(define non (make-GaussianRandomVar 'non))
(define multivariate-standard-gaussian (make-standard-gaussian (set x y z)))

(define standard-gaussian
  (partial-application-factor multivariate-standard-gaussian
                              (set (cons 'y 0.0) (cons 'x 0.0))))
 
(define another-standard-gaussian
  (partial-application-factor multivariate-standard-gaussian
                              (set (cons 'y 0.0)
                                   (cons 'non 0.0)
                                   (cons 'x 0.0))))
(define f (make-gaussian (set x y) (vector 1.0 5.0) (matrix [[0.1 0.0] [0.0 0.5]])))

(define Canonical-Factor-Evidence-Reduction
  (test-suite "Evidence Reduction on Canonical Factors"
              
              (test-= "Evidence reduction for multivariate standard gaussian."
                      (standard-gaussian (set (cons 'z 0.0)))
                      (multivariate-standard-gaussian (set '(z . 0.0) '(y . 0.0) '(x . 0.0)))
                      +TOLERANCE+)

              (test-= "Evidence reduction for evidence that includes vars outside of scope"
                      (another-standard-gaussian (set (cons 'z 0.0)))
                      (multivariate-standard-gaussian (set '(z . 0.0) '(y . 0.0) '(x . 0.0)))
                      +TOLERANCE+)
              (test-= "Evidence reduction for non-standard gaussions."
                      (f (set '(x . 1.0) '(y . 5.0)))
                      ((partial-application-factor f (set '(y . 5.0))) (set '(x . 1.0)))
                      +TOLERANCE+)

         
              (test-exn "Evidence reduction actually removes vars from scope."
                        exn:fail?
                        (thunk (standard-gaussian (set (cons 'x 1.0)))))))


(define difficulty (make-DiscreteRandomVar 'difficulty '("hard" "easy")))
(define grade (make-DiscreteRandomVar 'grade '("a" "b" "f")))
(define intelligence (make-DiscreteRandomVar 'iq '("high" "low")))
(define SAT (make-DiscreteRandomVar 'SAT '("high" "low")))
(add-dependency! SAT intelligence)
(add-dependency! grade difficulty)
(add-dependency! grade intelligence)
(define grade-hash
  #hash(
        (((difficulty . "easy") (iq . "low") (grade . "a")) . 0.3)
        (((grade . "b") (difficulty . "easy") (iq . "low")) . 0.4)
        (((grade . "f") (difficulty . "easy") (iq . "low")) . 0.3)
        (((difficulty . "hard") (iq . "low") (grade . "a")) . 0.05)
        (((difficulty . "hard") (grade . "b") (iq . "low")) . 0.25)
        (((difficulty . "hard") (grade . "f") (iq . "low")) . 0.7)
        (((iq . "high") (difficulty . "easy") (grade . "a")) . 0.9)
        (((iq . "high") (grade . "b") (difficulty . "easy")) . 0.08)
        (((iq . "high") (grade . "f") (difficulty . "easy")) . 0.02)
        (((iq . "high") (difficulty . "hard") (grade . "a")) . 0.5)
        (((iq . "high") (difficulty . "hard") (grade . "b")) . 0.3)
        (((iq . "high") (difficulty . "hard") (grade . "f")) . 0.2)
        ))

(define sat-hash
  #hash(
        (((SAT . "low") (iq . "low")) . 0.95)
        (((SAT . "high") (iq . "low")) . 0.05)
        (((iq . "high") (SAT . "low")) . 0.2)
        (((iq . "high") (SAT . "high")) . 0.8)
        ))

(define grade-factor
  (Factor (set grade difficulty intelligence)
          (for/hash : FactorData ([(k v) grade-hash])
            (values (list->set k) v))
          null))
(define sat-factor
  (Factor (set SAT intelligence)
          (for/hash : FactorData ([(k v) sat-hash])
            (values (list->set k) v))
          null))
  
(define grade-iq-high (partial-application-factor grade-factor (set '(iq . "high"))))         

(define Discrete-Factor-Evidence-Reduction
  (test-suite "Evidence Reduction on Purely Discrete Factors"
              
              (test-= "Evidence reduction for discrete variable."
                      (grade-iq-high (set '(difficulty . "hard") '(grade . "f")))
                      (grade-factor (set '(iq . "high") '(difficulty . "hard") '(grade . "f")))
                      0.0)
         
              (test-exn "Evidence reduction actually removes vars from scope."
                        exn:fail?
                        (thunk (grade-iq-high (set '(iq . "high") '(difficulty . "hard") '(grade . "f")))))))

(define continuous-grade (make-GaussianRandomVar 'grade))
(define continuous-SAT (make-GaussianRandomVar 'SAT))
(define health (make-GaussianRandomVar 'health))
(define continuous-intelligence (make-GaussianRandomVar 'iq))

(add-dependency! continuous-grade health)
(add-dependency! continuous-grade continuous-intelligence)
(add-dependency! continuous-grade difficulty)

(add-dependency! continuous-SAT continuous-intelligence)

(define continuous-SAT-factor
  (Factor (set continuous-SAT continuous-intelligence)
          (hash (ann (set) (Setof TableCPDIndex))
                (list (cons 1.0
                            (make-standard-gaussian (set continuous-intelligence continuous-SAT)))))
          null))

(define continuous-grade-factor
  (Factor (set continuous-grade continuous-intelligence health difficulty)
          (hash (set '(difficulty . "easy"))
                (list (cons 1.0 (make-gaussian (set continuous-grade health continuous-intelligence)
                                               (vector 85.0 0.95 100.0)
                                               (identity-matrix 3 1.0 0.0))))
                (set '(difficulty . "hard"))
                (list (cons 1.0 (make-gaussian (set continuous-grade health continuous-intelligence)
                                               (vector 65.0 0.85 115.0)
                                               (identity-matrix 3 1.0 0.0)))))
          null))
                
(define cgrade-when-hard
  (partial-application-factor continuous-grade-factor
                              (set '(difficulty . "hard"))))
(define cgrade-when-low-health
  (partial-application-factor continuous-grade-factor
                              (set '(health . 0.5))))

(define cgrade-when-easy-and-intelligent
  (partial-application-factor continuous-grade-factor
                              (set '(difficulty . "easy")
                                   '(iq . 130.0))))

(define (Factor-Evidence-Reduction)
  (define crime (make-GaussianRandomVar 'crime))
  (define crime-factor (Factor (set crime)
                               (hash (ann (set) (Setof TableCPDIndex))
                                     (list
                                      (cons 1.
                                            (make-gaussian (set crime)
                                                           (vector 3.59)
                                                           (matrix [[73.76]])))))
                               null))
  (test-suite "Evidence Reduction on general factors"
              (test-= "Reducing a continuos variable to itself"
                      (crime-factor (set '(crime . 3.59)))
                      ((partial-application-factor crime-factor (set '(crime . 3.59))) (ann (set) (Setof AnyIndex)))
                      +TOLERANCE+)
              (test-= "Reducing on discrete variable."
                      (cgrade-when-hard (set '(grade . 65.0) '(health . 0.85) '(iq . 115.0)))
                      (multivariate-standard-gaussian (set '(x . .0) '(y . .0) '(z . .0)))
                      +TOLERANCE+)
              (test-= "Reducing on continuous variable"
                      (cgrade-when-low-health (set '(iq . 110.0) '(difficulty . "hard") '(grade . 96.0)))
                      (continuous-grade-factor (set '(health . 0.5) '(iq . 110.0) '(difficulty . "hard") '(grade . 96.0)))
                      +TOLERANCE+)
              (test-= "Reducing on a mix of continuous and discrete variables"
                      (cgrade-when-easy-and-intelligent (set '(grade . 96.0) '(health . 0.99)))
                      (continuous-grade-factor (set '(health . 0.99) '(iq . 130.0) '(difficulty . "easy") '(grade . 96.0)))
                      +TOLERANCE+)))

(define (Canonical-Factor-Product)
  (define f1 (CanonicalFactor (set x y) (matrix [[1. 2.] [0. 3.]]) (vector 0. 0.) 0.))
  (define f2 (CanonicalFactor (set y z) (matrix [[0. 0.] [1. 0.]]) (vector 0. 0.) 0.))
  (define f3 (CanonicalFactor (set x z) (matrix [[-1. -1.] [-1. -1.]]) (vector 0. 0.) 0.))
  (define f4 (make-standard-gaussian (set x)))
  (define f5 (make-standard-gaussian (set y)))
  (define f6 (make-standard-gaussian (set x y z)))
  
  (define f7 (CanonicalFactor (set x z) (matrix [[-1. -1.] [-1. -1.]]) (vector 3. 3.) 0.))
  (define f8 (CanonicalFactor (set x y) (matrix [[-1. -1.] [-1. -1.]]) (vector 2. 2.) 0.))
  (test-suite "Products on Canonical Factors"
              ;(test-check "Matrix join 1" matrix= (sum-joint-K f1 f3) (matrix [[0. 2. -1.] [0. 3. 0.] [-1. 0. -1.]]))
              ;(test-check "Matrix join 2" matrix= (sum-joint-K f1 f2) (matrix [[1. 2. 0.] [0. 3. 0.] [0. 1. 0.]]))
              (test-= "Products on Canonical Factors"
                      ((product-factor f4 f5) (set '(x . 0.0) '(y . 0.)))
                      (fl* (f4 (set '(x . .0))) (f5 (set '(y . 0.))))
                      +TOLERANCE+)
              (test-true "K-Matrix subtraction on Canonical Factors"
                         (matrix= (K- f6 f5)
                                  (matrix [[1. 0. 0.] [0. 0. 0.] [0. 0. 1.]])))
              (test-equal? "h-vector subtraction on Canonical Factors"
                           (h- f7 f8)
                           (vector 1. -2. 3.))))

(define Factor-Product
  (test-suite "Products on Factors"
              (test-= "Product on two discrete factors"
                      ((product-factor grade-factor sat-factor) (set '(SAT . "low") '(iq . "high") '(difficulty . "hard") '(grade . "f")))
                      (fl* (sat-factor (set '(SAT . "low") '(iq . "high"))) (grade-factor (set '(iq . "high") '(difficulty . "hard") '(grade . "f"))))
                      +TOLERANCE+)
              (test-= "Product on a discrete factor and a mixed/continuous factor"
                      ((product-factor continuous-SAT-factor cgrade-when-low-health)
                       (set '(grade . 65.0) '(SAT . 0.0)
                            '(iq . 0.0) '(difficulty . "hard")))
                      (fl* (cgrade-when-low-health (set '(grade . 65.0) '(iq . 0.0) '(difficulty . "hard")))
                           (continuous-SAT-factor (set '(iq . 0.0) '(SAT . 0.0))))
                      +TOLERANCE+)))

(define (Factor-Marginalization)
  (define f (make-standard-gaussian (set x y)))
  (define diff (make-DiscreteRandomVar 'diff (list "hard" "easy")))
  (define iq (make-GaussianRandomVar 'iq))
  (define cgrade (make-GaussianRandomVar 'grade))
  (add-dependency! cgrade diff)
  (add-dependency! cgrade iq)
  (define f-grade
    (Factor (set diff cgrade iq)
            (hash (set '(diff . "hard"))
                  (list (cons 1. (make-gaussian (set cgrade iq) (vector 80. 110.) (matrix [[10. 0.] [0. 16.]]))))
                  (set '(diff . "easy"))
                  (list (cons 1. (make-gaussian (set cgrade iq) (vector 92. 100.) (matrix [[8. 0.] [0. 25.]])))))
            null))
  (define f-grade-no-iq
    (Factor (set diff cgrade)
            (hash (set '(diff . "hard"))
                  (list (cons 1. (make-gaussian (set cgrade) (vector 80.) (matrix [[10.]]))))
                  (set '(diff . "easy"))
                  (list (cons 1. (make-gaussian (set cgrade) (vector 92.) (matrix [[8.]])))))
            null))
  (define f-grade-no-iq* (factor-marginalization f-grade iq))
  (define f-grade-no-iq-diff (factor-marginalization f-grade-no-iq diff))

  (define grade (make-DiscreteRandomVar 'grade '("a" "b" "f")))
  (define letter (make-DiscreteRandomVar 'letter '("good" "poor")))
  (add-dependency! letter grade)
  (define letter-factor (make-factor-from-table* letter
                                                 #hash(
                                                       (((letter . "poor") (grade . "a")) . 0.1)
                                                       (((letter . "good") (grade . "a")) . 0.9)
                                                       (((grade . "b") (letter . "poor")) . 0.4)
                                                       (((grade . "b") (letter . "good")) . 0.6)
                                                       (((grade . "f") (letter . "poor")) . 0.99)
                                                       (((grade . "f") (letter . "good")) . 0.01)
                                                       )))
  (define letter-no-grade (factor-marginalization letter-factor grade))
  
   (test-suite "Marginalizing a single variable from a factor"
               (test-= "Marginalizing a continuous random var across a CanonicalForm"
                       ((canonical-factor-marginalization f y) (set '(x . 0.)))
                       (flnormal-pdf 0. 1.0 0. #f)
                       +TOLERANCE+)
               (test-= "Marginalizing a continuous random var across a Factor"
                       (f-grade-no-iq* (set '(diff . "hard") '(grade . 80.)))
                       (f-grade-no-iq (set '(diff . "hard") '(grade . 80.)))
                       +TOLERANCE+)
               (test-= "Marginalizing a discrete random var across a Factor"
                       (f-grade-no-iq-diff (set '(grade . 80.)))
                       (fl+ (f-grade-no-iq (set '(diff . "easy") '(grade . 80.)))
                            (f-grade-no-iq (set '(diff . "hard") '(grade . 80.))))
                       +TOLERANCE+)))
                           
                           

  (run-tests Canonical-Factor-Evidence-Reduction)
  (run-tests Discrete-Factor-Evidence-Reduction)
  (run-tests (Factor-Evidence-Reduction))
  (run-tests (Canonical-Factor-Product))
  (run-tests Factor-Product)
  (run-tests (Factor-Marginalization))

  