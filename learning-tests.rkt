#lang typed/racket
(require "learning.rkt"
         "student-data-complete.rkt"
         "continuous-data.rkt"
         "continuous-data2.rkt"
         "random-var.rkt"
         "model.rkt"
         "factors.rkt"
         threading
         math/matrix
         math/flonum
         math/distributions
         typed/rackunit
         typed/rackunit/text-ui)

(define +TOLERANCE+ 0.000001)
(define (Discrete-Model-training)
  (define difficulty (make-DiscreteRandomVar 'difficulty '("hard" "easy")))
  (define grade (make-DiscreteRandomVar 'grade '("a" "b" "c")))
  (define intelligence (make-DiscreteRandomVar 'intelligence '("avg" "high")))
  (define SAT (make-DiscreteRandomVar 'SAT '("noncollege" "college")))
  (define letter (make-DiscreteRandomVar 'letter '("good" "bad")))
  (add-dependency! grade difficulty)
  (add-dependency! grade intelligence)
  (add-dependency! SAT intelligence)
  (add-dependency! letter grade)

  (define model (make-naive-model (list difficulty grade SAT letter intelligence)))
  (mle-train-model! model student-data)
  (test-suite "Training a purely discrete model with no missing data."
              (test-equal? "Creating a sufficient Statistic"
                           (create-discrete-sufficient-statistic difficulty student-data)
                           (hash (set '(difficulty . "hard")) '(114. . 286.) (set '(difficulty . "easy")) '(172. . 286.)))
              (test-= "Checking learned params"
                      ((hash-ref (Model-cpds model) grade) (set '(grade . "a") '(difficulty . "easy") '(intelligence . "high")))
                      (exact->inexact (/ 47 52))
                      +TOLERANCE+)))

(run-tests (Discrete-Model-training))

            
(define +BIG-TOLERANCE+ 0.1)
(define +BIGGER-TOLERANCE+ 0.5)

(define (Continuous-Model-learning)
  (define x (make-GaussianRandomVar 'x))
  (define z (make-DiscreteRandomVar 'z '("cat1" "cat2")))
  (define y (make-GaussianRandomVar 'y))
  (add-dependency! y x)
  (add-dependency! y z)               
  (define x-statistic (create-gaussian-sufficient-statistic x continuous-data))
  (define y-statistic (create-gaussian-sufficient-statistic y continuous-data))
  (define model1 (make-naive-model (list x y z)))
  ;(create-factor-from-data x continuous-data)
  ;((create-factor-from-data y continuous-data) (set (cons 'x 1.) (cons 'y 1.) (cons 'z "cat1")))
  ;(create-factor-from-data z continuous-data)
  (mle-train-model! model1 continuous-data)
  (define y-factor (hash-ref (Model-cpds model1) y))
  
  (define y-statistic-matrix-cat1 (~> y-statistic
                                      (hash-ref _ (set '(z . "cat1")))
                                      cdr))
  (define y-statistic-matrix-cat2 (~> y-statistic
                                      (hash-ref _ (set '(z . "cat2")))
                                      cdr))

  (define a (make-GaussianRandomVar 'a))
  (define b (make-GaussianRandomVar 'b))
  (define c (make-GaussianRandomVar 'c))
  (add-dependency! c a)
  (add-dependency! c b)

  (define model2 (make-naive-model (list a b c)))
  (mle-train-model! model2 continuous-data2)
  (define c-factor (hash-ref (Model-cpds model2) c))

  (test-suite "Training a model for conditional linear regression"
              (test-= "Checking x mean is 25 in x-statistics"
                      25.
                      (~> x-statistic
                          (hash-ref _ (ann (set) (Setof TableCPDIndex)))
                          car
                          (vector-ref _ 0))
                      +BIG-TOLERANCE+)
              (test-= "Checking x mean is 25 in y-statistics, cat1"
                      25.
                      (~> y-statistic
                          (hash-ref _ (set '(z . "cat1")))
                          car
                          (vector-ref _ 0))
                      +BIG-TOLERANCE+)
              (test-= "Checking x mean is 25 in y-statistics, cat2"
                      25.
                      (~> y-statistic
                          (hash-ref _ (set '(z . "cat2")))
                          car
                          (vector-ref _ 0))
                      +BIG-TOLERANCE+)
              (test-= "Checking x variance is 9 in x-statistics"
                      9.
                      (~> x-statistic
                          (hash-ref _ (ann (set) (Setof TableCPDIndex)))
                          cdr
                          (matrix-ref _ 0 0))
                      +BIGGER-TOLERANCE+)
              (test-= "Checking x variance is 9 in y-statistics, cat1"
                      9.
                      (~> y-statistic
                          (hash-ref _ (set '(z . "cat1")))
                          cdr
                          (matrix-ref _ 0 0))
                      +BIGGER-TOLERANCE+)
              (test-= "Checking x variance is 9 in y-statistics, cat2"
                      9.
                      (~> y-statistic
                          (hash-ref _ (set '(z . "cat2")))
                          cdr
                          (matrix-ref _ 0 0))
                      +BIGGER-TOLERANCE+)
              (test-= "Checking y-mean-slope is 4 in y-statistics, cat1"
                      4.
                      (fl* (fl/ 1. (matrix-ref y-statistic-matrix-cat1 0 0))
                           (matrix-ref y-statistic-matrix-cat1 1 0))
                      +BIGGER-TOLERANCE+)
              (test-= "Checking y-mean-slope is 2 in y-statistics, cat2"
                      2.
                      (fl* (fl/ 1. (matrix-ref y-statistic-matrix-cat2 0 0))
                           (matrix-ref y-statistic-matrix-cat2 1 0))
                      +BIGGER-TOLERANCE+)
              (test-case "Checking pdf values for y, cat1"
                         (for ([i (flvector->list (flnormal-sample 25. 3. 10))])
                           (define j (fl- (fl* 4. i) 9.))
                           (check-= (y-factor (set '(z . "cat1") (cons 'x i) (cons 'y j)))
                                    (flnormal-pdf j 2. j #f)
                                    +BIGGER-TOLERANCE+)))
              (test-case "Checking pdf values for y, cat2"
                         (for ([i (flvector->list (flnormal-sample 25. 3. 10))])
                           (define j (fl+ (fl* 2. i) 25.))
                           (check-= (y-factor (set '(z . "cat2") (cons 'x i) (cons 'y j)))
                                    (flnormal-pdf j 4. j #f)
                                    +BIGGER-TOLERANCE+)))
              (test-case "Checking pdf values for c"
                         (for ([i (flvector->list (flnormal-sample 25. 3. 10))]
                               [j (flvector->list (flnormal-sample 100. 5. 10))])
                           (define k (fl+ -5. (fl+ (fl* 0.5 i) 
                                                   (fl* 19. j))))
                                               
                           (check-= (c-factor (set (cons 'c k) (cons 'a i) (cons 'b j)))
                                    (flnormal-pdf k 4. k #f)
                                    +BIGGER-TOLERANCE+)))))

(run-tests (Continuous-Model-learning))
