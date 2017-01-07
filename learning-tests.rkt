#lang typed/racket
(require "learning.rkt"
         "student-data-complete.rkt"
         "random-var.rkt"
         "model.rkt"
         "factors.rkt" 
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
                           (hash (set '(difficulty . "hard")) '(114 . 286) (set '(difficulty . "easy")) '(172 . 286)))
              (test-= "Checking learned params"
                      ((hash-ref (Model-cpds model) grade) (set '(grade . "a") '(difficulty . "easy") '(intelligence . "high")))
                      (exact->inexact (/ 47 52))
                      +TOLERANCE+)))

(run-tests (Discrete-Model-training))
                      
                      
                           

