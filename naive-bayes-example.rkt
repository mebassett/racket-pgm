#lang typed/racket
(require "random-var.rkt"
         "model.rkt"
         "learning.rkt"
         "factors.rkt"
         math/flonum
         threading
         math/flonum
         "utils.rkt")

; iris in a list of vectors.
(require "car-data.rkt")

; defining our variables
(define price (make-DiscreteRandomVar 'price '("vhigh" "high" "med" "low")))
(define maint (make-DiscreteRandomVar 'maint '("vhigh" "high" "med" "low")))
(define doors (make-DiscreteRandomVar 'doors '("2" "3" "4" "5more")))
(define persons (make-DiscreteRandomVar 'persons '("2" "3" "more")))
(define boot (make-DiscreteRandomVar 'boot '("small" "med" "big")))
(define safety (make-DiscreteRandomVar 'safety '("low" "med" "high")))
(define class (make-DiscreteRandomVar 'class '("unacc" "acc" "good" "vgood")))

; and their dependencies.
(define dep-vars (list price maint doors persons boot safety))
(for ([v dep-vars])
  (add-dependency! v class))

; we create a model with all of our vars.  naive-model sets everything to a standard
; gaussian or uniform distribution for continuous or discrete vars, respectively.
; this is basically just a hash between each var and its conditional prob. distro.
(define model (make-naive-model (cons class dep-vars)))
(displayln (Factor-data (hash-ref (Model-cpds model) class)))

(define training-set (take car-data 1500))
(define testing-set (drop car-data 1500))

; parameter learning/training for our model
(mle-train-model! model training-set)

; this model deals only in probability distributions.
; eventually we will have "factors" which correspond to some
; distribution.  these factors act on sets of pairs.  The car
; corresponds to the var name and the cdr the actual value.
;
; this is a helper function that turns our vectors in the dataset
; to a set that can be passed to a factor.
(define (row->index-set [row : Row]) : (Setof AnyIndex)
  (for/set : (Setof AnyIndex) ([var dep-vars]
                               [val row])
    ; more typed/racket inference annoyances.  it does not know that
    ; (Pairof Symbol (U Float String)) = (U (Pairof Symbol Float) (Pairof Symbol String))
    ; so this doesn't work.
    ; (cons (RandomVar-name dep-vars) val)
    (cond [(string? val)
           (cons (RandomVar-name var) val)]
          [(flonum? val)
           (cons (RandomVar-name var) val)])))

(printf "actual\tprob\tpredicting\n\n")

(define sum-correct
(for/sum : Real ([vec testing-set])

  ; these two lines are more typed/racket announces.  The first target* is a (U String Float).
  (define target* (vector-ref vec 6))
  (define target : String (if (string? target*) target* "no class"))
  
  ; this is the probability distribution.  Bayes net inference means we try to calculate
  ; a distribution of our target var conditioned on the evidence that we can see.
  ; this is a Factor from factors.rkt.
  ; So it both contains FactorData (A hashtable containing mixtures of CanonicalFactors)
  ; and its an actual function that consumes the sets we talked about.

  (define factor (probability-conditioned-on-evidence model (set class) (row->index-set vec)))

  ; this code searches the Factor-data for the index with the highest prob, which we take
  ; to be the predicting value.
  (define predicted (~> (for/list : (Listof (Pairof (Setof TableCPDIndex) Float)) ([(k v) (Factor-data factor)])
                          (cons k (if (flonum? v) v 0.)))
                        (sort _  (λ ([x : (Pairof (Setof TableCPDIndex) Float)]
                                     [y : (Pairof (Setof TableCPDIndex) Float)])
                                   (< (cdr x) (cdr y))))
                        last
                        car
                        set->list
                        cdar))
  (printf "~a\t~a\t~a\n" target (factor (set (cons 'class target))) predicted)
  (if (equal? predicted target) 1. 0.0001)))

(printf "\nAccuracy : ~a" (/ sum-correct (exact->inexact (length testing-set))))

  
