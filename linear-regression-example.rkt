#lang typed/racket
(require "random-var.rkt"
         "model.rkt"
         "learning.rkt"
         "factors.rkt"
         math/flonum
         "utils.rkt")

; Boston house prices dataset in a list of vectors.
(require "boston-data.rkt")

; first we have to define all of our variables.
(define crime (make-GaussianRandomVar 'crime))
(define ZN (make-GaussianRandomVar 'ZN))
(define INDUS (make-GaussianRandomVar 'INDUS))
(define River (make-GaussianRandomVar 'River)); '("0" "1")))
(define nox (make-GaussianRandomVar 'nox))
(define rm (make-GaussianRandomVar 'rm))
(define age (make-GaussianRandomVar 'age))
(define dis (make-GaussianRandomVar 'dis))
(define rad (make-GaussianRandomVar 'rad))
(define tax (make-GaussianRandomVar 'tax))
(define ptratio (make-GaussianRandomVar 'ptratio))
(define blacks (make-GaussianRandomVar 'blacks))
(define lstat (make-GaussianRandomVar 'lstat))
(define medv (make-GaussianRandomVar 'medv))

; and then we set their dependencies.  
(define dep-vars (list crime ZN INDUS River nox rm age dis rad tax ptratio blacks lstat))
(for ([v dep-vars])
  (add-dependency! v medv)) ; to get basic least-square regression we take the naive
                            ; assumption that every non-target var depends on the target.

; we create a model with all of our vars.  naive-model sets everything to a standard
; gaussian or uniform distribution for continuous or discrete vars, respectively.
; this is basically just a hash between each var and its conditional prob. distro.
(define model (make-naive-model (cons medv dep-vars)))

(define training-set (take boston-house-data 400))
(define testing-set (drop boston-house-data 400))


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

(printf "error\tmean prediction\tactual\n")

; okay, we're going to get some predictions now.  We're accumulating the errors here, too.
(define errors
  (for/list : (Listof Float) ([vec testing-set])

    ; these two lines are more typed/racket announces.  The first target* is a (U String Float).
    (define target* (vector-ref vec 13))
    (define target : Float (if (flonum? target*) target* 0.0))

    ; this is the probability distribution.  Bayes net inference means we try to calculate
    ; a distribution of our target var conditioned on the evidence that we can see.
    ; this is a Factor from factors.rkt.
    ; So it both contains FactorData (A hashtable containing mixtures of CanonicalFactors)
    ; and its an actual function that consumes the sets we talked about.
    (define factor (probability-conditioned-on-evidence model (set medv) (row->index-set vec)))

    ; the factor/function will return the probability of seeing a particular value.
    ; since our bayes nets are all conditional gaussians, its actually a normal distribution.
    ; we will call the mean of that distro our actual value.  To extract it, we need to pull
    ; out the canonical factor hidden in the hash table.  The hashtable contains
    ; (U Float CanonicalMixture), so we are doing more weirdness around the type checker.
    (define cf (let ([x (hash-ref (Factor-data factor)
                                  (ann (set) (Setof TableCPDIndex)))])
                 (if (list? x)
                     (cdar x)
                     (make-standard-gaussian (set medv)))))
    
    ; the formulas here are to turn the values associated to a canonical factor back into a normal.
    (define mean (fl* (fl/ 1. (unpack (CanonicalFactor-K cf))) (vector-ref (CanonicalFactor-h cf) 0)))
    
    (printf "~a\t~a\t~a\n" (abs (- mean target)) mean target)

    ; finally, we accum the square of the error.
    (fl* mean mean)))

; mean squared error.  scikit's linear_regression has 40.84.
; this should come out to 174.54.  
(displayln (flmean errors))

;(displayln ((probability-conditioned-on-evidence model (set medv) (row->index-set vec))
;            (set (cons 'medv target)))))