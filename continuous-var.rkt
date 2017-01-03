#lang typed/racket
(require "random-var.rkt"
         math/matrix
         math/flonum)

(struct GaussianRandomVar RandomVar ())

(define-type ContinuousIndex (Pairof Symbol Float))

(define (unpack [x : (Matrix Float)]) : Float
  (matrix-ref x 0 0))

; Canonical Factors is the space where Normal Distributions
; can live while doing factor products and factor marginalisations.
; you can do a standard gaussian like so:
;
; (define standard-gaussian
;    (CanonicalFactor (set Normal)
;                     (->row-matrix (vector 1.0))
;                     (vector 0.0)
;                     (- (fllog (flsqrt (* 2.0 pi))))))
; > (standard-gaussian (set '(normal . 0.0)))
; - : Flonum
; 0.3989422804014327
;
; compare with:
;
; (require math/distributions)
; > (flnormal-pdf 0.0 1.0 0.0 #f)
; - : Flonum
; 0.39894228040143265
(define-struct/exec CanonicalFactor
  ([scope : (Setof GaussianRandomVar)]
   [K : (Matrix Flonum)]
   [h : (Vectorof Flonum)]
   [g : Flonum])
  [(λ (self values)
     (define scope-labels
       (list->set
        (set-map (CanonicalFactor-scope self)
                 (λ ([x : GaussianRandomVar]) (RandomVar-name x)))))
     (define var-labels
       (list->set
        (set-map values (λ ([x : ContinuousIndex]) (car x)))))
     (unless (equal? scope-labels var-labels)
       (error "Factor ~a called with variables ~a which are outside of my scope: ~a"
              self values (CanonicalFactor-scope self)))
     (define X : (Matrix Float)
       (->col-matrix (set-map values (λ ([x : ContinuousIndex]) (cdr x)))))
     (flexp (+ (CanonicalFactor-g self)
               (unpack (matrix* (->row-matrix (CanonicalFactor-h self)) X))
               (* -0.5
                  (unpack (matrix* (->row-matrix X)
                                   (CanonicalFactor-K self)
                                   X))))))
    : (-> CanonicalFactor (Setof ContinuousIndex) Float)])

(define-type CanonicalMixture (Listof (Pairof Float CanonicalFactor))) 

(define-struct/exec Factor
  ([scope : (Setof RandomVar)]
   [data : (HashTable (Setof (Setof TableCPDIndex))
                      CanonicalMixture)])
  [(λ (self values)
     0.0) : (-> Factor
                (Setof (Setof (U ContinuousIndex TableCPDIndex)))                
                Float)])

