#lang typed/racket
(require "random-var.rkt"
         "utils.rkt"
         math/matrix
         math/flonum)

(provide ContinuousIndex
         ContinuousIndex?
                  TableCPDIndex
         TableCPDIndex?
         TableCPD
         generate-uniform-table-cpd
         is-valid-cpd?
         (struct-out DiscreteFactor)
         ;(struct-out CanonicalFactor)
         (struct-out Factor))


(define-type ContinuousIndex (Pairof Symbol Float))
(define-predicate ContinuousIndex? ContinuousIndex)

(define-type TableCPDIndex (Pairof Symbol String))
(define-type TableCPD (HashTable (Setof TableCPDIndex)
                                 Float))
(define-predicate TableCPDIndex? TableCPDIndex)

(define (make-TableCPD-labels [vars : (Listof DiscreteRandomVar)]) : (Listof (Listof TableCPDIndex))
  (apply cartesian-product (map (λ ([var : DiscreteRandomVar])
                                  (map (λ ([label : String])
                                         (cons (RandomVar-name var) label))
                                       (DiscreteRandomVar-labels var)))
                                vars)))
(define (get-cpd-labels [var : DiscreteRandomVar]) : (Setof (Setof TableCPDIndex))
  (list->set
   (for/list : (Listof (Setof TableCPDIndex))
     ([i (make-TableCPD-labels (filter DiscreteRandomVar?
                                       (cons var (RandomVar-depends-on var))))])
     (list->set i))))

(define (get-variable-labels [var : DiscreteRandomVar]) : (Setof (Setof TableCPDIndex))
  (for/set : (Setof (Setof TableCPDIndex)) ([i (DiscreteRandomVar-labels var)])
    (set (cons (RandomVar-name var) i))))

(define (get-variable-labels* [var : DiscreteRandomVar]) : (Setof TableCPDIndex)
    (foldl (λ ([x : (Setof TableCPDIndex)]
               [y : (Setof TableCPDIndex)]) : (Setof TableCPDIndex)
             (set-union x y))
           (ann (set) (Setof TableCPDIndex))
           (set->list (get-variable-labels var))))

(define (generate-uniform-table-cpd [var : DiscreteRandomVar]) : TableCPD
  (for/hash : TableCPD ([k (make-TableCPD-labels (filter DiscreteRandomVar?
                                                         (cons var (RandomVar-depends-on var))))])
    (values (list->set k)
            (exact->inexact (/ 1 (length (DiscreteRandomVar-labels var)))))))

(define (is-valid-cpd? [var : DiscreteRandomVar]
                       [cpd : TableCPD]) : Boolean
  (define (all-cells-exists?) : Boolean
    (equal? (list->set (hash-keys cpd))
            (list->set (for/list : (Listof (Setof (Pair Symbol String)))
                         ([l (make-TableCPD-labels (filter DiscreteRandomVar?
                                                           (cons var (RandomVar-depends-on var))))])
                         (list->set l)))))
  (define (rows-sum-to-one?) : Boolean
    (define rows (make-TableCPD-labels (filter DiscreteRandomVar? (RandomVar-depends-on var))))
    (define (sum-cell-row [pre-cellkeys : (Listof TableCPDIndex)]) : Float
      (flsum (map (λ ([cellkey : (Listof TableCPDIndex)]) (hash-ref cpd (list->set cellkey)))
                  (map (λ ([col : String])
                         (cons (cons (RandomVar-name var) col) pre-cellkeys))
                       (DiscreteRandomVar-labels var)))))
    (andmap (λ ([row : (Listof TableCPDIndex)]) (= (sum-cell-row row) 1)) rows))
  (and (all-cells-exists?)
       (rows-sum-to-one?)))

(define-struct/exec DiscreteFactor
  ([scope : (Setof (Setof TableCPDIndex))]
   [func : (-> (Setof TableCPDIndex) Float)])
  [(λ (self labels)
     (unless (set-member? (DiscreteFactor-scope self) labels)
        (error "The variables ~a are outside of the scope of Factor ~a" labels self (DiscreteFactor-scope self)))
     ((DiscreteFactor-func self) labels)) : (-> DiscreteFactor (Setof TableCPDIndex) Float)])









  




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
       (error "Canonical Factor ~a called with variables ~a which are outside of my scope: ~a"
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
   [data : (HashTable (Setof TableCPDIndex)
                      CanonicalMixture)])
  [(λ (self var-values)
     (define scope-vars
       (for/hash : (HashTable Symbol RandomVar) ([v (Factor-scope self)])
         (values (RandomVar-name v) v)))
     (define discrete-values : (Setof TableCPDIndex)
       (set-filter TableCPDIndex? var-values))
     (define continuous-values : (Setof ContinuousIndex)
       (set-filter ContinuousIndex? var-values))
     (define var-labels
       (list->set
        (set-map var-values (λ ([x : (U TableCPDIndex ContinuousIndex)]) (car x)))))
     (unless (equal? var-labels (list->set (hash-keys scope-vars)))
       (error "Factor ~a called with variables ~a which are outside of my scope: ~a"
              self var-values (Factor-scope self)))
     (foldl fl*
            1.0
            (map (λ ([weighted-cfactor : (Pairof Float CanonicalFactor)])
                   (fl* (car weighted-cfactor)
                        ((cdr weighted-cfactor) continuous-values)))
                 (hash-ref (Factor-data self) discrete-values))))
     : (-> Factor
                (Setof (U ContinuousIndex TableCPDIndex))
                Float)])

