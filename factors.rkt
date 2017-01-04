#lang typed/racket
(require "random-var.rkt"
         "utils.rkt"
         math/matrix
         math/flonum
         threading)

(provide ContinuousIndex
         ContinuousIndex?
         TableCPDIndex
         TableCPDIndex?
         AnyIndex
         AnyIndex?
         TableCPD
         generate-uniform-table-cpd
         is-valid-cpd?
         (struct-out DiscreteFactor)
         (struct-out CanonicalFactor)
         partial-application-factor
         (struct-out Factor))


(define-type ContinuousIndex (Pairof Symbol Float))
(define-predicate ContinuousIndex? ContinuousIndex)

(define-type TableCPDIndex (Pairof Symbol String))
(define-type TableCPD (HashTable (Setof TableCPDIndex)
                                 Float))
(define-predicate TableCPDIndex? TableCPDIndex)

(define-type AnyIndex (U ContinuousIndex TableCPDIndex))
(define-predicate AnyIndex? AnyIndex)

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



(: canonical-order (case-> (-> (Setof GaussianRandomVar) (Listof GaussianRandomVar))
                           (-> (Setof ContinuousIndex) (Listof ContinuousIndex))
                           (-> (Setof Symbol) (Listof Symbol))))
(define (canonical-order items)
  (define-predicate set-var? (Setof GaussianRandomVar))
  (define-predicate set-index? (Setof ContinuousIndex))
  (define-predicate set-sym? (Setof Symbol))
  (sort (set->list items)
        (cond [(set-var? items)
               (λ ([ x : GaussianRandomVar ]
                   [ y : GaussianRandomVar ])
                 (symbol<? (RandomVar-name x) (RandomVar-name y)))]
              [(set-index? items)
               (λ ([ x : ContinuousIndex ]
                   [ y : ContinuousIndex ])
                 (symbol<? (car x) (car y)))]
              [(set-sym? items) symbol<?])))

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
;
; we order the variable-set used in ((CanonicalFactor ...) variable-set)
; in alphabetical order by the RandomVar-name.  We asusme that the Matrix K
; and Vector h obey this convention.
(define-struct/exec CanonicalFactor
  ([scope : (Setof GaussianRandomVar)]
   [K : (Matrix Flonum)]
   [h : (Vectorof Flonum)]
   [g : Flonum])
  [(λ (self values)
     (define scope-labels
       (get-var-labels (CanonicalFactor-scope self)))
       
     (define var-labels
       (list->set
        (set-map values (λ ([x : ContinuousIndex]) (car x)))))
     (unless (equal? scope-labels var-labels)
       (error "Canonical Factor ~a called with variables ~a which are outside of my scope: ~a"
              self values (CanonicalFactor-scope self)))
     (define X : (Matrix Float)
       (->col-matrix (map (λ ([x : ContinuousIndex]) (cdr x))
                          (canonical-order values))))
     (flexp (+ (CanonicalFactor-g self)
               (unpack (matrix* (->row-matrix (CanonicalFactor-h self)) X))
               (* -0.5
                  (unpack (matrix* (->row-matrix X)
                                   (CanonicalFactor-K self)
                                   X))))))
   : (-> CanonicalFactor (Setof ContinuousIndex) Float)])

(define-type CanonicalMixture (Listof (Pairof Float CanonicalFactor)))
(define-type FactorData (HashTable (Setof TableCPDIndex)
                                   (U Float CanonicalMixture)))

(define-struct/exec Factor
  ([scope : (Setof RandomVar)]
   [data : FactorData]
   [func : (U Null (-> (Setof AnyIndex) Float))])
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
     (if (null? (Factor-func self))
         (let ([factor (hash-ref (Factor-data self) discrete-values)])
           (cond [(list? factor)
                  (foldl fl*
                         1.0
                         (map (λ ([weighted-cfactor : (Pairof Float CanonicalFactor)])
                                (fl* (car weighted-cfactor)
                                     ((cdr weighted-cfactor) continuous-values)))
                              factor))]
                 [(flonum? factor)
                  factor]))
         ((Factor-func self) var-values)))
   : (-> Factor
         (Setof AnyIndex)
         Float)])


(: partial-application-factor (case-> (-> CanonicalFactor
                                          (Setof ContinuousIndex)
                                          CanonicalFactor)
                                      (-> Factor (Setof AnyIndex) Factor)))
(define (partial-application-factor factor evidence)
  (define evidence-labels (set-map evidence (λ ([x : AnyIndex]) (car x))))
  (cond [(CanonicalFactor? factor)
         (define reduced-scope
           (set-filter (λ ([v : GaussianRandomVar]) (false? (member (RandomVar-name v) evidence-labels)))
                       (CanonicalFactor-scope factor)))
         (define scope-index (order-of-elems (canonical-order reduced-scope)
                                             (canonical-order (CanonicalFactor-scope factor))))
         (define evidence-names (~> evidence
                                    set->list
                                    (map (λ ([ x : ContinuousIndex ]) (car x)) _)
                                    list->set
                                    canonical-order))
         (define evidence-index         
           (order-of-elems (~> (set-subtract (CanonicalFactor-scope factor)
                                             reduced-scope)
                               (set-map _ (λ ([v : RandomVar]) (RandomVar-name v)))
                               list->set
                               canonical-order)
                           evidence-names))
                                                
         (define K_scope (submatrix (CanonicalFactor-K factor)
                                    scope-index
                                    scope-index))
         (define K_evidence (submatrix (CanonicalFactor-K factor)
                                       evidence-index
                                       evidence-index))
         (define K_scope/evidence (submatrix (CanonicalFactor-K factor)
                                             scope-index
                                             evidence-index))
         (define h (->col-matrix (CanonicalFactor-h factor)))
         (define h_scope (~> h
                             (submatrix _ scope-index (list 0))
                             ->col-matrix))
           
         (define h_evidence (~> h
                             (submatrix _ evidence-index (list 0))
                             ->col-matrix))
         (define evidence-vector (~> evidence
                                     canonical-order
                                     (map (λ ([x : ContinuousIndex]) (cdr x)) _)
                                     ->col-matrix
                                     (submatrix _ evidence-index (list 0))
                                     ->col-matrix))
         (CanonicalFactor reduced-scope
                          K_scope
                          (matrix->vector (matrix- h_scope
                                                   (matrix* K_scope/evidence evidence-vector)))
                          (+ (CanonicalFactor-g factor)
                             (unpack (matrix* (matrix-transpose h_evidence) evidence-vector))
                             (fl* -0.5 (unpack (matrix* (matrix-transpose evidence-vector) K_evidence evidence-vector)))))]
                                    
                                    
        [(Factor? factor) factor]))
  
                                    

(define (partial-application-factor-old [factor : Factor]
                                    [evidence : (Setof AnyIndex)]) : Factor
  (define evidence-labels (set-map evidence (λ ([x : AnyIndex]) (car x))))
  (define reduced-scope
    (set-filter (λ ([v : RandomVar]) (false? (member (RandomVar-name v) evidence-labels)))
                (Factor-scope factor)))
  (define discrete-scope : (Setof DiscreteRandomVar)
    (set-filter DiscreteRandomVar? reduced-scope))
  (define continuous-scope : (Setof GaussianRandomVar)
    (set-filter GaussianRandomVar? reduced-scope))
  (define reduced-evidence
    (set-filter (λ ([x : AnyIndex])
                  (set-member? (get-var-labels (Factor-scope factor)) (car x)))
                evidence))
  (Factor reduced-scope
          (ann (make-hash) FactorData)
          (λ ([var-values : (Setof AnyIndex)]) : Float
            (factor (set-union var-values reduced-evidence)))))

(define (labels-within-factor [f : Factor]
                              [labels : (Setof AnyIndex)]) : (Setof AnyIndex)
  (define (labels-within-scope [symbols : (Setof Symbol)]
                               [labels : (Setof AnyIndex)]) : (Setof AnyIndex)
     (set-filter (λ ([x : AnyIndex]) (set-member? symbols (car x))) labels))
  (define (get-factor-symbols) : (Setof Symbol)
    (list->set (set-map (Factor-scope f) RandomVar-name)))
  (labels-within-scope (get-factor-symbols) labels))
  

(: product-factor (-> Factor * Factor))
(define (product-factor . factors)
  (define (helper [f1 : Factor]
                  [f2 : Factor]) : Factor
    (Factor (set-union (Factor-scope f1) (Factor-scope f2))
            (ann (make-hash) FactorData)
            (λ ([var-values : (Setof AnyIndex)]) : Float
              (fl* (f1 (labels-within-factor f1 var-values))
                   (f2 (labels-within-factor f2 var-values))))))
  (foldl helper (car factors) (cdr factors)))
  
(define (factor-marginalization [factor : Factor]
                                [var : RandomVar]) : Factor
  (unless (set-member? (Factor-scope factor) var)
    (error "Variable ~a is not in the scope of Factor ~a" var factor))
  ;(cond [(DiscreteRandomVar? var)
  ;       (Factor (set-subtract (Factor-scope var))
  ;               (ann (make-hash) FactorData)
  ;               (λ ([var-values : (Setof AnyIndex)]) : Float
                   
                      
  
  factor)

                                

