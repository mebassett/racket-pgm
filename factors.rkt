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
         make-standard-gaussian
         make-gaussian 
         partial-application-factor
         product-factor
         canonical-factor-marginalization
         factor-marginalization
         FactorData
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

(define (make-standard-gaussian [vars : (Setof GaussianRandomVar)]) : CanonicalFactor
  (CanonicalFactor vars
                   (identity-matrix (set-count vars) 1.0 0.0)
                   (list->vector (set-map vars (λ (x) 0.0)))
                   (- (fllog (flexpt (* 2.0 pi) (/ (->fl (set-count vars)) 2))))))

(define (make-gaussian [vars : (Setof GaussianRandomVar)]
                       [means : (Vectorof Float)]
                       [Σ : (Matrix Float)]) : CanonicalFactor
  (define K (matrix-inverse Σ))
  (define μ (->col-matrix means))
  (define h (matrix* K μ))
  (define g (- (- (* 0.5 (unpack (matrix* (matrix-transpose μ) K μ))))
               (fllog (* (flsqrt (matrix-determinant Σ))
                         (flexpt (* 2.0 pi) (/ (->fl (set-count vars)) 2.0))))))
  (CanonicalFactor vars K (matrix->vector h) g))
                  

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
         (define reduced-evidence-labels (filter (λ ([label : Symbol])
                                                   (member label (set-map (CanonicalFactor-scope factor) RandomVar-name)))
                                                 evidence-labels))
         (define reduced-evidence-index
           (let* ([var-labels (list->set (set-map (CanonicalFactor-scope factor) RandomVar-name))]
                  [reduced-evidence-labels (list->set
                                            (filter (λ ([label : Symbol])
                                                      (set-member? var-labels label))
                                                    evidence-labels))])
             (order-of-elems (canonical-order reduced-evidence-labels)
                             (canonical-order var-labels))))
                                             
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
                                       reduced-evidence-index
                                       reduced-evidence-index))
         (define K_scope/evidence (submatrix (CanonicalFactor-K factor)
                                             scope-index
                                             reduced-evidence-index))
         (define h (->col-matrix (CanonicalFactor-h factor)))
         (define h_scope (~> h
                             (submatrix _ scope-index (list 0))
                             ->col-matrix))
           
         (define h_evidence (~> h
                                (submatrix _ reduced-evidence-index (list 0))
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
                                    
        [(Factor? factor)
         (define reduced-scope
           (set-filter (λ ([v : RandomVar]) (false? (member (RandomVar-name v) evidence-labels)))
                       (Factor-scope factor)))
         (define discrete-scope : (Setof DiscreteRandomVar)
           (set-filter DiscreteRandomVar? reduced-scope))
         (define continuous-scope : (Setof GaussianRandomVar)
           (set-filter GaussianRandomVar? reduced-scope))
         (define reduced-evidence
           (set-filter (λ ([evi : AnyIndex])
                         (set-member? (list->set (set-map (Factor-scope factor)
                                                          RandomVar-name))
                                      (car evi)))
                       evidence))
         (define discrete-evidence : (Setof TableCPDIndex)
           (set-filter TableCPDIndex? reduced-evidence))
         (define continuous-evidence : (Setof ContinuousIndex)
           (set-filter ContinuousIndex? reduced-evidence))
         (define new-data
           (for/hash : FactorData
             ([list-labels (make-TableCPD-labels (set->list discrete-scope))])
             (let ([datum : (U Float CanonicalMixture)
                          (hash-ref (Factor-data factor)
                                    (set-union (list->set list-labels)
                                               discrete-evidence))])
               (values (list->set list-labels)
                       (cond [(flonum? datum) datum]
                             [(list? datum)
                              (if (> (set-count continuous-evidence) 0)
                                  (map (λ ([ x : (Pairof Float CanonicalFactor)])
                                         (cons (car x) (partial-application-factor (cdr x) continuous-evidence)))
                                       datum)
                                  datum)])))))
         (Factor reduced-scope new-data null)]))
             

(define (labels-within-factor [f : Factor]
                              [labels : (Setof AnyIndex)]) : (Setof AnyIndex)
  (define (labels-within-scope [symbols : (Setof Symbol)]
                               [labels : (Setof AnyIndex)]) : (Setof AnyIndex)
    (set-filter (λ ([x : AnyIndex]) (set-member? symbols (car x))) labels))
  (define (get-factor-symbols) : (Setof Symbol)
    (list->set (set-map (Factor-scope f) RandomVar-name)))
  (labels-within-scope (get-factor-symbols) labels))

(: get-joint-scope (case-> (-> (Listof CanonicalFactor)  (Setof GaussianRandomVar))
                           (-> (Listof Factor)  (Setof RandomVar))))
(define (get-joint-scope factors)
  (cond [(CanonicalFactor? (car factors))
         (apply (curry (inst set-union GaussianRandomVar) 
                       (ann (set) (Setof GaussianRandomVar)))
                (map CanonicalFactor-scope factors))]
        [(Factor? (car factors))
         (apply (curry (inst set-union RandomVar) 
                       (ann (set) (Setof RandomVar)))
                (map Factor-scope factors))]))
         

(: joint-var-order (-> CanonicalFactor * (HashTable Integer GaussianRandomVar))) 
(define (joint-var-order . factors)
  (define joint-scope (get-joint-scope factors))
  (for/hash : (HashTable Integer GaussianRandomVar)
    ([i (in-range (set-count joint-scope))]
     [var joint-scope])
    (values i var)))

(define (get-var-index [v : GaussianRandomVar]
                       [f : CanonicalFactor]) : Integer
  (unless (set-member? (CanonicalFactor-scope f) v)
    (error "RandomVar ~a is not in scope ~a of Factor ~a" v (CanonicalFactor-scope f) f))
  (car (order-of-elems (list v) (canonical-order (CanonicalFactor-scope f)))))

(: sum-joint-K (-> CanonicalFactor * (Matrix Float)))
(define (sum-joint-K . factors)
  (define joint-scope (get-joint-scope factors))
  (define index-var-lookup : (HashTable Integer GaussianRandomVar)
    (apply joint-var-order factors))
  (build-matrix (set-count joint-scope)
                (set-count joint-scope)
                (λ ([r : Integer] [c : Integer]) : Float
                  (define r-var (hash-ref index-var-lookup r))
                  (define c-var (hash-ref index-var-lookup c))
                  (define elem-factors
                    (set-filter (λ ([f : CanonicalFactor])
                                  (and (set-member? (CanonicalFactor-scope f) r-var)
                                       (set-member? (CanonicalFactor-scope f) c-var)))
                                (list->set factors)))
                  (flsum (set-map elem-factors (λ ([f : CanonicalFactor])
                                                 (matrix-ref (CanonicalFactor-K f)
                                                             (get-var-index r-var f)
                                                             (get-var-index c-var f))))))))

(: sum-joint-h (-> CanonicalFactor * (Vectorof Float)))
(define (sum-joint-h . factors)
  (define joint-scope (get-joint-scope factors))
  (define index-var-lookup : (HashTable Integer GaussianRandomVar)
    (apply joint-var-order factors))
  (for/vector : (Vectorof Float)
    ([var joint-scope])
    (define elem-factors
      (set-filter (λ ([f : CanonicalFactor])
                    (set-member? (CanonicalFactor-scope f) var))
                  (list->set factors)))
    (flsum (set-map elem-factors (λ ([f : CanonicalFactor])
                                   (vector-ref (CanonicalFactor-h f)
                                               (get-var-index var f)))))))

(define (mixed-product [ fl : Float ]
                       [ mix : CanonicalMixture ]) : CanonicalMixture
  (map (λ ([ cm : (Pairof Float CanonicalFactor)])
           (let ([f (cdr cm)])
             (cons (car cm)
                   (CanonicalFactor (CanonicalFactor-scope f)
                                    (CanonicalFactor-K f)
                                    (CanonicalFactor-h f)
                                    (fl+ (fllog fl) (CanonicalFactor-g f))))))
       mix))
  
(define (canonical-mixture-product [cm1 : CanonicalMixture]
                                   [cm2 : CanonicalMixture]) : CanonicalMixture
  (define normalization-factor
    (flsum (for*/list : (Listof Float) ([pair1 cm1]
                                        [pair2 cm2])
             (fl* (car pair1) (car pair2)))))
  (for*/list : CanonicalMixture ([pair1 cm1]
                                 [pair2 cm2])
    (cons (fl/ (fl* (car pair1) (car pair2)) normalization-factor)
          (product-factor (cdr pair1) (cdr pair2)))))    
  

(: product-factor (case-> (-> CanonicalFactor * CanonicalFactor)
                          (-> Factor * Factor)))
(define (product-factor . factors)
  (cond [(CanonicalFactor? (car factors))
         (define joint-scope (get-joint-scope factors))
         (CanonicalFactor joint-scope
                          (apply sum-joint-K factors)
                          (apply sum-joint-h factors)
                          (flsum (map CanonicalFactor-g factors)))]
        [(Factor? (car factors))
         (define (single-product [f1 : Factor]
                                 [f2 : Factor]) : Factor
           (define joint-scope (get-joint-scope factors))
           (define discrete-scope : (Setof DiscreteRandomVar) (set-filter DiscreteRandomVar? joint-scope))
           (define continuous-scope : (Setof GaussianRandomVar) (set-filter GaussianRandomVar? joint-scope))
           (define new-data
             (for/hash : FactorData ([label-list (make-TableCPD-labels (set->list discrete-scope))])
               (let* ([label (list->set label-list)]
                      [data1 : FactorData (Factor-data f1)]
                      [data2 : FactorData (Factor-data f2)]
                      [label1 : (Setof TableCPDIndex) (set-filter TableCPDIndex?
                                                                  (labels-within-factor f1 label))]
                      [label2 : (Setof TableCPDIndex) (set-filter TableCPDIndex?
                                                                  (labels-within-factor f2 label))]
                      [val1 : (U Float CanonicalMixture) (hash-ref data1
                                                                   label1)]
                      [val2 : (U Float CanonicalMixture) (hash-ref data2
                                                                   label2)])
                 (values label (cond [(and (flonum? val1) (flonum? val2))
                                      (fl* val1 val2)]
                                     [(and (list? val1) (flonum? val2))
                                      (mixed-product val2 val1)]
                                     [(and (list? val2) (flonum? val1))
                                      (mixed-product val1 val2)]
                                     [(and (list? val1) (list? val2))
                                      (canonical-mixture-product val1 val2)]
                                     [else 0.] ;unreachable, but makes type-checker happy.
                                     )))))
           (Factor joint-scope new-data null))
         (foldl single-product (car factors) (cdr factors))]))

(define (canonical-factor-marginalization [factor : CanonicalFactor]
                                          [var : GaussianRandomVar]) : CanonicalFactor
  (unless (set-member? (CanonicalFactor-scope factor) var)
    (error "Variable ~a is not in the scope of Factor ~a" var factor))
  (define reduced-scope (set-subtract (CanonicalFactor-scope factor) (set var)))
  (define reduced-order (order-of-elems (canonical-order reduced-scope)
                                        (canonical-order (CanonicalFactor-scope factor))))
  (define var-index (order-of-elems (list var)
                                    (canonical-order (CanonicalFactor-scope factor))))
  (define K (CanonicalFactor-K factor))
  (define K_reduced (submatrix K reduced-order reduced-order))
  (define K_reduced/var (submatrix K reduced-order var-index))
  (define K_var-inv (matrix-inverse (submatrix K var-index var-index)))
  (define K_var/reduced (submatrix K var-index reduced-order))
  (define h_reduced (~> (CanonicalFactor-h factor)
                        ->col-matrix
                        (submatrix _ reduced-order (list 0))
                        ->col-matrix))
  (define h_var (~> (CanonicalFactor-h factor)
                    ->col-matrix
                    (submatrix _ var-index (list 0))
                    ->col-matrix))
  (displayln K_var-inv)
  (displayln h_var)
  (CanonicalFactor reduced-scope
                   (matrix- K_reduced
                            (matrix* K_reduced/var K_var-inv K_var/reduced))
                   (matrix->vector (matrix- h_reduced
                                            (matrix* K_reduced/var K_var-inv h_var)))
                   (fl+ (CanonicalFactor-g factor)
                        (fl* 0.5
                             (fl+ (fllog (fl* (fl* 2. pi) (unpack K_var-inv)))
                                  (unpack (matrix* (matrix-transpose h_var)
                                                   K_var-inv
                                                   h_var)))))))
      
(define (factor-marginalization [factor : Factor]
                                [var : (U DiscreteRandomVar GaussianRandomVar)]) : Factor
  (unless (set-member? (Factor-scope factor) var)
    (error "Variable ~a is not in the scope of Factor ~a" var factor))
  (define reduced-scope (set-subtract (Factor-scope factor) (set var)))
  (cond [(GaussianRandomVar? var)
         (define discrete-scope : (Setof DiscreteRandomVar) (set-filter DiscreteRandomVar? reduced-scope))
         (Factor reduced-scope
                 (for/hash : FactorData ([label-list (make-TableCPD-labels (set->list discrete-scope))])
                   (let* ([label (list->set label-list)]
                          [datum (hash-ref (Factor-data factor) label)])
                     (values label
                             (if (list? datum)
                                 (map (λ ([pair : (Pairof Float CanonicalFactor)])
                                        (cons (car pair) (canonical-factor-marginalization (cdr pair) var)))
                                      datum)
                                 datum ; unreachable but the type-checker doesn't know this. :(
                                       ; since we are ensured by the unless clause that this
                                       ; continuous var is in the scope of the factor, we know
                                       ; datum is a CanonicalMixture and not a Float.
                                 ))))
                 null)]
        [(DiscreteRandomVar? var)
         factor]))