#lang typed/racket
(require "factors.rkt"
         "random-var.rkt"
         "utils.rkt"
         threading
         math/flonum
         math/matrix
         math/array
         ; math/statistics
         ; psssh.  these are all Real and not Flonum.  
         "model.rkt")
(provide (all-defined-out))


(define-type Row (Vectorof (U Float String Null)))
(define-type Dataset (Listof Row))

;(define-type Incomplete-Row (Vectorof (U Float String Null)))
;(define-type Incomplete-Dataset (Listof Incomplete-Row))

(define-type Dataset-Header (HashTable Symbol Integer))

(define (make-dataset-header [dataset : Dataset]) : Dataset-Header
  (define var-header (vector-filter string? (car dataset)))
  (for/hash : Dataset-Header ([var-string-name var-header]
                              [index (in-range (vector-length var-header))])
    (values (string->symbol var-string-name)
            index)))     
(: get-var (case-> (-> Dataset-Header Row GaussianRandomVar (U Null Float))
                   (-> Dataset-Header Row DiscreteRandomVar (U Null String))
                   (-> Dataset-Header Row Symbol (U Null String))))
(define (get-var var-index-lookup row var)
  (let ([ret (vector-ref row (hash-ref var-index-lookup 
                                       (if (RandomVar? var)
                                           (RandomVar-name var)
                                           var)))])
    (cond [(or (DiscreteRandomVar? var) (symbol? var))
           (if (or (null? ret) (string? ret))
               ret
               (begin
                 (error "String expected at row ~a for var ~a." row var)
                 "null."))]
          [(GaussianRandomVar? var)
           (if (or (null? ret) (flonum? ret))
               ret
               (begin
                 (error "Flonum expected at row ~a for var ~a." row (RandomVar-name var))
                 "null."))])))

; takes a dataset and returns a new dataset partitioned by the set of values of its discrete 
; parents.  mostly used to find a sufficient statistic.
; it assumes that the first Vector is a list of variable names in the order that
; they appear in the dataset.
(define (filter-var-deps [var : RandomVar]
                         [dataset : Dataset]) : (HashTable (Setof TableCPDIndex) 
                                                           Dataset)
  (define (has-header? [dataset : Dataset]) : Boolean
    (andmap string? (vector->list (car dataset))))
  (define (right-shape? [dataset : Dataset]) : Boolean
    (define vlength (vector-length (car dataset)))
    (foldl (λ ([v1 : VectorTop] [last-status : Boolean])
             (and last-status (= (vector-length v1) vlength)))
           #t
           dataset))

  (unless (has-header? dataset)
    (error "Cannot understand dataset - first row should be vector of var names."))
  (unless (right-shape? dataset)
    (error "Dataset must have same number of columns throughout."))

  (define discrete-deps : (Listof DiscreteRandomVar)
    (filter DiscreteRandomVar? (RandomVar-depends-on var)))
  (define continuous-deps : (Listof GaussianRandomVar)
    (filter GaussianRandomVar? (RandomVar-depends-on var)))
  (define return-vars (cons var continuous-deps))
  (define var-index-lookup (make-dataset-header dataset))
  (define (filter-func? [labels : (Setof TableCPDIndex)]
                        [row : Row]) : Boolean
    (andmap (λ ([label : TableCPDIndex])
              (let ([val (get-var var-index-lookup row (car label))])
                (or (equal? (cdr label) val)
                    (null? val))))
            (set->list labels)))
  (define header (car dataset))
  (for/hash : (HashTable (Setof TableCPDIndex)
                         Dataset)
    ([labels-list (make-TableCPD-labels discrete-deps)])
    (define labels (list->set labels-list))

    (values labels
            (cons header
                  (filter (curry filter-func? labels)
                          (cdr dataset))))))

(define-type Gaussian-Statistic (HashTable (Setof TableCPDIndex)
                                           (Pairof (Vectorof Float)
                                                   (Matrix Float))))

; returns the mean and covariance matrix associated to a random var, partitioned by the
; values of its discrete parents.
; we return stats related to the variables in canonical order.

(define (create-gaussian-sufficient-statistic [var : GaussianRandomVar]
                                              [ds : Dataset]) : Gaussian-Statistic
  (define ds-header (make-dataset-header ds))
  (define continuous-var-names (map RandomVar-name
                                    (filter GaussianRandomVar?
                                            (cons var (RandomVar-depends-on var)))))
  (for/hash : Gaussian-Statistic ([(k dataset) (filter-var-deps var ds)])
    (define raw-header continuous-var-names)         
    (define ordered-header (~> raw-header list->set canonical-order))
    (define var-lookup (for*/hash : (HashTable Integer Integer)
                         ([i (in-range (length ordered-header))]
                          [(k v) ds-header]
                          #:when (equal? (list-ref ordered-header i) k))
                         (values i v)))

    (define (get-column [i : Integer]) : (Listof Float)
      (filter flonum?
              (map (λ ([row : Row]) (vector-ref row (hash-ref var-lookup i)))
                   (cdr dataset))))
    (define means
      (build-vector (length raw-header)
                    (λ ([i : Integer])
                      (flmean (get-column i)))))
    
    ; we should memoize and symmetrize this
    (define (matrix-elem [i : Integer]
                         [j : Integer]) : Float
      (flmean (array->list (array* (array- (list->array (get-column j))
                                           (array (vector-ref means j)))
                                   (array- (list->array (get-column i))
                                           (array (vector-ref means i)))))))
    
    (values k
            (cons means
                  (build-matrix (length raw-header)
                                (length raw-header)
                                matrix-elem)))))

(define-type Discrete-Statistic (HashTable (Setof TableCPDIndex)
                                           (Pairof Float Float)))
(define (create-discrete-sufficient-statistic [var : DiscreteRandomVar]
                                              [dataset : Dataset]) : Discrete-Statistic
  (define ds-header (make-dataset-header dataset))
  (define datasets (filter-var-deps var dataset))
  (define vals (list->set (DiscreteRandomVar-labels var)))
  (for*/hash : Discrete-Statistic
    ([(k v) datasets]
     [val vals])
    (define new-label : TableCPDIndex (cons (RandomVar-name var) val))
    (define total-length (->fl (length (cdr v))))
    (values (set-union k (set new-label))
            (cons (->fl (length (filter (λ ([row : Row])
                                          (equal? (get-var ds-header row var) val))
                                        (cdr v))))
                  total-length))))

(: estimate-discrete-sufficient-statistic (-> DiscreteRandomVar
                                              Dataset;Incomplete-Dataset
                                              Model
                                              (-> Model (Setof RandomVar) (Setof AnyIndex) Factor)
                                              Discrete-Statistic))
(define (estimate-discrete-sufficient-statistic var dataset model inferer)
  (define discrete-vars (list->set (cons var (filter DiscreteRandomVar?
                                                     (RandomVar-depends-on var)))))
  (define (row->any-index [excluding-vars : (Setof DiscreteRandomVar)]
                          [row : Row]) : (Setof AnyIndex)
    (for/set : (Setof AnyIndex) ([val row]
                                 [var-sym (~> dataset
                                              car
                                              vector->list
                                              (filter string? _)
                                              (map string->symbol _))]
                                 #:when (and (not (null? val))
                                             (not (set-member? excluding-vars
                                                               (hash-ref (Model-vars model)
                                                                         var-sym)))))
      ; (U (Pairof Symbol Float) (Pairof Symbol String)) != (Pairof Symbol (U String Float))
      ; in typed/racket, hence this craziness.
      (cond [(flonum? val)
             (cons var-sym val)]
            [(string? val)
             (cons var-sym val)]
            [else (cons var-sym 0.)] ;unreachable
            )))
  (define ds-header (make-dataset-header dataset))
  (define datasets (filter-var-deps var dataset))
  (define vals (list->set (DiscreteRandomVar-labels var)))
  (for*/hash : Discrete-Statistic 
    ([(k v) datasets]
     [val vals])
    (define new-label : TableCPDIndex (cons (RandomVar-name var) val))
    (define all-labels (set-union k (set new-label)))
    (define total-length (length (cdr v)))
    (define (row->statistic-var [row : Row]) : Float     
                         
      (cond [(equal? (get-var ds-header row var) val) 1.]
            [(null? (get-var ds-header row var))
             ;(ormap null? (set-map discrete-vars (λ ([v : DiscreteRandomVar])
             ;                                       (get-var ds-header row v))))
             (define p (inferer model
                                (list->set (cons var (filter DiscreteRandomVar?
                                                             (RandomVar-depends-on var))))
                                (row->any-index discrete-vars row)))             
             (p all-labels)]
            [else 0.]))
    
    (define (row->statistic-conditions [row : Row]) : Float     
                         
      (cond [(andmap (λ ([index : TableCPDIndex])
                       (equal? (get-var ds-header row (car index)) (cdr index)))
                     (set->list k))
                     1.]
            [(ormap null? (set-map k (λ ([index : TableCPDIndex])
                                       (get-var ds-header row (car index)))))
             (define p (inferer model
                                (list->set (filter DiscreteRandomVar?
                                                   (RandomVar-depends-on var)))
                                (row->any-index (set-subtract discrete-vars (set var))
                                                row))) 
             (p k)]
            [else 0.]))
    
    (values all-labels
            (cons (flsum (map row->statistic-var (cdr v)))
                  (flsum (map row->statistic-conditions (cdr v)))))))

(: discrete-statistic->factor (-> DiscreteRandomVar Discrete-Statistic Factor))
(define (discrete-statistic->factor var statistic) 
  (make-factor-from-table var
                          (for/hash : (HashTable (Setof TableCPDIndex) Float)
                            ([(k v) statistic])
                            (values k (fl/ (car v)
                                           (cdr v))))))

;(: create-factor-from-data (case-> (-> DiscreteRandomVar Dataset Factor)
;                                   (-> GaussianRandomVar Dataset Factor)))
;
; another typed/racket inference complaint here.
; > (: my-test (case-> (-> String Integer)
;                      (-> Symbol Integer)))
; > (define (my-test string/sym)
;     1)
; > (define what-am-I? : (U String Symbol) 'a)
; > (my-test what-am-I?)
; Type Checker: No function domains matched in function application:
; Domains: Symbol
;          String
; Arguments: (U Symbol String)
;  in: (my-test what-am-I?)
;
; I'm not too experienced with stongly typed systems, but I think that case is pretty clear to handle.
; have not raised it to the racket guys, but hopefully in 6.8 it'll go away.  In the meantime...

(: create-factor-from-data (-> RandomVar Dataset Factor))
(define (create-factor-from-data var dataset)
  (cond [(DiscreteRandomVar? var)
         (define statistic (create-discrete-sufficient-statistic var dataset))
         (discrete-statistic->factor var statistic)]
        [(GaussianRandomVar? var)
         (define statistic (create-gaussian-sufficient-statistic var dataset))
         (define full-scope (list->set (cons var (RandomVar-depends-on var))))
         (define deps-scope (list->set (filter GaussianRandomVar?
                                               (RandomVar-depends-on var))))
         (define continuous-scope (set-add deps-scope var))
         (define deps-index (order-of-elems (canonical-order deps-scope)
                                            (canonical-order continuous-scope)))
         (define data
           (for/hash : FactorData ([(k v) statistic])
         
             (define μ (car v))
             (define Σ (cdr v))
             (if (set-empty? deps-scope)
                 (values k (list (cons 1. (make-gaussian continuous-scope μ Σ))))
                 (let* ([Sigma_x (submatrix Σ deps-index deps-index)]
                        [mu_x (~> μ
                                  ->col-matrix
                                  (submatrix _ deps-index 0)
                                  array->vector)]
                          
             
                        [cf_XY (make-gaussian continuous-scope μ Σ)]
                        [cf_X (make-gaussian deps-scope mu_x Sigma_x)])
                   (values k (list (cons 1. (CanonicalFactor continuous-scope
                                                             (K- cf_XY cf_X)
                                                             (h- cf_XY cf_X)
                                                             (fl- (CanonicalFactor-g cf_XY)
                                                                  (CanonicalFactor-g cf_X))))))))))
         (Factor full-scope
                 data
                 null)]
        [else ;unreachable.  why am I doing this? see complaint.
         (Factor (set var)
                 (ann (hash) FactorData)
                 null)]))
         
; more typed racket complaints
; let's say we want to filter out null elements.  Let's start with a "not-null?"!
; (: not-null? (-> Any Boolean : #:- Null))
; (define not-null? (compose not null?))
; oops!  that doesn't work, it doesn't match our negative predicate.  ugh.
; this does though:
; (define (not-null? x) (not (null? x)))
; A tad bit annoying, but no biggie.  Now let's filter:
; (define l : (Listof (U String Null))
;    (list "adsf" null "bad"))
; (filter not-null? l)
; - : (Listof (U Null String))
; of course! because if you check the type signature for filter, it says *nothing*
; about negative predicates.
;
; so in a perfect would, we would do
; (filter (compose not null?) list-row)
; here, we have to do
; (: string-or-float? (-> (U String Float Null) Boolean : (U String Float))
; (of course, (define string-or-float? (compose not null?)) isn't going to work, it gets
; stupider
; (define (string-or-float? x) (or (string? x) (flonum? x)))
; (filter string-or-float? list-row)



(define (mle-train-model! [model : Model]
                          [dataset : Dataset]) : Void
  (for ([v (hash-values (Model-vars model))])
    (update-cpd! model v (create-factor-from-data v dataset))))
  
