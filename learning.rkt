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


(define-type Row (Vectorof (U Float String)))
(define-type Dataset (Listof Row))


; takes a dataset and returns a new dataset containing only this var along with 
; its continuous parents partitioned by the set of values of its discrete 
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
  
  (define (make-var-header [dataset : Dataset]) : (HashTable Symbol Integer)
    (define var-header (vector-filter string? (car dataset)))
    (for/hash : (HashTable Symbol Integer) ([var-string-name var-header]
                                            [index (in-range (vector-length var-header))])
      (values (string->symbol var-string-name)
              index)))     


  (unless (has-header? dataset)
    (error "Cannot understand dataset - first row should be vector of var names."))
  (unless (right-shape? dataset)
    (error "Dataset must have same number of columns throughout."))


  (define discrete-deps : (Listof DiscreteRandomVar)
    (filter DiscreteRandomVar? (RandomVar-depends-on var)))
  (define continuous-deps : (Listof GaussianRandomVar)
    (filter GaussianRandomVar? (RandomVar-depends-on var)))
  (define return-vars (cons var continuous-deps))
  (define var-index-lookup (make-var-header dataset))

  (: get-var (case-> (-> Row GaussianRandomVar Float)
                     (-> Row DiscreteRandomVar String)
                     (-> Row Symbol String)))
  (define (get-var row var)
    (let ([ret (vector-ref row (hash-ref var-index-lookup 
                                         (if (RandomVar? var)
                                             (RandomVar-name var)
                                             var)))])
      (cond [(or (DiscreteRandomVar? var) (symbol? var))
             (if (string? ret)
                 ret
                 (begin
                   (error "String expected at row ~a for var ~a." row var)
                   "null."))]
            [(GaussianRandomVar? var)
             (if (flonum? ret)
                 ret
                 (begin
                   (error "Flonum expected at row ~a for var ~a." row (RandomVar-name var))
                   "null."))])))
  (define (filter-func? [labels : (Setof TableCPDIndex)]
                        [row : Row]) : Boolean
    (andmap (λ ([label : TableCPDIndex]) 
              (equal? (cdr label) (get-var row (car label))))
            (set->list labels)))

  (define (restrict-to-vars [row : Row]) : Row
    (for/vector : Row ([v return-vars])
      (cond [(DiscreteRandomVar? v) (get-var row v)]
            [(GaussianRandomVar? v) (get-var row v)]
            [else 0.0] ; unreachable
            )))
  (define header
    (for/vector : Row ([v return-vars])
      (symbol->string (RandomVar-name v))))
  (for/hash : (HashTable (Setof TableCPDIndex)
                         Dataset)
    ([labels-list (make-TableCPD-labels discrete-deps)])
    (define labels (list->set labels-list))
    (values labels
            (cons header
                  (map restrict-to-vars
                       (filter (curry filter-func? labels)
                               (cdr dataset)))))))

(define-type Gaussian-Statistic (HashTable (Setof TableCPDIndex)
                                           (Pairof (Vectorof Float)
                                                   (Matrix Float))))

; returns the mean and covariance matrix associated to a random var, partitioned by the
; values of its discrete parents.
; we return stats related to the variables in canonical order.
; (one could leave the var in fron and put its deps in canonical order.  that might make more sense.
;  but I feel its a lot of work either way.)
(define (create-gaussian-sufficient-statistic [var : GaussianRandomVar]
                                              [ds : Dataset]) : Gaussian-Statistic
  (for/hash : Gaussian-Statistic ([(k dataset) (filter-var-deps var ds)])
    (define raw-header (~> (car dataset)
                           vector->list
                           (filter string? _)))
         
    (define ordered-header (~> raw-header
                               (map string->symbol _)
                               list->set canonical-order))
    (define var-lookup (for*/hash : (HashTable Integer Integer)
                         ([i (in-range (length ordered-header))]
                          [j (in-range (length raw-header))]
                          #:when (equal? (list-ref ordered-header i)
                                         (~> j (list-ref raw-header _) string->symbol)))
                         (values j i)))

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

(define (create-discrete-sufficient-statistic [var : DiscreteRandomVar]
                                              [dataset : Dataset]) : (HashTable (Setof TableCPDIndex)
                                                                                (Pairof Integer Integer))
  (define datasets (filter-var-deps var dataset))
  (define all-vals (apply append (map (λ ([ds : Dataset])
                                        (filter string?
                                                (map (λ ([row : Row]) (vector-ref row 0))
                                                     (cdr ds))))
                                      (hash-values datasets))))
  (define vals (list->set all-vals))
  (for*/hash : (HashTable (Setof TableCPDIndex)
                          (Pairof Integer Integer)) 
    ([(k v) datasets]
     [val vals])
    (define new-label : TableCPDIndex (cons (RandomVar-name var) val))
    (define total-length (length (cdr v)))
    (values (set-union k (set new-label))
            (cons (length (filter (λ ([row : Row])
                                    (equal? (vector-ref row 0) val))
                                  v))
                  total-length))))



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
  
         (make-factor-from-table var
                                 (for/hash : (HashTable (Setof TableCPDIndex) Float)
                                   ([(k v) statistic])
                                   (values k (fl/ (->fl (car v))
                                                  (->fl (cdr v))))))]
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
             (displayln (RandomVar-name var))
             (displayln (map RandomVar-name (RandomVar-depends-on var)))
             (displayln deps-scope)
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
         
         

(define (mle-train-model! [model : Model]
                          [dataset : Dataset]) : Void
  (for ([v (hash-values (Model-vars model))])
    (update-cpd! model v (create-factor-from-data v dataset))))
  
