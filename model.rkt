#lang typed/racket
(require "random-var.rkt"
         "factors.rkt"
         "utils.rkt"
         math/flonum)
(provide (struct-out DiscreteModel)
         is-valid-model?
         make-naive-model
         update-cpd!
         probability-conditioned-on-evidence)

(struct DiscreteModel ([vars : (HashTable Symbol DiscreteRandomVar)]
                       [cpds : (HashTable DiscreteRandomVar TableCPD)]) #:mutable)

(define (get-variable-labels [var : DiscreteRandomVar]) : (Setof (Setof TableCPDIndex))
  (for/set : (Setof (Setof TableCPDIndex)) ([i (DiscreteRandomVar-labels var)])
    (set (cons (RandomVar-name var) i))))

(define (get-variable-labels* [var : DiscreteRandomVar]) : (Setof TableCPDIndex)
    (foldl (λ ([x : (Setof TableCPDIndex)]
               [y : (Setof TableCPDIndex)]) : (Setof TableCPDIndex)
             (set-union x y))
           (ann (set) (Setof TableCPDIndex))
           (set->list (get-variable-labels var))))

(define (is-valid-model? [model : DiscreteModel]) : Boolean
  (define (all-cpds-exist?)
    (for/and : Boolean ([v (hash-values (DiscreteModel-vars model))]) (not (false? (hash-ref (DiscreteModel-cpds model) v #f)))))
  (define (all-cpds-valid?)
    (for/and : Boolean ([(var cpd) (DiscreteModel-cpds model)]) (is-valid-cpd? var cpd)))
  (and (all-cpds-exist?)
       (all-cpds-valid?)))

(define (make-naive-model [vars : (Listof DiscreteRandomVar)]) : DiscreteModel
  (define cpds : (HashTable DiscreteRandomVar TableCPD) (make-hash))
  (for ([var vars])
    (hash-set! cpds var (generate-uniform-table-cpd var)))
  (DiscreteModel (for/hash : (HashTable Symbol DiscreteRandomVar) ([v vars])
                   (values (RandomVar-name v) v))
                 cpds))

(define (update-cpd! [model : DiscreteModel]
                     [var : DiscreteRandomVar]
                     [pre-cpd : (HashTable (Listof TableCPDIndex) Float)]) : Void
  (define cpd
    (for/hash : TableCPD ([([k : (Listof TableCPDIndex)]
                            [v : Float])          pre-cpd])
      (values (list->set k) v)))
  (unless (member var (hash-values (DiscreteModel-vars model)))
    (error "Random Variable ~a not in Model ~a" var model))
  (unless (is-valid-cpd? var cpd)
    (error "TableCPD ~a is not valid for random variable ~a.  Please check its dependencies and that it sums to 1." cpd var))
  (hash-set! (DiscreteModel-cpds model) var cpd))

; as from Friedman & Koller, Probabilistic Graphical Models Alg9.1
(define (sum-product-eliminate-var [factors : (Setof DiscreteFactor)]
                                   [var : DiscreteRandomVar]
                                   [model : DiscreteModel]) : (Setof DiscreteFactor)
  (define factors-with-var (list->set (filter (λ ([factor : DiscreteFactor]) (var-in-scope? var factor))
                                              (set->list factors))))
  (define factors-without-var (set-subtract factors factors-with-var))
  (define new-factor (factor-marginalization (apply (curry product-factor model) (set->list factors-with-var)) var model))
  (set-add factors-without-var new-factor))

; ibid
; we need to define a sensible order here.  set->list gives a different order each time and not all these orders make sense.

(define (eliminate-variable [factors : (Setof DiscreteFactor)]
                          [vars : (Setof DiscreteRandomVar)]
                          [model : DiscreteModel]) : DiscreteFactor
  (apply (curry product-factor model) (set->list (foldl (λ ([var : DiscreteRandomVar]
                                              [fs : (Setof DiscreteFactor)]) : (Setof DiscreteFactor)
                                            (sum-product-eliminate-var fs var model))
                                          factors
                                          (set->list vars)))))
(: partial-application-factor (-> DiscreteFactor (Setof TableCPDIndex) DiscreteFactor))
  (define (partial-application-factor factor evidence)
    (define evidence-labels (set-map evidence (λ ([x : TableCPDIndex]) (car x))))
    (define reduced-scope (for/set : (Setof (Setof TableCPDIndex)) ([labels (DiscreteFactor-scope factor)])
                            (list->set (filter (λ ([index : TableCPDIndex]) (false? (member (car index) evidence-labels))) (set->list labels)))))
    (define reduced-evidence (set-filter (λ ([x : TableCPDIndex]) (set-member? (get-scope-symbols (DiscreteFactor-scope factor)) (car x)))
                                         evidence))
      
    (DiscreteFactor reduced-scope (λ (labels) (factor (set-union labels reduced-evidence)))))


(define (var-in-scope? [var : DiscreteRandomVar] [factor : DiscreteFactor]) : Boolean
  ; we should probably change factors to have a list of their random vars, rather than the entire domain.
  (set-member? (get-scope-symbols (DiscreteFactor-scope factor)) (RandomVar-name var)))
  ;(for/and : Boolean
  ;  ([i (DiscreteFactor-scope factor)]
  ;   [j (get-variable-labels var)])
  ;  (subset? j i)))

(define (factor-marginalization [factor : DiscreteFactor]
                                [var : DiscreteRandomVar]
                                [model : DiscreteModel]) : DiscreteFactor
  ;(define var-scope (list->set (hash-keys (hash-ref (DiscreteModel-cpds model) var)))) 
  (unless (var-in-scope? var factor)
    (error "Variable ~a is not in the scope of Factor ~a" var factor))
  (define variable-labels (get-variable-labels* var))
 
  (DiscreteFactor (for/set : (Setof (Setof TableCPDIndex))
                    ([i (DiscreteFactor-scope factor)])
                    (set-filter (λ ([x : TableCPDIndex])
                                  (false? (set-member? variable-labels x)))
                                i))

                  (λ (factor-labels)
                    (flsum (map (λ ([var-label : (Setof TableCPDIndex)]) : Float
                                  (factor (set-union factor-labels var-label)))
                                (set->list (get-variable-labels var)))))))

(define (get-scope-symbols [scope : (Setof (Setof TableCPDIndex))]) : (Setof Symbol)
    (list->set (set-map (first (set->list scope)) (λ ([x : TableCPDIndex]) (car x)))))

(define (labels-within-scope [scope-symbols : (Setof Symbol)]
                             [labels : (Setof TableCPDIndex)]) : (Setof TableCPDIndex)

      (list->set
       (filter (λ ([x : TableCPDIndex]) (set-member? scope-symbols (car x)))
               (set->list labels))))

(define (associated-factor [var : DiscreteRandomVar]
                             [model : DiscreteModel]) : DiscreteFactor
    (define cpd (hash-ref (DiscreteModel-cpds model) var))
    (DiscreteFactor (list->set (hash-keys cpd)) (λ (labels) (hash-ref cpd labels))))

(: product-factor (-> DiscreteModel DiscreteFactor * DiscreteFactor))
(define (product-factor model . factors)
  (define (helper [f1 : DiscreteFactor]
                  [f2 : DiscreteFactor]) : DiscreteFactor
    (define f1-symbols (get-scope-symbols (DiscreteFactor-scope f1)))
    (define f2-symbols (get-scope-symbols (DiscreteFactor-scope f2)))
    (define joint-scope 
      (list->set (map (λ ([labels : (Listof TableCPDIndex)]) (list->set labels))
                      (apply cartesian-product
                             (set-map (set-union f1-symbols f2-symbols)
                                      (λ ([x : Symbol])
                                        (set->list (get-variable-labels* (hash-ref (DiscreteModel-vars model) x)))))))))
    (define (factor-func [labels : (Setof TableCPDIndex)]) : Float
      (* (f1 (labels-within-scope f1-symbols labels))
         (f2 (labels-within-scope f2-symbols labels))))
    (DiscreteFactor joint-scope factor-func))

  (foldl helper (car factors) (cdr factors)))


(define (probability-conditioned-on-evidence [model : DiscreteModel]
                                             [vars : (Setof DiscreteRandomVar)]
                                             [evidence : (Setof TableCPDIndex)]) : DiscreteFactor
    (define factors (list->set (map (λ ([var : DiscreteRandomVar]) (partial-application-factor (associated-factor var model) evidence))
                                    (hash-values (DiscreteModel-vars model)))))
    
    (define evidence-labels (set-map evidence (λ ([x : TableCPDIndex]) (car x))))
    (define evidence-vars (list->set (map (λ ([x : Symbol]) (hash-ref (DiscreteModel-vars model) x)) evidence-labels)))
    (define unnorm-factor (eliminate-variable factors
                                              (set-subtract (list->set (hash-values (DiscreteModel-vars model))) vars evidence-vars)
                                              model))
    (define alpha (flsum (set-map (DiscreteFactor-scope unnorm-factor) (λ ([x : (Setof TableCPDIndex)]) (unnorm-factor x)))))
    (DiscreteFactor (DiscreteFactor-scope unnorm-factor)
                    (λ (label)
                      (/ (unnorm-factor label) alpha))))
