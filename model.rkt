#lang typed/racket
(require "random-var.rkt"
         "factors.rkt"
         "utils.rkt"
         math/flonum)
(provide (struct-out Model)
         ;is-valid-model?
         make-naive-model
         update-cpd!
         probability-conditioned-on-evidence)
         

(struct Model ([vars : (HashTable Symbol RandomVar)]
               [cpds : (HashTable RandomVar Factor)]) #:mutable)

(define (make-naive-model [vars : (Listof (U GaussianRandomVar
                                             DiscreteRandomVar))]) : Model
  (define cpds : (HashTable RandomVar Factor) (make-hash))
  (for ([var vars])
    (hash-set! cpds var (cond [(DiscreteRandomVar? var)
                               (generate-uniform-table-cpd var)]
                              [(GaussianRandomVar? var)
                               (standard-gaussian-factor (set var))])))
  (Model (for/hash : (HashTable Symbol RandomVar) ([v vars])
           (values (RandomVar-name v) v))
         cpds))



(define (update-cpd! [model : Model]
                     [var : RandomVar]
                     [factor : Factor]) : Void
  (unless (is-valid-cpd? var factor)
    (error "Factor ~a is not valid for random variable ~a.  Please check its dependencies and that it sums to 1." factor var))
  (hash-set! (Model-cpds model) var factor))

(define (associated-factor [var : RandomVar]
                           [model : Model]) : Factor
  (hash-ref (Model-cpds model) var))
;
;(define (get-scope-symbols [scope : (Setof (Setof TableCPDIndex))]) : (Setof Symbol)
;    (list->set (set-map (first (set->list scope)) (λ ([x : TableCPDIndex]) (car x)))))
;
;(define (labels-within-scope [scope-symbols : (Setof Symbol)]
;                             [labels : (Setof TableCPDIndex)]) : (Setof TableCPDIndex)
;
;      (list->set
;       (filter (λ ([x : TableCPDIndex]) (set-member? scope-symbols (car x)))
;               (set->list labels))))
;
;
;
(define (probability-conditioned-on-evidence [model : Model]
                                             [vars : (Setof RandomVar)]
                                             [evidence : (Setof AnyIndex)]) : Factor
  (define factors (list->set
                   (map (λ ([var : RandomVar])
                          (partial-application-factor (associated-factor var model)
                                                      evidence))
                        (hash-values (Model-vars model)))))    
  (define evidence-labels (set-map evidence
                                   (λ ([x : AnyIndex]) (car x))))
  (define evidence-vars (list->set (map (λ ([x : Symbol])
                                          (hash-ref (Model-vars model) x))
                                        evidence-labels)))
  (define unnorm-factor
    (eliminate-variable factors
                        (set-subtract (list->set (hash-values (Model-vars model)))
                                      vars
                                      evidence-vars)))  
  (normalize-factor unnorm-factor))