#lang typed/racket
(require (prefix-in g: "graph-lib.rkt"))

(struct RandomVar ([name : Symbol]) #:transparent)
(struct DiscreteRandomVar RandomVar ([labels : (Listof Symbol)]) #:transparent)

(define-type DiscreteCPD (HashTable (Listof Symbol)
                                    Float))

(struct GraphicalModel ([g : (g:Graph DiscreteRandomVar)]
                        [cpds : (HashTable DiscreteRandomVar DiscreteCPD)]) #:transparent)

; (define (graphical-model-ready? [pgm : GraphicalModel]) : Boolean
;  (not 
;   (for/or ([v : RandomVar (g:get-vertices (g:transpose (GraphicalModel-g pgm)))])
;    (false? (hash-ref (GraphicalModel-cpds pgm) v #f)))))

;(define (create-naive-model [graph: (Graph RandomVar)]) : GraphicalModel

;)

(define (naive-table-pd [var : DiscreteRandomVar]) : DiscreteCPD
(for/hash : DiscreteCPD ([v (DiscreteRandomVar-labels var)])
 (values (list v)
         (exact->inexact (/ 1 (length (DiscreteRandomVar-labels var)))))))

(define (naive-table-cpd [var : DiscreteRandomVar]
                         [deps : (Listof DiscreteRandomVar)]) : DiscreteCPD
 (for/hash : DiscreteCPD ([k (apply cartesian-product
                                    (map DiscreteRandomVar-labels
                                         (cons var deps)))])
  (values k 
          (exact->inexact (/ 1 (length (DiscreteRandomVar-labels var)))))))

(define (make-naive-model [g : (g:Graph DiscreteRandomVar)]) : GraphicalModel
 (define g^T (g:transpose g))
 (define cpds : (HashTable DiscreteRandomVar DiscreteCPD) (make-hash))
 (define roots (g:get-leaves g^T))

 (define (travel-root [v : DiscreteRandomVar]) : Void
  (unless (hash-ref cpds v #f)
   (displayln v)
   (hash-set! cpds v (naive-table-cpd v (g:get-neighbors g^T v)))
   (unless (null? (g:get-neighbors g v))
    (for ([n (g:get-neighbors g v)])
     (travel-root n)))))

 (for ([root : DiscreteRandomVar roots])
  (hash-set! cpds root (naive-table-pd root))
  (for-each travel-root (g:get-neighbors g root)))

 (GraphicalModel g cpds))

(define (get-cpd [model : GraphicalModel]
                 [var : DiscreteRandomVar]) : DiscreteCPD
 (hash-ref (GraphicalModel-cpds model) var))
