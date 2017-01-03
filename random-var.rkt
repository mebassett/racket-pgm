#lang typed/racket
; not released under any license as this really isn't usable.
; will think about that if it gets that far.

(require math/flonum "utils.rkt")
(provide (struct-out RandomVar)
         (struct-out DiscreteRandomVar)
         (struct-out GaussianRandomVar)
         add-dependency!
         make-DiscreteRandomVar
         make-GaussianRandomVar
         get-var-labels)
         

(struct RandomVar ([name : Symbol]
                   [depends-on : (Listof RandomVar)]
                   [influences : (Listof RandomVar)])
  #:mutable)
(struct DiscreteRandomVar RandomVar ([labels : (Listof String)]))
(struct GaussianRandomVar RandomVar ())

(define (add-dependency! [var : RandomVar]
                         [dependency : RandomVar]) : Void
  (set-RandomVar-influences! dependency (cons var (RandomVar-influences dependency)))
  (set-RandomVar-depends-on! var (cons dependency (RandomVar-depends-on var))))

(define (make-DiscreteRandomVar [name : Symbol]
                                [labels : (Listof String)]) : DiscreteRandomVar
  (DiscreteRandomVar name '() '() labels))

(define (make-GaussianRandomVar [name : Symbol]) : GaussianRandomVar
  (GaussianRandomVar name '() '()))
  
(define (get-var-labels [vars : (Setof RandomVar)]) : (Setof Symbol)
  (list->set (set-map vars (λ ([v : RandomVar]) (RandomVar-name v)))))







                                   


                                   