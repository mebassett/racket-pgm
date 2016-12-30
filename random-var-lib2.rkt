#lang typed/racket

(struct RandomVar ([name : Symbol]
                   [depends-on : (Listof RandomVar)]
                   [influences : (Listof RandomVar)])
  #:mutable)
(struct DiscreteRandomVar RandomVar ([labels : (Listof Symbol)]))

(define (add-dependency! [var : RandomVar]
                         [dependency : RandomVar]) : Void
  (set-RandomVar-influences! dependency (cons var (RandomVar-influences dependency)))
  (set-RandomVar-depends-on! var (cons dependency (RandomVar-depends-on var))))

(define (make-DiscreteRandomVar [name : Symbol]
                                [labels : (Listof Symbol)]) : DiscreteRandomVar
  (DiscreteRandomVar name '() '() labels))

(define difficulty (make-DiscreteRandomVar 'difficulty '(hard easy)))
(define grade (make-DiscreteRandomVar 'grade '(a b f)))
(define SAT (make-DiscreteRandomVar 'SAT '(high low)))
(define letter (make-DiscreteRandomVar 'letter '(good poor)))
(define intelligence (make-DiscreteRandomVar 'iq '(high low)))

(add-dependency! grade difficulty)
(add-dependency! grade intelligence)
(add-dependency! SAT intelligence)
(add-dependency! letter grade)