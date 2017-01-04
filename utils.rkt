#lang typed/racket

(require math/matrix)
(provide set-filter unpack order-of-elems)

(: set-filter (All (T R)
                   (case->
                    (-> (-> T Boolean) (Setof T) (Setof T))
                    (-> (-> T Boolean : R) (Setof T) (Setof R)))))
(define (set-filter valid? set)
  (list->set (filter valid? (set->list set))))

(define (unpack [x : (Matrix Float)]) : Float
  (matrix-ref x 0 0))

(: order-of-elems (All (A) (-> (Listof A)
                               (Listof A)
                               (Listof Integer))))
(define (order-of-elems small-list big-list)
  (define (helper [small-list : (Listof A)]
                  [big-list : (Listof A)]
                  [count : Integer]
                  [ret : (Listof Integer)]) : (Listof Integer)
    (cond [(or (null? big-list) (null? small-list)) ret]
          [(equal? (car small-list) (car big-list))
           (helper (cdr small-list) (cdr big-list) (add1 count) (append ret (list count)))]
          [else (helper small-list (cdr big-list) (add1 count) ret)]))
  (helper small-list big-list 0 '()))
    
                               