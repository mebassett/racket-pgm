#lang typed/racket

(require math/matrix math/flonum)
(provide set-filter unpack order-of-elems flmean)

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

; the following are from math/statistics, but the types are changed to Flonum

(: check-lengths! (All (A B) (Symbol String A B Index Index -> Void)))
(define (check-lengths! name what xs ys m n)
  (unless (= m n) (error name "~a must be the same length; given ~e (length ~a) and ~e (length ~a)"
                         what xs m ys n)))

(: sequences->weighted-samples
   (All (A) (Symbol (Sequenceof A) (Sequenceof Float)
                    -> (Values (Listof A) (Listof Float)))))
(define (sequences->weighted-samples name x-seq w-seq)
  (define xs (sequence->list x-seq))
  (define ws
    (for/list: : (Listof Float) ([w w-seq])
      (cond [(w . >= . 0)  w]
            [else  (raise-argument-error name "(Sequenceof Nonnegative-Real)" 1 x-seq w-seq)])))
  (check-lengths! name "values and weights" xs ws (length xs) (length ws))
(values xs ws))

(: flmean (case-> ((Sequenceof Float) -> Float)
                ((Sequenceof Float) (Option (Sequenceof Float)) -> Float)))
(define (flmean xs [ws #f])
  (cond [ws  (let-values ([(xs ws)  (sequences->weighted-samples 'mean xs ws)])
               (define n (flsum ws))
               (cond [(zero? n)  +nan.0]
                     [else  (/ (flsum (map * xs ws)) n)]))]
        [else  (let ([xs  (sequence->list xs)])
                 (define n (length xs))
                 (cond [(zero? n)  +nan.0]
                       [else (/ (flsum xs) n)]))]))                               