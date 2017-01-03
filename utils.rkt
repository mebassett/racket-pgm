#lang typed/racket

(require math/matrix)
(provide set-filter unpack)

(: set-filter (All (T R)
                   (case->
                    (-> (-> T Boolean) (Setof T) (Setof T))
                    (-> (-> T Boolean : R) (Setof T) (Setof R)))))
(define (set-filter valid? set)
  (list->set (filter valid? (set->list set))))

(define (unpack [x : (Matrix Float)]) : Float
  (matrix-ref x 0 0))