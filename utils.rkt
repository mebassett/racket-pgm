#lang typed/racket

(provide set-filter)

(: set-filter (All (T R)
                   (case->
                    (-> (-> T Boolean) (Setof T) (Setof T))
                    (-> (-> T Boolean : R) (Setof T) (Setof R)))))
(define (set-filter valid? set)
  (list->set (filter valid? (set->list set))))